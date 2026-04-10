#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x
r_string_literal <- function(x) encodeString(x, quote = "\"")

timestamp_tag <- function(x = Sys.time()) {
  format(as.POSIXct(x, tz = Sys.timezone()), "%Y%m%d-%H%M%S")
}

parse_args <- function(args) {
  out <- list(
    safe_main = "--safe-main" %in% args,
    output_dir = file.path("tools", "benchmark-results", paste0("mlx-native-rsvd-", timestamp_tag()))
  )

  for (arg in args) {
    if (!startsWith(arg, "--")) next
    key <- sub("^--", "", arg)
    if (!grepl("=", key, fixed = TRUE)) next
    pieces <- strsplit(key, "=", fixed = TRUE)[[1L]]
    name <- pieces[[1L]]
    value <- paste(pieces[-1L], collapse = "=")
    if (identical(name, "output-dir")) out$output_dir <- value
  }

  out
}

relaunch_safe_master_if_needed <- function(args) {
  raw_args <- commandArgs(trailingOnly = FALSE)
  direct_file_arg <- grep("^--file=", raw_args, value = TRUE)
  direct_file_entry <- length(direct_file_arg) > 0L

  if (!direct_file_entry || isTRUE(args$safe_main)) {
    return(invisible(FALSE))
  }

  script_path <- normalizePath(sub("^--file=", "", direct_file_arg[[1L]]), winslash = "/", mustWork = TRUE)
  repo_root <- normalizePath(dirname(dirname(script_path)), winslash = "/", mustWork = TRUE)
  expr <- sprintf(
    "setwd(%s); source(%s, local = globalenv())",
    r_string_literal(repo_root),
    r_string_literal(script_path)
  )
  relaunch_args <- c("-e", expr, "--args", "--safe-main", commandArgs(trailingOnly = TRUE))
  quoted_args <- vapply(relaunch_args, shQuote, character(1), USE.NAMES = FALSE)
  status <- system2(file.path(R.home("bin"), "Rscript"), quoted_args)
  quit(save = "no", status = status)
}

initialize_context <- function() {
  if (file.exists(file.path("tools", "benchmark-helpers.R"))) {
    source(file.path("tools", "benchmark-helpers.R"), local = FALSE)
  }
  load_benchmark_amatrix()
  if (!requireNamespace("irlba", quietly = TRUE)) {
    stop("Package 'irlba' is required for this benchmark", call. = FALSE)
  }

  specs <- .benchmark_optional_backend_specs(include_arrayfire = FALSE)
  if (!isTRUE(.benchmark_enable_backend(specs[["mlx"]]))) {
    stop("MLX backend is not available", call. = FALSE)
  }

  options(
    amatrix.mlx.available = TRUE,
    amatrix.mlx.safe_spectral = FALSE,
    amatrix.mlx.rsvd.engine = "resident"
  )
  Sys.unsetenv("AMATRIX_MLX_SAFE_SPECTRAL")
  Sys.setenv(AMATRIX_MLX_NATIVE_SPECTRAL = "1")
  invisible(TRUE)
}

default_cases <- function() {
  list(
    list(id = "500x400-k20", n = 500L, p = 400L, k = 20L, n_oversamples = 10L, n_iter = 2L),
    list(id = "1000x800-k20", n = 1000L, p = 800L, k = 20L, n_oversamples = 10L, n_iter = 2L),
    list(id = "3000x1200-k40", n = 3000L, p = 1200L, k = 40L, n_oversamples = 12L, n_iter = 2L)
  )
}

make_host_case <- function(case) {
  set.seed(20260409L + case$n + case$p + case$k)
  matrix(rnorm(case$n * case$p), nrow = case$n, ncol = case$p)
}

benchmark_elapsed <- function(fn, reps = 3L, warmup = NULL) {
  if (is.function(warmup)) warmup()
  timings <- numeric(reps)
  last <- NULL
  for (idx in seq_len(reps)) {
    gc()
    timings[[idx]] <- system.time(last <- fn())[["elapsed"]]
  }
  list(elapsed = median(timings), result = last)
}

relative_sv_error <- function(actual, expected) {
  max(abs(actual - expected) / pmax(abs(expected), 1e-12))
}

new_row <- function(...) {
  defaults <- list(
    case = NA_character_,
    algorithm = "rsvd",
    backend = NA_character_,
    precision = NA_character_,
    status = "ok",
    reason = NA_character_,
    elapsed = NA_real_,
    rel_sv_err = NA_real_,
    iter = NA_integer_,
    mprod = NA_integer_,
    selected_backend = NA_character_,
    dispatch_state = NA_character_,
    stringsAsFactors = FALSE
  )
  as.data.frame(modifyList(defaults, list(...)), stringsAsFactors = FALSE)
}

benchmark_case <- function(case, reps = 3L) {
  host <- make_host_case(case)
  ref_sv <- base::svd(host, nu = case$k, nv = case$k)$d[seq_len(case$k)]

  cpu <- benchmark_elapsed(
    function() irlba::svdr(host, k = case$k, extra = case$n_oversamples, it = case$n_iter),
    reps = reps
  )

  x <- adgeMatrix(host, preferred_backend = "mlx", precision = "fast")
  mlx <- benchmark_elapsed(
    function() rsvd(x, k = case$k, n_oversamples = case$n_oversamples, n_iter = case$n_iter),
    reps = reps,
    warmup = function() invisible(rsvd(x, k = case$k, n_oversamples = case$n_oversamples, n_iter = case$n_iter))
  )

  rbind(
    new_row(
      case = case$id,
      backend = "cpu",
      precision = "strict",
      elapsed = cpu$elapsed,
      rel_sv_err = relative_sv_error(cpu$result$d[seq_len(case$k)], ref_sv),
      iter = as.integer(cpu$result$iter %||% NA_integer_),
      mprod = as.integer(cpu$result$mprod %||% NA_integer_),
      selected_backend = "cpu",
      dispatch_state = "cpu_baseline"
    ),
    new_row(
      case = case$id,
      backend = "mlx",
      precision = "fast",
      elapsed = mlx$elapsed,
      rel_sv_err = relative_sv_error(mlx$result$d[seq_len(case$k)], ref_sv),
      iter = as.integer(mlx$result$iter %||% NA_integer_),
      mprod = as.integer(mlx$result$mprod %||% NA_integer_),
      selected_backend = "mlx",
      dispatch_state = "accelerated",
      reason = "native top-level resident rsvd"
    )
  )
}

write_outputs <- function(rows, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cpu <- rows[rows$backend == "cpu", c("case", "elapsed")]
  names(cpu)[names(cpu) == "elapsed"] <- "cpu_reference_elapsed"
  summary <- merge(rows, cpu, by = "case", all.x = TRUE)
  summary$speedup_vs_cpu <- ifelse(
    is.na(summary$cpu_reference_elapsed) | is.na(summary$elapsed),
    NA_real_,
    summary$cpu_reference_elapsed / summary$elapsed
  )
  summary <- summary[order(summary$case, summary$backend), ]

  raw_path <- file.path(output_dir, "raw-results.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  metadata_path <- file.path(output_dir, "metadata.rds")
  write.csv(rows, raw_path, row.names = FALSE)
  write.csv(summary, summary_path, row.names = FALSE)
  saveRDS(
    list(
      created_at = Sys.time(),
      hostname = Sys.info()[["nodename"]],
      r_version = paste(R.version$major, R.version$minor, sep = "."),
      platform = R.version$platform,
      note = "Run as a top-level Rscript -e source(...) process to avoid MLX subprocess initialization crashes."
    ),
    metadata_path
  )

  list(raw = raw_path, summary = summary_path, metadata = metadata_path)
}

run_benchmark <- function(output_dir) {
  initialize_context()
  rows <- do.call(rbind, lapply(default_cases(), benchmark_case, reps = 3L))
  row.names(rows) <- NULL
  paths <- write_outputs(rows, output_dir)

  print(read.csv(paths$summary, stringsAsFactors = FALSE), row.names = FALSE)
  cat("\nArtifacts:\n")
  cat("  raw: ", paths$raw, "\n", sep = "")
  cat("  summary: ", paths$summary, "\n", sep = "")
  cat("  metadata: ", paths$metadata, "\n", sep = "")
  invisible(paths)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
relaunch_safe_master_if_needed(args)

if (!identical(Sys.getenv("AMATRIX_BENCHMARK_NO_AUTORUN", unset = ""), "1")) {
  run_benchmark(args$output_dir)
}
