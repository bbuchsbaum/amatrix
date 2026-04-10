#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x
r_string_literal <- function(x) encodeString(x, quote = "\"")

parse_args <- function(args) {
  out <- list(
    worker = "--worker" %in% args,
    mlx_native_spectral = "--mlx-native-spectral" %in% args ||
      identical(Sys.getenv("AMATRIX_MLX_NATIVE_SPECTRAL", unset = ""), "1"),
    mlx_native_inline = "--mlx-native-inline" %in% args ||
      identical(Sys.getenv("AMATRIX_MLX_NATIVE_INLINE", unset = ""), "1"),
    plan = NULL,
    out = NULL,
    entry_id = NULL,
    output_dir = file.path("tools", "benchmark-results", paste0("svd-backends-", format(Sys.time(), "%Y%m%d-%H%M%S")))
  )

  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }
    key <- sub("^--", "", arg)
    if (!grepl("=", key, fixed = TRUE)) {
      next
    }
    pieces <- strsplit(key, "=", fixed = TRUE)[[1L]]
    name <- pieces[[1L]]
    value <- paste(pieces[-1L], collapse = "=")
    if (identical(name, "plan")) out$plan <- value
    if (identical(name, "out")) out$out <- value
    if (identical(name, "entry-id")) out$entry_id <- value
    if (identical(name, "output-dir")) out$output_dir <- value
  }

  out
}

available_svd_backends <- function() {
  available_benchmark_backends(
    include_cpu = TRUE,
    include_mlx = TRUE,
    include_metal = FALSE,
    include_opencl = TRUE,
    include_arrayfire = .benchmark_arrayfire_requested()
  )
}

initialize_svd_benchmark_context <- local({
  initialized <- FALSE

  function() {
    if (initialized) {
      return(invisible(TRUE))
    }

    suppressPackageStartupMessages({
      if (file.exists(file.path("tools", "benchmark-helpers.R"))) {
        source(file.path("tools", "benchmark-helpers.R"), local = FALSE)
      }
      load_benchmark_amatrix()
    })

    if (!requireNamespace("irlba", quietly = TRUE)) {
      stop("Package 'irlba' is required for this benchmark", call. = FALSE)
    }

    initialized <<- TRUE
    invisible(TRUE)
  }
})

activate_worker_backend <- function(backend_name) {
  if (identical(backend_name, "cpu")) {
    return(TRUE)
  }

  specs <- .benchmark_optional_backend_specs(include_arrayfire = .benchmark_arrayfire_requested())
  spec <- specs[[backend_name]]
  if (is.null(spec)) {
    return(FALSE)
  }

  isTRUE(.benchmark_enable_backend(spec))
}

default_cases <- function() {
  list(
    list(id = "500x400-k20", n = 500L, p = 400L, k = 20L, n_oversamples = 10L, n_iter = 2L),
    list(id = "1000x800-k20", n = 1000L, p = 800L, k = 20L, n_oversamples = 10L, n_iter = 2L),
    list(id = "3000x1200-k40", n = 3000L, p = 1200L, k = 40L, n_oversamples = 12L, n_iter = 2L)
  )
}

benchmark_elapsed <- function(fn, reps = 3L, warmup = NULL) {
  if (is.function(warmup)) {
    warmup()
  }

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
    algorithm = NA_character_,
    backend = NA_character_,
    precision = NA_character_,
    status = "ok",
    reason = NA_character_,
    elapsed = NA_real_,
    rel_sv_err = NA_real_,
    iter = NA_integer_,
    mprod = NA_integer_,
    selected_backend = NA_character_,
    stringsAsFactors = FALSE
  )
  as.data.frame(modifyList(defaults, list(...)), stringsAsFactors = FALSE)
}

make_host_case <- function(case) {
  set.seed(20260409L + case$n + case$p + case$k)
  matrix(rnorm(case$n * case$p), nrow = case$n, ncol = case$p)
}

benchmark_reference_case <- function(case, reps = 3L) {
  host <- make_host_case(case)
  ref_exact <- benchmark_elapsed(function() base::svd(host, nu = case$k, nv = case$k), reps = 1L)
  ref_sv <- ref_exact$result$d[seq_len(case$k)]

  ref_rsvd <- benchmark_elapsed(
    function() irlba::svdr(host, k = case$k, extra = case$n_oversamples, it = case$n_iter),
    reps = reps
  )

  rbind(
    new_row(
      case = case$id,
      algorithm = "base::svd",
      backend = "cpu",
      precision = "strict",
      elapsed = ref_exact$elapsed,
      rel_sv_err = 0,
      selected_backend = "cpu"
    ),
    new_row(
      case = case$id,
      algorithm = "irlba::svdr",
      backend = "cpu",
      precision = "strict",
      elapsed = ref_rsvd$elapsed,
      rel_sv_err = relative_sv_error(ref_rsvd$result$d[seq_len(case$k)], ref_sv),
      iter = as.integer(ref_rsvd$result$iter %||% NA_integer_),
      mprod = as.integer(ref_rsvd$result$mprod %||% NA_integer_),
      selected_backend = "cpu"
    )
  )
}

svd_support_status <- function(x, backend_name, algorithm) {
  activate_worker_backend(backend_name)
  backend <- tryCatch(amatrix:::.amatrix_get_backend(backend_name), error = function(e) NULL)
  if (is.null(backend) || !isTRUE(backend$available())) {
    return(list(ok = FALSE, reason = "backend unavailable"))
  }

  if (identical(algorithm, "svd")) {
    has_svd <- "svd" %in% backend$capabilities() && is.function(backend$svd)
    if (!has_svd) {
      return(list(ok = FALSE, reason = "backend lacks svd capability"))
    }
    if (!isTRUE(backend$supports("svd", x))) {
      return(list(ok = FALSE, reason = "backend supports('svd') returned FALSE"))
    }
    return(list(ok = TRUE, reason = NA_character_))
  }

  has_rsvd <- "rsvd" %in% backend$capabilities() && is.function(backend$rsvd)
  if (!has_rsvd) {
    return(list(ok = FALSE, reason = "backend lacks rsvd capability"))
  }

  list(ok = TRUE, reason = NA_character_)
}

benchmark_backend_case <- function(case, backend_spec, reps = 3L) {
  host <- make_host_case(case)
  ref_sv <- base::svd(host, nu = case$k, nv = case$k)$d[seq_len(case$k)]

  backend_name <- backend_spec$name
  precision <- backend_spec$precision
  x <- adgeMatrix(host, preferred_backend = backend_name, precision = precision)
  rows <- list()

  exact_status <- svd_support_status(x, backend_name, "svd")
  if (isTRUE(exact_status$ok)) {
    exact_bench <- benchmark_elapsed(
      function() svd(x, nu = case$k, nv = case$k),
      reps = reps,
      warmup = if (!identical(backend_name, "cpu")) function() invisible(svd(x, nu = case$k, nv = case$k)) else NULL
    )
    plan <- tryCatch(amatrix_backend_plan(x, "svd"), error = function(e) NULL)
    selected_backend <- plan$chosen %||% backend_name
    reason <- NA_character_
    if (identical(backend_name, "opencl")) {
      selected_backend <- "cpu"
      reason <- "OpenCL exact svd currently uses host LAPACK"
    }
    rows[[length(rows) + 1L]] <- new_row(
      case = case$id,
      algorithm = "svd",
      backend = backend_name,
      precision = precision,
      elapsed = exact_bench$elapsed,
      rel_sv_err = relative_sv_error(exact_bench$result$d[seq_len(case$k)], ref_sv),
      selected_backend = selected_backend,
      reason = reason
    )
  } else {
    rows[[length(rows) + 1L]] <- new_row(
      case = case$id,
      algorithm = "svd",
      backend = backend_name,
      precision = precision,
      status = "unsupported",
      reason = exact_status$reason
    )
  }

  rsvd_status <- svd_support_status(x, backend_name, "rsvd")
  if (isTRUE(rsvd_status$ok)) {
    rsvd_bench <- benchmark_elapsed(
      function() rsvd(x, k = case$k, n_oversamples = case$n_oversamples, n_iter = case$n_iter),
      reps = reps,
      warmup = if (!identical(backend_name, "cpu")) function() invisible(rsvd(x, k = case$k, n_oversamples = case$n_oversamples, n_iter = case$n_iter)) else NULL
    )
    selected_backend <- if (identical(backend_name, "cpu")) {
      "cpu"
    } else {
      tryCatch(amatrix:::.amatrix_svd_factor_rsvd_backend(x), error = function(e) NULL) %||% "cpu"
    }
    rows[[length(rows) + 1L]] <- new_row(
      case = case$id,
      algorithm = "rsvd",
      backend = backend_name,
      precision = precision,
      elapsed = rsvd_bench$elapsed,
      rel_sv_err = relative_sv_error(rsvd_bench$result$d[seq_len(case$k)], ref_sv),
      iter = as.integer(rsvd_bench$result$iter %||% NA_integer_),
      mprod = as.integer(rsvd_bench$result$mprod %||% NA_integer_),
      selected_backend = selected_backend
    )
  } else {
    rows[[length(rows) + 1L]] <- new_row(
      case = case$id,
      algorithm = "rsvd",
      backend = backend_name,
      precision = precision,
      status = "unsupported",
      reason = rsvd_status$reason
    )
  }

  do.call(rbind, rows)
}

mlx_single_cell_expr <- function(case, algorithm, out_path, reps = 3L, native_spectral = FALSE) {
  helper_path <- normalizePath(file.path("tools", "benchmark-helpers.R"), winslash = "/", mustWork = TRUE)
  out_path <- normalizePath(out_path, winslash = "/", mustWork = FALSE)
  native_spectral <- isTRUE(native_spectral)
  case_args <- sprintf(
    "id = %s, n = %dL, p = %dL, k = %dL, n_oversamples = %dL, n_iter = %dL",
    r_string_literal(case$id),
    case$n,
    case$p,
    case$k,
    case$n_oversamples,
    case$n_iter
  )

  if (identical(algorithm, "svd")) {
    run_expr <- sprintf(
      "bench_eval(function() svd(x, nu = case$k, nv = case$k), reps = %dL)",
      reps
    )
    warmup_expr <- "invisible(svd(x, nu = case$k, nv = case$k))"
    iter_expr <- "NA_integer_"
    mprod_expr <- "NA_integer_"
  } else {
    run_expr <- sprintf(
      "bench_eval(function() rsvd(x, k = case$k, n_oversamples = case$n_oversamples, n_iter = case$n_iter), reps = %dL)",
      reps
    )
    warmup_expr <- "invisible(rsvd(x, k = case$k, n_oversamples = case$n_oversamples, n_iter = case$n_iter))"
    iter_expr <- "as.integer(bench$result$iter %||% NA_integer_)"
    mprod_expr <- "as.integer(bench$result$mprod %||% NA_integer_)"
  }

  paste(
    sprintf("source(%s, local = FALSE)", r_string_literal(helper_path)),
    "load_benchmark_amatrix()",
    "if (!requireNamespace(\"irlba\", quietly = TRUE)) stop(\"Package 'irlba' is required for this benchmark\", call. = FALSE)",
    "`%||%` <- function(x, y) if (is.null(x)) y else x",
    "relative_sv_error <- function(actual, expected) max(abs(actual - expected) / pmax(abs(expected), 1e-12))",
    "bench_eval <- function(fn, reps) { timings <- numeric(reps); last <- NULL; for (idx in seq_len(reps)) { gc(); timings[[idx]] <- system.time(last <- fn())[[\"elapsed\"]] }; list(elapsed = median(timings), result = last) }",
    "new_row <- function(...) { defaults <- list(case = NA_character_, algorithm = NA_character_, backend = NA_character_, precision = NA_character_, status = \"ok\", reason = NA_character_, elapsed = NA_real_, rel_sv_err = NA_real_, iter = NA_integer_, mprod = NA_integer_, selected_backend = NA_character_, stringsAsFactors = FALSE); as.data.frame(modifyList(defaults, list(...)), stringsAsFactors = FALSE) }",
    sprintf("case <- list(%s)", case_args),
    sprintf("algorithm <- %s", r_string_literal(algorithm)),
    sprintf("out_path <- %s", r_string_literal(out_path)),
    sprintf("native_spectral <- %s", if (native_spectral) "TRUE" else "FALSE"),
    "ns <- ensure_optional_backend_namespace(\"amatrix.mlx\", repo_dir = \"backends/amatrix.mlx\")",
    "if (is.null(ns)) { saveRDS(new_row(case = case$id, algorithm = algorithm, backend = \"mlx\", precision = \"fast\", status = \"unsupported\", reason = \"backend namespace unavailable\"), out_path); quit(save = \"no\", status = 0L) }",
    "options(amatrix.mlx.available = TRUE, amatrix.mlx.safe_spectral = !native_spectral, amatrix.mlx.rsvd.engine = \"resident\")",
    "if (isTRUE(native_spectral)) { Sys.unsetenv(\"AMATRIX_MLX_SAFE_SPECTRAL\"); Sys.setenv(AMATRIX_MLX_NATIVE_SPECTRAL = \"1\") } else { Sys.setenv(AMATRIX_MLX_SAFE_SPECTRAL = \"1\"); Sys.unsetenv(\"AMATRIX_MLX_NATIVE_SPECTRAL\") }",
    "get(\"amatrix_mlx_register\", envir = ns, inherits = FALSE)(overwrite = TRUE)",
    "set.seed(20260409L + case$n + case$p + case$k)",
    "host <- matrix(rnorm(case$n * case$p), nrow = case$n, ncol = case$p)",
    "ref_sv <- base::svd(host, nu = case$k, nv = case$k)$d[seq_len(case$k)]",
    "x <- adgeMatrix(host, preferred_backend = \"mlx\", precision = \"fast\")",
    warmup_expr,
    sprintf("bench <- %s", run_expr),
    "safe_spectral <- isTRUE(getOption(\"amatrix.mlx.safe_spectral\", FALSE)) || identical(Sys.getenv(\"AMATRIX_MLX_SAFE_SPECTRAL\", unset = \"\"), \"1\")",
    "selected_backend <- if (safe_spectral) \"cpu\" else if (identical(algorithm, \"svd\")) \"cpu\" else \"mlx\"",
    "reason <- if (safe_spectral) \"safe spectral fallback to cpu\" else if (identical(algorithm, \"svd\")) \"MLX exact svd uses CPU stream\" else NA_character_",
    "row <- new_row(case = case$id, algorithm = algorithm, backend = \"mlx\", precision = \"fast\", elapsed = bench$elapsed, rel_sv_err = relative_sv_error(bench$result$d[seq_len(case$k)], ref_sv), iter = NULL, mprod = NULL, selected_backend = selected_backend, reason = reason)",
    sprintf("row$iter <- %s", iter_expr),
    sprintf("row$mprod <- %s", mprod_expr),
    "saveRDS(row, out_path)",
    sep = "; "
  )
}

run_mlx_single_cell <- function(case, algorithm, reps = 3L, native_spectral = FALSE) {
  out_path <- tempfile(sprintf("%s-%s-mlx-", case$id, algorithm), fileext = ".rds")
  on.exit(unlink(out_path), add = TRUE)

  launch <- benchmark_system2_capture(
    file.path(R.home("bin"), "Rscript"),
    c("-e", mlx_single_cell_expr(case, algorithm, out_path, reps = reps, native_spectral = native_spectral))
  )

  if (launch$status == 0L && file.exists(out_path)) {
    return(readRDS(out_path))
  }

  new_row(
    case = case$id,
    algorithm = algorithm,
    backend = "mlx",
    precision = "fast",
    status = "crash",
    reason = sprintf(
      "worker exited with status %s%s",
      launch$status,
      if (length(launch$output) > 0L) paste0(": ", tail(launch$output, 1L)) else ""
    )
  )
}

benchmark_mlx_inline_rsvd <- function(case, reps = 3L) {
  activate_worker_backend("mlx")
  old_options <- options(
    amatrix.mlx.available = TRUE,
    amatrix.mlx.safe_spectral = FALSE,
    amatrix.mlx.rsvd.engine = "resident"
  )
  old_safe <- Sys.getenv("AMATRIX_MLX_SAFE_SPECTRAL", unset = NA_character_)
  old_native <- Sys.getenv("AMATRIX_MLX_NATIVE_SPECTRAL", unset = NA_character_)
  on.exit(options(old_options), add = TRUE)
  on.exit({
    if (is.na(old_safe)) Sys.unsetenv("AMATRIX_MLX_SAFE_SPECTRAL") else Sys.setenv(AMATRIX_MLX_SAFE_SPECTRAL = old_safe)
    if (is.na(old_native)) Sys.unsetenv("AMATRIX_MLX_NATIVE_SPECTRAL") else Sys.setenv(AMATRIX_MLX_NATIVE_SPECTRAL = old_native)
  }, add = TRUE)
  Sys.unsetenv("AMATRIX_MLX_SAFE_SPECTRAL")
  Sys.setenv(AMATRIX_MLX_NATIVE_SPECTRAL = "1")

  tryCatch({
    host <- make_host_case(case)
    ref_sv <- base::svd(host, nu = case$k, nv = case$k)$d[seq_len(case$k)]
    x <- adgeMatrix(host, preferred_backend = "mlx", precision = "fast")
    bench <- benchmark_elapsed(
      function() rsvd(x, k = case$k, n_oversamples = case$n_oversamples, n_iter = case$n_iter),
      reps = reps,
      warmup = function() invisible(rsvd(x, k = case$k, n_oversamples = case$n_oversamples, n_iter = case$n_iter))
    )

    new_row(
      case = case$id,
      algorithm = "rsvd",
      backend = "mlx",
      precision = "fast",
      elapsed = bench$elapsed,
      rel_sv_err = relative_sv_error(bench$result$d[seq_len(case$k)], ref_sv),
      iter = as.integer(bench$result$iter %||% NA_integer_),
      mprod = as.integer(bench$result$mprod %||% NA_integer_),
      selected_backend = "mlx",
      reason = "native inline resident rsvd"
    )
  }, error = function(e) {
    new_row(
      case = case$id,
      algorithm = "rsvd",
      backend = "mlx",
      precision = "fast",
      status = "error",
      reason = conditionMessage(e)
    )
  })
}

crash_rows <- function(case, backend_spec, message_text) {
  rbind(
    new_row(
      case = case$id,
      algorithm = "svd",
      backend = backend_spec$name,
      precision = backend_spec$precision,
      status = "crash",
      reason = message_text
    ),
    new_row(
      case = case$id,
      algorithm = "rsvd",
      backend = backend_spec$name,
      precision = backend_spec$precision,
      status = "crash",
      reason = message_text
    )
  )
}

worker_main <- function(args) {
  if (is.null(args$plan) || is.null(args$out) || is.null(args$entry_id)) {
    stop("worker mode requires --plan, --entry-id, and --out", call. = FALSE)
  }

  entries <- readRDS(args$plan)
  entry <- Filter(function(x) identical(x$id, args$entry_id), entries)[[1L]]
  result <- benchmark_backend_case(entry$case, entry$backend, reps = 3L)
  saveRDS(result, args$out)
  invisible(NULL)
}

master_main <- function(args = parse_args(commandArgs(trailingOnly = TRUE))) {
  cases <- default_cases()
  backends <- available_svd_backends()

  rows <- list()
  for (case in cases) {
    rows[[length(rows) + 1L]] <- benchmark_reference_case(case, reps = 3L)
  }

  entries <- list()
  for (case in cases) {
    for (backend in backends) {
      entries[[length(entries) + 1L]] <- list(
        id = sprintf("%s-%s", case$id, backend$name),
        case = case,
        backend = backend
      )
    }
  }

  plan_path <- tempfile("benchmark-svd-plan-", fileext = ".rds")
  saveRDS(entries, plan_path)
  on.exit(unlink(plan_path), add = TRUE)

  script_path <- normalizePath(sys.frame(1)$ofile %||% "tools/benchmark-svd-backends.R", mustWork = FALSE)

  for (entry in entries) {
    message(sprintf("  %s", entry$id))
    if (identical(entry$backend$name, "mlx")) {
      # Exact MLX SVD is CPU-stream only and remains on the safe path. Native
      # MLX spectral probing is useful for RSVD, where the resident algorithm is
      # dominated by GPU matmul/crossprod.
      rows[[length(rows) + 1L]] <- run_mlx_single_cell(entry$case, "svd", reps = 3L, native_spectral = FALSE)
      rows[[length(rows) + 1L]] <- if (isTRUE(args$mlx_native_spectral) && isTRUE(args$mlx_native_inline)) {
        benchmark_mlx_inline_rsvd(entry$case, reps = 3L)
      } else {
        run_mlx_single_cell(entry$case, "rsvd", reps = 3L, native_spectral = args$mlx_native_spectral)
      }
      next
    }

    out_path <- tempfile(sprintf("%s-", entry$id), fileext = ".rds")
    launch <- benchmark_system2_capture(
      file.path(R.home("bin"), "Rscript"),
      benchmark_rscript_source_args(
        script_path,
        args = c("--worker", paste0("--plan=", plan_path), paste0("--entry-id=", entry$id), paste0("--out=", out_path))
      )
    )
    log <- launch$output
    status <- launch$status

    if (status == 0L && file.exists(out_path)) {
      rows[[length(rows) + 1L]] <- readRDS(out_path)
    } else {
      rows[[length(rows) + 1L]] <- crash_rows(
        entry$case,
        entry$backend,
        sprintf("worker exited with status %s%s", status, if (length(log) > 0L) paste0(": ", tail(log, 1L)) else "")
      )
    }

    unlink(out_path)
  }

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

add_cpu_reference <- function(rows) {
  cpu <- rows[
    rows$status == "ok" & rows$backend == "cpu",
    c("case", "algorithm", "elapsed")
  ]
  names(cpu)[names(cpu) == "elapsed"] <- "cpu_reference_elapsed"
  merge(rows, cpu, by = c("case", "algorithm"), all.x = TRUE)
}

summarize_svd_benchmark <- function(rows) {
  out <- add_cpu_reference(rows)
  out$speedup_vs_cpu <- ifelse(
    is.na(out$cpu_reference_elapsed) | is.na(out$elapsed),
    NA_real_,
    out$cpu_reference_elapsed / out$elapsed
  )
  out$dispatch_state <- ifelse(
    out$status != "ok",
    out$status,
    ifelse(
      out$backend == "cpu",
      "cpu_baseline",
      ifelse(
        is.na(out$selected_backend),
        "unknown_dispatch",
        ifelse(out$selected_backend == out$backend, "accelerated", ifelse(out$selected_backend == "cpu", "cpu_fallback", "rerouted"))
      )
    )
  )
  out[order(out$case, out$algorithm, out$backend), ]
}

write_svd_outputs <- function(rows, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  summary <- summarize_svd_benchmark(rows)
  incidents <- summary[summary$status != "ok", , drop = FALSE]
  fallbacks <- summary[
    summary$status == "ok" &
      summary$backend != "cpu" &
      summary$dispatch_state == "cpu_fallback",
    ,
    drop = FALSE
  ]
  raw_path <- file.path(output_dir, "raw-results.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  incidents_path <- file.path(output_dir, "incidents.csv")
  fallbacks_path <- file.path(output_dir, "fallbacks.csv")
  metadata_path <- file.path(output_dir, "metadata.rds")

  write.csv(rows, raw_path, row.names = FALSE)
  write.csv(summary, summary_path, row.names = FALSE)
  write.csv(incidents, incidents_path, row.names = FALSE)
  write.csv(fallbacks, fallbacks_path, row.names = FALSE)
  saveRDS(
    list(
      created_at = Sys.time(),
      hostname = Sys.info()[["nodename"]],
      r_version = paste(R.version$major, R.version$minor, sep = "."),
      platform = R.version$platform,
      timezone = Sys.timezone()
    ),
    metadata_path
  )

  list(
    raw = raw_path,
    summary = summary_path,
    incidents = incidents_path,
    fallbacks = fallbacks_path,
    metadata = metadata_path
  )
}

print_svd_benchmark <- function(output_dir) {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  rows <- master_main(args)
  paths <- write_svd_outputs(rows, output_dir)
  summary <- summarize_svd_benchmark(rows)

  cat("Notes:\n")
  cat("- `svd` benchmarks exact singular-value decomposition through the public amatrix surface.\n")
  cat("- `rsvd` benchmarks the standalone truncated randomized SVD surface.\n")
  cat("- MLX exact SVD stays on the safe CPU fallback path because `mlx_linalg_svd` is CPU-stream only in the current bridge.\n")
  cat("- MLX RSVD uses safe CPU fallback by default; set `AMATRIX_MLX_NATIVE_SPECTRAL=1` or pass `--mlx-native-spectral` to crash-probe native RSVD in isolated workers.\n")
  cat("- For actual native MLX RSVD timings in the top-level process, add `AMATRIX_MLX_NATIVE_INLINE=1` or `--mlx-native-inline`; this is intentionally opt-in because it does not protect the master process from native aborts.\n")
  cat("- `selected_backend` records the backend actually used by the path; `unsupported` and `crash` rows make backend gaps explicit.\n")
  cat("- `dispatch_state` labels accelerated rows separately from CPU fallbacks.\n")
  cat("- `rel_sv_err` is measured against leading singular values from `base::svd()` on the same host matrix.\n\n")

  print(summary, row.names = FALSE)
  cat("\nArtifacts:\n")
  cat("  raw: ", paths$raw, "\n", sep = "")
  cat("  summary: ", paths$summary, "\n", sep = "")
  cat("  incidents: ", paths$incidents, "\n", sep = "")
  cat("  fallbacks: ", paths$fallbacks, "\n", sep = "")
  cat("  metadata: ", paths$metadata, "\n", sep = "")
  invisible(summary)
}

relaunch_safe_master_if_needed <- function(args) {
  raw_args <- commandArgs(trailingOnly = FALSE)
  direct_file_arg <- grep("^--file=", raw_args, value = TRUE)
  direct_file_entry <- length(direct_file_arg) > 0L
  safe_main <- identical(Sys.getenv("AMATRIX_SVD_BENCH_SAFE_MAIN", unset = ""), "true")

  if (!direct_file_entry || isTRUE(args$worker) || safe_main) {
    return(invisible(FALSE))
  }

  script_path <- normalizePath(sub("^--file=", "", direct_file_arg[[1L]]), winslash = "/", mustWork = TRUE)
  repo_root <- normalizePath(dirname(dirname(script_path)), winslash = "/", mustWork = TRUE)
  expr <- sprintf(
    "Sys.setenv(AMATRIX_SVD_BENCH_SAFE_MAIN = \"true\"); setwd(%s); source(%s, local = globalenv())",
    r_string_literal(repo_root),
    r_string_literal(script_path)
  )
  trailing_args <- commandArgs(trailingOnly = TRUE)
  if (length(trailing_args) > 0L && identical(trailing_args[[1L]], "--args")) {
    trailing_args <- trailing_args[-1L]
  }
  relaunch_args <- c("-e", expr, if (length(trailing_args) > 0L) c("--args", trailing_args))
  quoted_args <- vapply(relaunch_args, shQuote, character(1), USE.NAMES = FALSE)
  warned_status <- NULL
  relaunch_output <- withCallingHandlers(
    system2(file.path(R.home("bin"), "Rscript"), quoted_args, stdout = TRUE, stderr = TRUE),
    warning = function(w) {
      warned_status <<- attr(w, "status") %||% warned_status
      invokeRestart("muffleWarning")
    }
  )
  launch <- list(output = relaunch_output, status = attr(relaunch_output, "status") %||% warned_status %||% 0L)

  if (length(launch$output) > 0L) {
    cat(paste(launch$output, collapse = "\n"), sep = "\n")
    if (!grepl("\n$", paste(launch$output, collapse = "\n"))) {
      cat("\n")
    }
  }

  quit(save = "no", status = launch$status)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
relaunch_safe_master_if_needed(args)
initialize_svd_benchmark_context()

if (isTRUE(args$worker)) {
  worker_main(args)
} else if (!identical(Sys.getenv("AMATRIX_BENCHMARK_NO_AUTORUN", unset = ""), "1")) {
  print_svd_benchmark(args$output_dir)
}
