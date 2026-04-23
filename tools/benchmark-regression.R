#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x
r_string_literal <- function(x) encodeString(x, quote = "\"")

benchmark_regression_cli_script_path <- function(script_path = NULL) {
  if (is.null(script_path)) {
    script_path <- getOption("amatrix.benchmark_regression.script_path", "tools/benchmark-regression.R")
  }

  repo_root <- normalizePath(dirname(dirname(script_path)), winslash = "/", mustWork = TRUE)
  cli_path <- normalizePath(file.path(repo_root, "tools", "benchmark-regression-cli.R"), winslash = "/", mustWork = TRUE)

  list(
    script_path = normalizePath(script_path, winslash = "/", mustWork = TRUE),
    repo_root = repo_root,
    cli_path = cli_path
  )
}

benchmark_regression_bootstrap_direct_entry <- function() {
  raw_args <- commandArgs(trailingOnly = FALSE)
  direct_file_paths <- sub("^--file=", "", grep("^--file=", raw_args, value = TRUE))
  if (length(direct_file_paths) == 0L) {
    return(invisible(FALSE))
  }
  direct_file_paths <- direct_file_paths[basename(direct_file_paths) == "benchmark-regression.R"]
  if (length(direct_file_paths) == 0L) {
    return(invisible(FALSE))
  }

  trailing_args <- commandArgs(trailingOnly = TRUE)
  if ("--safe-main" %in% trailing_args || "--worker" %in% trailing_args) {
    return(invisible(FALSE))
  }

  cli_info <- benchmark_regression_cli_script_path(direct_file_paths[[1L]])
  relaunch_args <- c(cli_info$cli_path, trailing_args)
  warned_status <- NULL
  relaunch_output <- withCallingHandlers(
    system2(file.path(R.home("bin"), "Rscript"), vapply(relaunch_args, shQuote, character(1), USE.NAMES = FALSE), stdout = TRUE, stderr = TRUE),
    warning = function(w) {
      warned_status <<- attr(w, "status") %||% warned_status
      invokeRestart("muffleWarning")
    }
  )
  status <- attr(relaunch_output, "status") %||% warned_status %||% 0L

  if (length(relaunch_output) > 0L) {
    cat(paste(relaunch_output, collapse = "\n"), sep = "\n")
    if (!grepl("\n$", paste(relaunch_output, collapse = "\n"))) {
      cat("\n")
    }
  }

  quit(save = "no", status = status)
}

benchmark_regression_bootstrap_direct_entry()

timestamp_tag <- function(x = Sys.time()) {
  format(as.POSIXct(x, tz = Sys.timezone()), "%Y%m%d-%H%M%S")
}

benchmark_launch_debug <- function(...) {
  path <- Sys.getenv("AMATRIX_BENCHMARK_LAUNCH_DEBUG", unset = "")
  if (!nzchar(path)) {
    return(invisible(FALSE))
  }

  line <- paste(..., collapse = "")
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), line), file = path, append = TRUE)
  invisible(TRUE)
}

benchmark_regression_usage <- function() {
  paste(
    "amatrix benchmark harness",
    "",
    "Convenience entry point:",
    "  bash tools/benchmark-regression.sh [options]",
    "",
    "Stable direct entry point:",
    "  Rscript tools/benchmark-regression-cli.R [options]",
    "",
    "Historical direct entry point:",
    "  Rscript tools/benchmark-regression.R [options]",
    "",
    "Options:",
    "  --help                    Show this help and exit",
    "  --update                  Refresh tools/baseline.csv from the current run",
    "  --against-baseline        Explicit no-op alias for compare-against-baseline mode",
    "  --warm-only               Shorten the console summary to warm/resident rows",
    "  --baseline=PATH           Baseline CSV to compare against",
    "  --output-dir=PATH         Result bundle directory (default: tools/benchmark-results/<timestamp>)",
    "  --suites=dense,sparse     Suites to run",
    "  --backends=cpu,opencl     Requested backends",
    "  --include-arrayfire=1     Include ArrayFire during backend discovery",
    "",
    "Result bundle:",
    "  raw-results.csv           Per-cell measurements",
    "  summary.csv               Decorated summary with routing/performance states",
    "  regressions.csv           Cells slower than baseline by >20%",
    "  warm-ratios.csv           Cold-vs-warm/resident comparisons",
    "  routing-summary.csv       Backend/status/routing counts",
    "  benchmark-summary.md      Human-readable markdown summary",
    "  benchmark-report.qmd      Quarto source report",
    "  benchmark-report.html     Rendered HTML report when Quarto is available",
    "  benchmark-report.pdf      Rendered PDF report when Quarto is available",
    "  plots/*.png               Baseline, warm-state, and routing charts",
    "",
    "Notes:",
    "  - Baseline comparison is on by default when the baseline file exists.",
    "  - On Apple Silicon, prefer cpu/mlx/metal for canonical runs.",
    "  - ArrayFire on Apple is diagnostic-only and may require the direct Rscript entry path.",
    sep = "\n"
  )
}

parse_args <- function(args) {
  out <- list(
    help = any(args %in% c("--help", "-h")),
    update = "--update" %in% args,
    against_baseline = "--against-baseline" %in% args,
    warm_only = "--warm-only" %in% args,
    worker = "--worker" %in% args,
    safe_main = "--safe-main" %in% args,
    baseline = "tools/baseline.csv",
    output_dir = file.path("tools", "benchmark-results", timestamp_tag()),
    plan = NULL,
    out = NULL,
    group_id = NULL,
    suites = c("dense", "sparse"),
    backends = NULL,
    include_arrayfire = identical(Sys.getenv("AMATRIX_BENCHMARK_ARRAYFIRE", unset = ""), "1") ||
      identical(Sys.getenv("AMATRIX_ARRAYFIRE_PROBE_GPU", unset = ""), "1")
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
    if (identical(name, "baseline")) out$baseline <- value
    if (identical(name, "output-dir")) out$output_dir <- value
    if (identical(name, "plan")) out$plan <- value
    if (identical(name, "out")) out$out <- value
    if (identical(name, "group-id")) out$group_id <- value
    if (identical(name, "suites")) out$suites <- strsplit(value, ",", fixed = TRUE)[[1L]]
    if (identical(name, "backends")) out$backends <- strsplit(value, ",", fixed = TRUE)[[1L]]
    if (identical(name, "include-arrayfire")) out$include_arrayfire <- identical(value, "1")
  }

  out$suites <- unique(out$suites[nzchar(out$suites)])
  out
}

relaunch_safe_master_if_needed <- function(args) {
  raw_args <- commandArgs(trailingOnly = FALSE)
  direct_file_paths <- sub("^--file=", "", grep("^--file=", raw_args, value = TRUE))
  direct_file_paths <- direct_file_paths[basename(direct_file_paths) == "benchmark-regression.R"]
  if (length(direct_file_paths) == 0L) {
    script_like_paths <- raw_args[basename(raw_args) == "benchmark-regression.R"]
    if (length(script_like_paths) > 0L) {
      direct_file_paths <- script_like_paths
    } else {
      existing_paths <- raw_args[file.exists(raw_args)]
      direct_file_paths <- existing_paths[basename(existing_paths) == "benchmark-regression.R"]
    }
  }
  direct_file_entry <- length(direct_file_paths) > 0L
  benchmark_launch_debug(
    "raw_args=", paste(raw_args, collapse = " | "),
    " ; direct_file_paths=", paste(direct_file_paths, collapse = " | "),
    " ; direct_file_entry=", direct_file_entry,
    " ; worker=", isTRUE(args$worker),
    " ; safe_main=", isTRUE(args$safe_main)
  )

  if (!direct_file_entry || isTRUE(args$worker) || isTRUE(args$safe_main)) {
    benchmark_launch_debug("skipping relaunch")
    return(invisible(FALSE))
  }

  cli_info <- benchmark_regression_cli_script_path(direct_file_paths[[1L]])
  relaunch_args <- c(cli_info$cli_path, commandArgs(trailingOnly = TRUE))
  quoted_relaunch_args <- vapply(relaunch_args, shQuote, character(1), USE.NAMES = FALSE)
  benchmark_launch_debug("relaunching via cli wrapper with script_path=", cli_info$script_path)
  warned_status <- NULL
  relaunch_output <- withCallingHandlers(
    system2(file.path(R.home("bin"), "Rscript"), quoted_relaunch_args, stdout = TRUE, stderr = TRUE),
    warning = function(w) {
      warned_status <<- attr(w, "status") %||% warned_status
      invokeRestart("muffleWarning")
    }
  )
  status <- attr(relaunch_output, "status") %||% warned_status %||% 0L

  if (length(relaunch_output) > 0L) {
    cat(paste(relaunch_output, collapse = "\n"), sep = "\n")
    if (!grepl("\n$", paste(relaunch_output, collapse = "\n"))) {
      cat("\n")
    }
  }

  quit(save = "no", status = status)
}

initialize_regression_benchmark_context <- local({
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
      if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Matrix package required for sparse benchmark suite", call. = FALSE)
      }
    })

    initialized <<- TRUE
    invisible(TRUE)
  }
})

canonical_backend_specs <- function(include_arrayfire = .benchmark_arrayfire_requested(), only = NULL) {
  requested <- if (is.null(only)) {
    NULL
  } else {
    unique(only[nzchar(only)])
  }

  available_benchmark_backends(
    include_cpu = is.null(requested) || "cpu" %in% requested,
    include_mlx = (is.null(requested) || "mlx" %in% requested) && benchmark_worker_mlx_enabled(),
    include_metal = is.null(requested) || "metal" %in% requested,
    include_opencl = is.null(requested) || "opencl" %in% requested,
    include_arrayfire = (is.null(requested) || "arrayfire" %in% requested) &&
      benchmark_worker_arrayfire_enabled(include_arrayfire)
  )
}

benchmark_running_on_apple_silicon <- function() {
  identical(Sys.info()[["sysname"]], "Darwin") &&
    grepl("arm64|aarch64", R.version$arch, ignore.case = TRUE)
}

benchmark_metadata_on_apple_silicon <- function(metadata) {
  platform <- tolower(as.character(metadata$platform %||% ""))
  sysname <- as.character(metadata$sysname %||% "")
  arch <- tolower(as.character(metadata$arch %||% ""))

  grepl("apple-darwin", platform, fixed = TRUE) ||
    (identical(sysname, "Darwin") && grepl("arm64|aarch64", arch))
}

benchmark_worker_mlx_enabled <- function() {
  flag <- Sys.getenv("AMATRIX_BENCHMARK_MLX_WORKERS", unset = "")
  if (identical(flag, "1")) {
    return(TRUE)
  }
  if (identical(flag, "0")) {
    return(FALSE)
  }

  benchmark_running_on_apple_silicon()
}

benchmark_worker_arrayfire_enabled <- function(include_arrayfire = .benchmark_arrayfire_requested()) {
  if (!isTRUE(include_arrayfire)) {
    return(FALSE)
  }

  if (benchmark_running_on_apple_silicon() &&
      !identical(Sys.getenv("AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE", unset = ""), "1")) {
    return(FALSE)
  }

  TRUE
}

benchmark_worker_backend_allowed <- function(backend_name, include_arrayfire = .benchmark_arrayfire_requested()) {
  if (identical(backend_name, "mlx")) {
    return(benchmark_worker_mlx_enabled())
  }

  if (identical(backend_name, "arrayfire")) {
    return(benchmark_worker_arrayfire_enabled(include_arrayfire))
  }

  TRUE
}

benchmark_worker_backend_block_reason <- function(backend_name, include_arrayfire = .benchmark_arrayfire_requested()) {
  if (benchmark_worker_backend_allowed(backend_name, include_arrayfire = include_arrayfire)) {
    return(NA_character_)
  }

  if (identical(backend_name, "mlx")) {
    if (benchmark_running_on_apple_silicon()) {
      return("MLX disabled via AMATRIX_BENCHMARK_MLX_WORKERS=0")
    }
    return("MLX canonical workers require Apple Silicon or AMATRIX_BENCHMARK_MLX_WORKERS=1")
  }

  if (identical(backend_name, "arrayfire") &&
      benchmark_running_on_apple_silicon() &&
      !identical(Sys.getenv("AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE", unset = ""), "1")) {
    return("ArrayFire is diagnostic-only on Apple benchmark workers unless AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE=1")
  }

  sprintf("backend '%s' is disabled in canonical regression workers", backend_name)
}

filter_backend_specs <- function(specs, names = NULL) {
  if (is.null(names)) {
    return(specs)
  }

  keep <- names[names %in% vapply(specs, `[[`, character(1), "name")]
  specs[vapply(specs, function(spec) spec$name %in% keep, logical(1))]
}

prime_requested_backend <- function(backend_name, include_arrayfire = .benchmark_arrayfire_requested()) {
  if (identical(backend_name, "cpu")) {
    return(invisible(TRUE))
  }

  if (!benchmark_worker_backend_allowed(backend_name, include_arrayfire = include_arrayfire)) {
    return(invisible(FALSE))
  }

  specs <- .benchmark_optional_backend_specs(include_arrayfire = include_arrayfire)
  spec <- specs[[backend_name]]
  if (is.null(spec)) {
    return(invisible(FALSE))
  }

  invisible(.benchmark_enable_backend(spec))
}

dense_sizes <- function() {
  list(
    small = list(n = 256L, p = 32L, sink_n = 128L),
    medium = list(n = 1024L, p = 128L, sink_n = 512L),
    large = list(n = 4096L, p = 128L, sink_n = 1024L),
    xlarge = list(n = 4096L, p = 1024L, sink_n = 2048L)
  )
}

sparse_sizes <- function() {
  list(
    medium = list(nrow = 4000L, ncol = 1000L),
    large = list(nrow = 8000L, ncol = 2000L)
  )
}

sparse_densities <- function() c(0.001, 0.005, 0.01, 0.05)

spmm_rhs_widths <- function() c(8L, 32L)

sparse_iterative_ops <- function() c("block_lanczos", "svd_factor_subspace")

sparse_iterative_densities <- function() c(0.05)

sparse_iterative_params <- function(op, size_label) {
  size <- sparse_sizes()[[size_label]]
  if (is.null(size)) {
    stop(sprintf("unknown sparse size '%s'", size_label), call. = FALSE)
  }

  switch(
    op,
    block_lanczos = list(
      k = 8L,
      block_size = min(16L, size$nrow, size$ncol),
      n_steps = 4L
    ),
    svd_factor_subspace = list(
      k = 8L,
      n_oversamples = 8L,
      n_iter = 1L
    ),
    stop(sprintf("unknown sparse iterative op '%s'", op), call. = FALSE)
  )
}

new_result_row <- function(...) {
  defaults <- list(
    suite = NA_character_,
    op = NA_character_,
    size_label = NA_character_,
    variant = NA_character_,
    requested_backend = NA_character_,
    dispatch_probe_op = NA_character_,
    requested_supported = NA,
    requested_support_reason = NA_character_,
    dispatch_backend = NA_character_,
    dispatch_path = NA_character_,
    status = "ok",
    error_message = NA_character_,
    nrow = NA_integer_,
    ncol = NA_integer_,
    rhs_width = NA_integer_,
    nnz = NA_integer_,
    density = NA_real_,
    density_bucket = NA_character_,
    reps = NA_integer_,
    median_ms = NA_real_,
    mean_ms = NA_real_,
    sd_ms = NA_real_,
    p05_ms = NA_real_,
    p95_ms = NA_real_,
    n_reps = NA_integer_,
    rel_err = NA_real_
  )
  values <- modifyList(defaults, list(...))
  as.data.frame(values, stringsAsFactors = FALSE)
}

benchmark_time_ms <- function(fn, reps = 7L, warmup = 1L) {
  reps <- as.integer(reps)
  warmup <- as.integer(warmup)

  if (warmup > 0L) {
    for (idx in seq_len(warmup)) {
      fn()
    }
  }

  timings <- vapply(seq_len(reps), function(idx) {
    gc()
    start <- proc.time()[["elapsed"]]
    fn()
    (proc.time()[["elapsed"]] - start) * 1000
  }, numeric(1))

  list(
    median  = stats::median(timings),
    mean    = mean(timings),
    sd      = if (length(timings) > 1L) stats::sd(timings) else 0,
    p05     = unname(stats::quantile(timings, 0.05, names = FALSE, na.rm = TRUE)),
    p95     = unname(stats::quantile(timings, 0.95, names = FALSE, na.rm = TRUE)),
    n_reps  = length(timings),
    samples = timings
  )
}

make_dense_operand <- function(host, backend) {
  precision <- if (identical(backend, "cpu")) "strict" else "fast"
  adgeMatrix(host, preferred_backend = backend, precision = precision)
}

make_sparse_operand <- function(host, backend) {
  precision <- if (identical(backend, "cpu")) "strict" else "fast"
  as_adgCMatrix(host, preferred_backend = backend, precision = precision)
}

with_sparse_benchmark_backend_options <- function(requested_backend, rhs_width, code) {
  rhs_width <- as.integer(rhs_width)
  old <- options()
  on.exit(options(old), add = TRUE)

  if (identical(requested_backend, "mlx")) {
    options(
      amatrix.mlx.spmv_min_nnz = 1L,
      amatrix.mlx.spmm_min_nnz = 1L
    )
  }

  if (identical(requested_backend, "metal")) {
    options(
      amatrix.metal.spmv_min_nnz = 1L,
      amatrix.metal.spmm_min_nnz = 1L
    )
  }

  if (identical(requested_backend, "opencl")) {
    options(
      amatrix.opencl.spmv_min_nnz = 1L,
      amatrix.opencl.spmm_min_nnz = 1L
    )
  }

  if (identical(requested_backend, "arrayfire")) {
    options(
      amatrix.arrayfire.spmv_min_nnz = 1L,
      amatrix.arrayfire.spmm_min_nnz = 1L
    )
  }

  force(code)
}

release_residency <- function(x) {
  if (!inherits(x, "aMatrix")) {
    return(invisible(FALSE))
  }

  entry <- amatrix:::.amatrix_resident_entry(x)
  if (is.null(entry)) {
    return(invisible(FALSE))
  }

  backend <- tryCatch(amatrix:::.amatrix_get_backend(entry$backend), error = function(e) NULL)
  if (is.null(backend)) {
    return(invisible(FALSE))
  }

  if (isTRUE(entry$sparse) && is.function(backend$sparse_resident_drop)) {
    try(backend$sparse_resident_drop(entry$resident_key), silent = TRUE)
  } else if (is.function(backend$resident_drop)) {
    try(backend$resident_drop(entry$resident_key), silent = TRUE)
  }
  amatrix:::.amatrix_drop_resident_binding(x)
  invisible(TRUE)
}

drop_svd_factor_cache <- function(x, k, n_oversamples, n_iter) {
  if (!inherits(x, "aMatrix")) {
    return(invisible(FALSE))
  }

  plan <- tryCatch(
    amatrix:::.amatrix_svd_factor_plan(
      X = x,
      k = as.integer(k),
      method = "subspace",
      n_oversamples = as.integer(n_oversamples),
      n_iter = as.integer(n_iter)
    ),
    error = function(e) NULL
  )
  if (is.null(plan)) {
    return(invisible(FALSE))
  }

  key <- tryCatch(amatrix:::.amatrix_svd_cache_key(x, as.integer(k), plan), error = function(e) NULL)
  if (is.null(key) || !nzchar(key)) {
    return(invisible(FALSE))
  }

  if (exists(key, envir = amatrix:::.amatrix_state$model_cache, inherits = FALSE)) {
    rm(list = key, envir = amatrix:::.amatrix_state$model_cache)
  }
  if (exists(key, envir = amatrix:::.amatrix_state$cache_atime, inherits = FALSE)) {
    rm(list = key, envir = amatrix:::.amatrix_state$cache_atime)
  }

  invisible(TRUE)
}

compute_dispatch_info <- function(x, op, y = NULL, requested_backend = NULL) {
  plan <- tryCatch(amatrix_backend_plan(x, op, y = y), error = function(e) NULL)
  if (is.null(plan)) {
    return(list(
      dispatch_probe_op = op,
      requested_supported = NA,
      requested_support_reason = "dispatch plan unavailable",
      dispatch_backend = NA_character_,
      dispatch_path = NA_character_
    ))
  }

  requested_backend <- requested_backend %||% x@preferred_backend %||% NA_character_
  requested_candidate <- NULL
  if (!is.na(requested_backend)) {
    candidate_idx <- match(requested_backend, vapply(plan$candidates, `[[`, character(1), "name"))
    if (!is.na(candidate_idx)) {
      requested_candidate <- plan$candidates[[candidate_idx]]
    }
  }

  requested_supported <- NA
  requested_support_reason <- NA_character_
  if (!is.null(requested_candidate)) {
    requested_supported <- isTRUE(requested_candidate$supported)
    requested_support_reason <- if (!isTRUE(requested_candidate$registered)) {
      "backend not registered"
    } else if (!isTRUE(requested_candidate$available)) {
      "backend unavailable"
    } else if (!isTRUE(requested_candidate$precision_compatible)) {
      "precision incompatible"
    } else if (isTRUE(requested_candidate$supported_resident)) {
      "resident supported"
    } else if (isTRUE(requested_candidate$supported_cold) && !isTRUE(requested_candidate$calibration_ok)) {
      "calibration rejected"
    } else if (isTRUE(requested_candidate$supported_cold)) {
      "cold supported"
    } else if (isTRUE(requested_candidate$resident_active)) {
      "resident op unsupported"
    } else {
      "op unsupported"
    }
  }

  list(
    dispatch_probe_op = op,
    requested_supported = requested_supported,
    requested_support_reason = requested_support_reason,
    dispatch_backend = plan$chosen %||% NA_character_,
    dispatch_path = plan$chosen_path %||% NA_character_
  )
}

run_dense_case <- function(requested_backend, op, size_label, variant) {
  size <- dense_sizes()[[size_label]]
  if (is.null(size)) {
    stop(sprintf("unknown dense size '%s'", size_label), call. = FALSE)
  }

  set.seed(1L)
  X <- matrix(rnorm(size$n * size$p), size$n, size$p)
  set.seed(2L)
  Y <- matrix(rnorm(size$n * 3L), size$n, 3L)
  set.seed(3L)
  B_host <- matrix(rnorm(size$p * size$p), size$p, size$p)
  set.seed(4L)
  Z <- matrix(rnorm(size$p * size$p), size$p, size$p)
  SPD <- crossprod(Z) + diag(size$p) * 2
  SPD_rhs <- matrix(rnorm(size$p * 3L), size$p, 3L)
  set.seed(5L)
  sink_host <- exp(matrix(rnorm(size$sink_n * size$sink_n), size$sink_n, size$sink_n))
  n_dist <- min(size$n, 512L)
  Xs <- X[seq_len(n_dist), , drop = FALSE]

  if (identical(op, "sinkhorn")) {
    probe_x <- make_dense_operand(sink_host, requested_backend)
    dispatch <- compute_dispatch_info(probe_x, "matmul", requested_backend = requested_backend)
    release_residency(probe_x)
  } else {
    probe_x <- make_dense_operand(if (op %in% c("chol", "solve_rhs", "eigen_sym")) SPD else X, requested_backend)
    probe_y <- switch(
      op,
      matmul = make_dense_operand(B_host, requested_backend),
      solve_rhs = make_dense_operand(SPD_rhs, requested_backend),
      NULL
    )
    dispatch_op <- switch(
      op,
      matmul = "matmul",
      crossprod = "crossprod",
      covariance = "crossprod",
      dist = "crossprod",
      chol = "chol",
      solve_rhs = "solve",
      eigen_sym = "eigen",
      many_lm = "matmul",
      rsvd = "rsvd",
      svd = "svd",
      sinkhorn = "matmul"
    )
    dispatch <- compute_dispatch_info(probe_x, dispatch_op, y = probe_y, requested_backend = requested_backend)
    release_residency(probe_x)
    if (inherits(probe_y, "aMatrix")) {
      release_residency(probe_y)
    }
  }

  result_meta <- list(
    suite = "dense",
    op = op,
    size_label = size_label,
    variant = variant,
    requested_backend = requested_backend,
    dispatch_probe_op = dispatch$dispatch_probe_op,
    requested_supported = dispatch$requested_supported,
    requested_support_reason = dispatch$requested_support_reason,
    dispatch_backend = dispatch$dispatch_backend,
    dispatch_path = dispatch$dispatch_path,
    nrow = if (identical(op, "sinkhorn")) size$sink_n else if (identical(op, "eigen_sym")) size$p else size$n,
    ncol = if (identical(op, "sinkhorn")) size$sink_n else size$p,
    rhs_width = switch(op, matmul = ncol(B_host), solve_rhs = ncol(SPD_rhs), many_lm = ncol(Y), sinkhorn = size$sink_n, 0L)
  )

  if (!identical(requested_backend, "cpu") && identical(dispatch$requested_supported, FALSE)) {
    return(do.call(
      new_result_row,
      c(
        result_meta,
        list(status = "unsupported", error_message = dispatch$requested_support_reason %||% "requested backend unsupported")
      )
    ))
  }

  reps <- 7L
  runner <- switch(
    variant,
    cold = switch(
      op,
      matmul = function() { aX <- make_dense_operand(X, requested_backend); aB <- make_dense_operand(B_host, requested_backend); on.exit({ release_residency(aX); release_residency(aB) }, add = TRUE); invisible(aX %*% aB) },
      crossprod = function() { aX <- make_dense_operand(X, requested_backend); on.exit(release_residency(aX), add = TRUE); invisible(crossprod(aX)) },
      covariance = function() { aX <- make_dense_operand(X, requested_backend); on.exit(release_residency(aX), add = TRUE); invisible(covariance(aX)) },
      dist = function() { aXs <- make_dense_operand(Xs, requested_backend); on.exit(release_residency(aXs), add = TRUE); invisible(dist_matrix(aXs)) },
      chol = function() { aS <- make_dense_operand(SPD, requested_backend); on.exit(release_residency(aS), add = TRUE); invisible(chol(aS)) },
      solve_rhs = function() { aS <- make_dense_operand(SPD, requested_backend); aR <- make_dense_operand(SPD_rhs, requested_backend); on.exit({ release_residency(aS); release_residency(aR) }, add = TRUE); invisible(solve(aS, aR)) },
      eigen_sym = function() { aS <- make_dense_operand(SPD, requested_backend); on.exit(release_residency(aS), add = TRUE); invisible(eigh(aS)) },
      many_lm = function() { aX <- make_dense_operand(X, requested_backend); on.exit(release_residency(aX), add = TRUE); invisible(many_lm(aX, Y, method = "qr", cache = FALSE)) },
      rsvd = function() { aX <- make_dense_operand(X, requested_backend); on.exit(release_residency(aX), add = TRUE); invisible(rsvd(aX, k = 10L)) },
      svd = function() { aX <- make_dense_operand(X, requested_backend); on.exit(release_residency(aX), add = TRUE); invisible(svd(aX)) },
      sinkhorn = function() { aSink <- make_dense_operand(sink_host, requested_backend); on.exit(release_residency(aSink), add = TRUE); invisible(sinkhorn(aSink, max_iter = 25L, tol = 0, return_info = FALSE)) }
    ),
    warm = {
      aX <- if (op %in% c("chol", "solve_rhs", "eigen_sym")) NULL else make_dense_operand(X, requested_backend)
      aB <- if (identical(op, "matmul")) make_dense_operand(B_host, requested_backend) else NULL
      aXs <- if (identical(op, "dist")) make_dense_operand(Xs, requested_backend) else NULL
      aS <- if (op %in% c("chol", "solve_rhs", "eigen_sym")) make_dense_operand(SPD, requested_backend) else NULL
      aR <- if (identical(op, "solve_rhs")) make_dense_operand(SPD_rhs, requested_backend) else NULL
      aSink <- if (identical(op, "sinkhorn")) make_dense_operand(sink_host, requested_backend) else NULL
      on.exit({
        for (obj in list(aX, aB, aXs, aS, aR, aSink)) {
          release_residency(obj)
        }
      }, add = TRUE)

      tryCatch({
        if (!is.null(aB)) invisible(aX %*% aB)
        if (identical(op, "crossprod")) invisible(crossprod(aX))
        if (identical(op, "covariance")) invisible(covariance(aX))
        if (identical(op, "dist")) invisible(dist_matrix(aXs))
        if (identical(op, "chol")) invisible(chol(aS))
        if (identical(op, "solve_rhs")) invisible(solve(aS, aR))
        if (identical(op, "eigen_sym")) invisible(eigh(aS))
        if (identical(op, "many_lm")) invisible(many_lm(aX, Y, method = "qr", cache = TRUE))
        if (identical(op, "rsvd")) invisible(rsvd(aX, k = 10L))
        if (identical(op, "svd")) invisible(svd(aX))
        if (identical(op, "sinkhorn")) invisible(sinkhorn(aSink, max_iter = 5L, tol = 0, return_info = FALSE))
      }, error = function(e) NULL)

      switch(
        op,
        matmul = function() invisible(aX %*% aB),
        crossprod = function() invisible(crossprod(aX)),
        covariance = function() invisible(covariance(aX)),
        dist = function() invisible(dist_matrix(aXs)),
        chol = function() invisible(chol(aS)),
        solve_rhs = function() invisible(solve(aS, aR)),
        eigen_sym = function() invisible(eigh(aS)),
        many_lm = function() invisible(many_lm(aX, Y, method = "qr", cache = TRUE)),
        rsvd = function() invisible(rsvd(aX, k = 10L)),
        svd = function() invisible(svd(aX)),
        sinkhorn = function() invisible(sinkhorn(aSink, max_iter = 25L, tol = 0, return_info = FALSE))
      )
    },
    stop(sprintf("unknown dense variant '%s'", variant), call. = FALSE)
  )

  timing <- benchmark_time_ms(runner, reps = reps)

  rel_err <- NA_real_
  accuracy_error <- NA_character_
  if (!identical(requested_backend, "cpu") && op %in% c("matmul", "crossprod", "rsvd", "chol")) {
    rel_err <- tryCatch({
      ref <- switch(
        op,
        matmul    = as.matrix(adgeMatrix(X, preferred_backend = "cpu") %*% adgeMatrix(B_host, preferred_backend = "cpu")),
        crossprod = as.matrix(crossprod(adgeMatrix(X, preferred_backend = "cpu"))),
        rsvd      = { r <- rsvd(adgeMatrix(X, preferred_backend = "cpu"), k = 10L); r$u %*% diag(r$d) %*% t(r$v) },
        chol      = as.matrix(chol(adgeMatrix(SPD, preferred_backend = "cpu")))
      )
      gpu_result <- switch(
        op,
        matmul    = as.matrix(adgeMatrix(X, preferred_backend = requested_backend) %*% adgeMatrix(B_host, preferred_backend = requested_backend)),
        crossprod = as.matrix(crossprod(adgeMatrix(X, preferred_backend = requested_backend))),
        rsvd      = { r <- rsvd(adgeMatrix(X, preferred_backend = requested_backend), k = 10L); r$u %*% diag(r$d) %*% t(r$v) },
        chol      = as.matrix(chol(adgeMatrix(SPD, preferred_backend = requested_backend)))
      )
      assert_backend_accuracy(ref, gpu_result, op)
    }, error = function(e) {
      accuracy_error <<- conditionMessage(e)
      NA_real_
    })
  }

  new_result_row(
    suite = result_meta$suite,
    op = result_meta$op,
    size_label = result_meta$size_label,
    variant = result_meta$variant,
    requested_backend = result_meta$requested_backend,
    dispatch_probe_op = result_meta$dispatch_probe_op,
    requested_supported = result_meta$requested_supported,
    requested_support_reason = result_meta$requested_support_reason,
    dispatch_backend = result_meta$dispatch_backend,
    dispatch_path = result_meta$dispatch_path,
    status = if (is.na(accuracy_error)) "ok" else "error",
    error_message = accuracy_error,
    nrow = result_meta$nrow,
    ncol = result_meta$ncol,
    rhs_width = result_meta$rhs_width,
    nnz = 0L,
    density = 0,
    density_bucket = "dense",
    reps = reps,
    median_ms = timing$median,
    mean_ms = timing$mean,
    sd_ms = timing$sd,
    p05_ms = timing$p05,
    p95_ms = timing$p95,
    n_reps = timing$n_reps,
    rel_err = rel_err
  )
}

make_sparse_case <- function(nrow, ncol, density, rhs_width, seed) {
  set.seed(seed)
  list(
    X_host = Matrix::rsparsematrix(nrow, ncol, density = density),
    rhs_host = matrix(rnorm(ncol * rhs_width), nrow = ncol, ncol = rhs_width)
  )
}

run_sparse_case <- function(requested_backend, op, size_label, density, rhs_width, variant) {
  size <- sparse_sizes()[[size_label]]
  if (is.null(size)) {
    stop(sprintf("unknown sparse size '%s'", size_label), call. = FALSE)
  }

  rhs_width <- as.integer(rhs_width)
  case <- make_sparse_case(
    nrow = size$nrow,
    ncol = size$ncol,
    density = density,
    rhs_width = rhs_width,
    seed = 1000L + rhs_width + as.integer(density * 1000)
  )

  sparse_result <- with_sparse_benchmark_backend_options(requested_backend, rhs_width, {
    probe_x <- make_sparse_operand(case$X_host, requested_backend)
    probe_y <- make_dense_operand(case$rhs_host, requested_backend)

    if (identical(variant, "resident") &&
        requested_backend %in% c("mlx", "metal", "arrayfire")) {
      probe_x <- amatrix::amatrix_bind_resident(probe_x, backend = requested_backend)
      probe_y <- amatrix::amatrix_bind_resident(probe_y, backend = requested_backend)
    }

    dispatch <- compute_dispatch_info(probe_x, "matmul", y = probe_y, requested_backend = requested_backend)
    release_residency(probe_x)
    release_residency(probe_y)

    if (!identical(requested_backend, "cpu") && identical(dispatch$requested_supported, FALSE)) {
      list(
        dispatch = dispatch,
        reps = NA_integer_,
        timing = NULL,
        unsupported = TRUE
      )
    } else {

    reps <- 5L
    timing <- switch(
      variant,
      cold = {
        runner <- function() {
          x <- make_sparse_operand(case$X_host, requested_backend)
          y <- make_dense_operand(case$rhs_host, requested_backend)
          on.exit({
            release_residency(x)
            release_residency(y)
          }, add = TRUE)
          invisible(x %*% y)
        }
        benchmark_time_ms(runner, reps = reps, warmup = 0L)
      },
      resident = {
        resident_x <- make_sparse_operand(case$X_host, requested_backend)
        resident_y <- make_dense_operand(case$rhs_host, requested_backend)

        if (requested_backend %in% c("mlx", "metal", "arrayfire")) {
          resident_x <- amatrix::amatrix_bind_resident(resident_x, backend = requested_backend)
          resident_y <- amatrix::amatrix_bind_resident(resident_y, backend = requested_backend)
        }

        on.exit({
          release_residency(resident_x)
          release_residency(resident_y)
        }, add = TRUE)

        runner <- function() invisible(resident_x %*% resident_y)
        # Explicit prime: run once outside the timer so the warm variant
        # measures a primed state. benchmark_time_ms is then called with
        # warmup = 0L so every timed rep is genuinely warm.
        runner()
        benchmark_time_ms(runner, reps = reps, warmup = 0L)
      },
      stop(sprintf("unknown sparse variant '%s'", variant), call. = FALSE)
    )

    list(
      dispatch = dispatch,
      reps = reps,
      timing = timing
    )
    } # end else (supported path)
  })

  dispatch <- sparse_result$dispatch
  if (isTRUE(sparse_result$unsupported)) {
    density_value <- amatrix:::.amatrix_sparse_density(case$X_host)
    return(new_result_row(
      suite = "sparse",
      op = op,
      size_label = size_label,
      variant = variant,
      requested_backend = requested_backend,
      dispatch_probe_op = dispatch$dispatch_probe_op,
      requested_supported = dispatch$requested_supported,
      requested_support_reason = dispatch$requested_support_reason,
      dispatch_backend = dispatch$dispatch_backend,
      dispatch_path = dispatch$dispatch_path,
      nrow = size$nrow,
      ncol = size$ncol,
      rhs_width = rhs_width,
      nnz = length(case$X_host@x),
      density = density_value,
      density_bucket = amatrix:::.amatrix_sparse_density_bucket(density_value),
      status = "unsupported",
      error_message = dispatch$requested_support_reason %||% "requested backend unsupported"
    ))
  }
  reps <- sparse_result$reps
  timing <- sparse_result$timing
  density_value <- amatrix:::.amatrix_sparse_density(case$X_host)

  new_result_row(
    suite = "sparse",
    op = op,
    size_label = size_label,
    variant = variant,
    requested_backend = requested_backend,
    dispatch_probe_op = dispatch$dispatch_probe_op,
    requested_supported = dispatch$requested_supported,
    requested_support_reason = dispatch$requested_support_reason,
    dispatch_backend = dispatch$dispatch_backend,
    dispatch_path = dispatch$dispatch_path,
    nrow = size$nrow,
    ncol = size$ncol,
    rhs_width = rhs_width,
    nnz = length(case$X_host@x),
    density = density_value,
    density_bucket = amatrix:::.amatrix_sparse_density_bucket(density_value),
    reps = reps,
    median_ms = timing$median %||% NA_real_,
    mean_ms = timing$mean %||% NA_real_,
    sd_ms = timing$sd %||% NA_real_,
    p05_ms = timing$p05 %||% NA_real_,
    p95_ms = timing$p95 %||% NA_real_,
    n_reps = timing$n_reps %||% NA_integer_
  )
}

run_sparse_iterative_case <- function(requested_backend, op, size_label, density, variant) {
  size <- sparse_sizes()[[size_label]]
  if (is.null(size)) {
    stop(sprintf("unknown sparse size '%s'", size_label), call. = FALSE)
  }

  params <- sparse_iterative_params(op, size_label)
  case <- make_sparse_case(
    nrow = size$nrow,
    ncol = size$ncol,
    density = density,
    rhs_width = 1L,
    seed = 2000L + as.integer(density * 1000)
  )

  iterative_result <- with_sparse_benchmark_backend_options(requested_backend, 32L, {
    probe_x <- make_sparse_operand(case$X_host, requested_backend)
    dispatch_backend <- if (identical(op, "svd_factor_subspace")) {
      plan <- tryCatch(amatrix:::.amatrix_svd_factor_plan(
        probe_x,
        k = params$k,
        method = "subspace",
        n_oversamples = params$n_oversamples %||% 0L,
        n_iter = params$n_iter %||% 0L
      ), error = function(e) NULL)
      plan$subspace_backend %||% requested_backend
    } else {
      tryCatch(amatrix::amatrix_resident_backend_for(probe_x, op = "matmul"), error = function(e) NULL) %||% requested_backend
    }
    release_residency(probe_x)
    requested_supported <- if (identical(requested_backend, "cpu")) {
      NA
    } else {
      isTRUE(!is.na(dispatch_backend) && identical(dispatch_backend, requested_backend))
    }
    requested_support_reason <- if (identical(requested_backend, "cpu")) {
      NA_character_
    } else if (isTRUE(requested_supported)) {
      "iterative supported"
    } else {
      "iterative rerouted"
    }

    reps <- 3L
    timing <- switch(
      variant,
      cold = {
        runner <- function() {
          x <- make_sparse_operand(case$X_host, requested_backend)
          on.exit(release_residency(x), add = TRUE)
          if (identical(op, "block_lanczos")) {
            invisible(block_lanczos(
              x,
              nv = params$k,
              nu = params$k,
              block_size = params$block_size,
              n_steps = params$n_steps
            ))
          } else {
            invisible(svd_factor(
              x,
              k = params$k,
              method = "subspace",
              n_oversamples = params$n_oversamples,
              n_iter = params$n_iter
            ))
          }
        }
        benchmark_time_ms(runner, reps = reps, warmup = 0L)
      },
      warm = {
        resident_x <- make_sparse_operand(case$X_host, requested_backend)
        if (requested_backend %in% c("mlx", "metal", "arrayfire")) {
          resident_x <- amatrix::amatrix_bind_resident(resident_x, backend = requested_backend, op = "matmul")
        }
        on.exit(release_residency(resident_x), add = TRUE)

        runner <- if (identical(op, "block_lanczos")) {
          function() invisible(block_lanczos(
            resident_x,
            nv = params$k,
            nu = params$k,
            block_size = params$block_size,
            n_steps = params$n_steps
          ))
        } else {
          function() {
            drop_svd_factor_cache(
              resident_x,
              k = params$k,
              n_oversamples = params$n_oversamples,
              n_iter = params$n_iter
            )
            invisible(svd_factor(
              resident_x,
              k = params$k,
              method = "subspace",
              n_oversamples = params$n_oversamples,
              n_iter = params$n_iter
            ))
          }
        }
        # Explicit prime: execute once outside the timer so every timed rep
        # runs against a primed resident state. benchmark_time_ms is called
        # with warmup = 0L so the prime is not double-counted.
        runner()
        benchmark_time_ms(runner, reps = reps, warmup = 0L)
      },
      stop(sprintf("unknown sparse iterative variant '%s'", variant), call. = FALSE)
    )

    list(
      requested_supported = requested_supported,
      requested_support_reason = requested_support_reason,
      dispatch_backend = dispatch_backend,
      reps = reps,
      timing = timing
    )
  })

  density_value <- amatrix:::.amatrix_sparse_density(case$X_host)

  new_result_row(
    suite = "sparse",
    op = op,
    size_label = size_label,
    variant = variant,
    requested_backend = requested_backend,
    dispatch_probe_op = op,
    requested_supported = iterative_result$requested_supported,
    requested_support_reason = iterative_result$requested_support_reason,
    dispatch_backend = iterative_result$dispatch_backend,
    dispatch_path = "iterative",
    nrow = size$nrow,
    ncol = size$ncol,
    rhs_width = params$k,
    nnz = length(case$X_host@x),
    density = density_value,
    density_bucket = amatrix:::.amatrix_sparse_density_bucket(density_value),
    reps = iterative_result$reps,
    median_ms = iterative_result$timing$median %||% NA_real_,
    mean_ms = iterative_result$timing$mean %||% NA_real_,
    sd_ms = iterative_result$timing$sd %||% NA_real_,
    p05_ms = iterative_result$timing$p05 %||% NA_real_,
    p95_ms = iterative_result$timing$p95 %||% NA_real_,
    n_reps = iterative_result$timing$n_reps %||% NA_integer_
  )
}

safe_case <- function(fun, cell_info) {
  tryCatch(
    fun(),
    error = function(e) {
      new_result_row(
        suite = cell_info$suite,
        op = cell_info$op,
        size_label = cell_info$size_label,
        variant = cell_info$variant,
        requested_backend = cell_info$requested_backend,
        nrow = cell_info$nrow %||% NA_integer_,
        ncol = cell_info$ncol %||% NA_integer_,
        rhs_width = cell_info$rhs_width %||% NA_integer_,
        density = cell_info$density %||% NA_real_,
        status = "error",
        error_message = conditionMessage(e)
      )
    }
  )
}

format_backend_diagnostics <- function(backend_name, include_arrayfire = .benchmark_arrayfire_requested()) {
  parts <- c(sprintf("requested_backend=%s", backend_name))

  spec <- .benchmark_optional_backend_specs(include_arrayfire = include_arrayfire)[[backend_name]]
  if (!is.null(spec)) {
    if (!is.null(spec$env)) {
      env_parts <- vapply(
        names(spec$env),
        function(name) sprintf("%s=%s", name, Sys.getenv(name, unset = "<unset>")),
        character(1)
      )
      parts <- c(parts, env_parts)
    }
    if (!is.null(spec$options)) {
      option_parts <- vapply(
        names(spec$options),
        function(name) sprintf("%s=%s", name, as.character(getOption(name))),
        character(1)
      )
      parts <- c(parts, option_parts)
    }
  }

  backend <- tryCatch(amatrix:::.amatrix_get_backend(backend_name), error = function(e) NULL)
  parts <- c(parts, sprintf("registered=%s", !is.null(backend)))
  if (!is.null(backend)) {
    available <- tryCatch(isTRUE(backend$available()), error = function(e) NA)
    parts <- c(parts, sprintf("backend_available=%s", as.character(available)))
  }

  if (identical(backend_name, "opencl")) {
    ns <- tryCatch(
      ensure_optional_backend_namespace("amatrix.opencl", repo_dir = "backends/amatrix.opencl"),
      error = function(e) NULL
    )
    if (!is.null(ns)) {
      info <- tryCatch(get("amatrix_opencl_bridge_info", envir = ns, inherits = FALSE)(), error = function(e) NULL)
      if (is.list(info)) {
        bridge_keys <- intersect(c("compiled", "clblast", "native", "probe_enabled", "engine", "available"), names(info))
        bridge_parts <- vapply(
          bridge_keys,
          function(name) sprintf("bridge_%s=%s", name, as.character(info[[name]])),
          character(1)
        )
        parts <- c(parts, bridge_parts)
      }
    }
  }

  paste(parts, collapse = ", ")
}

group_plan <- function(backends, suites = c("dense", "sparse")) {
  out <- list()

  if ("dense" %in% suites) {
    dense_ops <- c("matmul", "crossprod", "covariance", "dist", "chol", "solve_rhs", "eigen_sym", "many_lm", "rsvd", "svd", "sinkhorn")
    for (backend in backends) {
      for (op in dense_ops) {
        out[[length(out) + 1L]] <- list(group_id = sprintf("dense-%s-%s", backend, op), suite = "dense", requested_backend = backend, op = op)
      }
    }
  }

  if ("sparse" %in% suites) {
    sparse_ops <- c("spmv", "spmm", sparse_iterative_ops())
    for (backend in backends) {
      for (op in sparse_ops) {
        out[[length(out) + 1L]] <- list(group_id = sprintf("sparse-%s-%s", backend, op), suite = "sparse", requested_backend = backend, op = op)
      }
    }
  }

  out
}

expand_group_rows <- function(group) {
  rows <- list()
  if (identical(group$suite, "dense")) {
    for (size_label in names(dense_sizes())) {
      for (variant in c("cold", "warm")) {
        rows[[length(rows) + 1L]] <- new_result_row(
          suite = "dense",
          op = group$op,
          size_label = size_label,
          variant = variant,
          requested_backend = group$requested_backend,
          nrow = dense_sizes()[[size_label]]$n,
          ncol = dense_sizes()[[size_label]]$p
        )
      }
    }
  }

  if (identical(group$suite, "sparse")) {
    if (group$op %in% sparse_iterative_ops()) {
      for (size_label in names(sparse_sizes())) {
        for (density in sparse_iterative_densities()) {
          for (variant in c("cold", "warm")) {
            params <- sparse_iterative_params(group$op, size_label)
            rows[[length(rows) + 1L]] <- new_result_row(
              suite = "sparse",
              op = group$op,
              size_label = size_label,
              variant = variant,
              requested_backend = group$requested_backend,
              nrow = sparse_sizes()[[size_label]]$nrow,
              ncol = sparse_sizes()[[size_label]]$ncol,
              rhs_width = params$k,
              density = density
            )
          }
        }
      }
    } else {
      rhs_widths <- if (identical(group$op, "spmv")) 1L else spmm_rhs_widths()
      for (size_label in names(sparse_sizes())) {
        for (density in sparse_densities()) {
          for (rhs_width in rhs_widths) {
            for (variant in c("cold", "resident")) {
              rows[[length(rows) + 1L]] <- new_result_row(
                suite = "sparse",
                op = group$op,
                size_label = size_label,
                variant = variant,
                requested_backend = group$requested_backend,
                nrow = sparse_sizes()[[size_label]]$nrow,
                ncol = sparse_sizes()[[size_label]]$ncol,
                rhs_width = rhs_width,
                density = density
              )
            }
          }
        }
      }
    }
  }

  do.call(rbind, rows)
}

run_group <- function(group) {
  backend_name <- group$requested_backend
  if (!benchmark_worker_backend_allowed(backend_name)) {
    rows <- expand_group_rows(group)
    rows$status <- "unavailable"
    rows$error_message <- benchmark_worker_backend_block_reason(backend_name)
    return(rows)
  }
  backend <- tryCatch(amatrix:::.amatrix_get_backend(backend_name), error = function(e) NULL)
  if (is.null(backend) || !isTRUE(backend$available())) {
    rows <- expand_group_rows(group)
    rows$status <- "unavailable"
    rows$error_message <- sprintf(
      "backend '%s' is not available [%s]",
      backend_name,
      format_backend_diagnostics(backend_name)
    )
    return(rows)
  }

  rows <- list()
  if (identical(group$suite, "dense")) {
    for (size_label in names(dense_sizes())) {
      for (variant in c("cold", "warm")) {
        cell <- list(suite = "dense", op = group$op, size_label = size_label, variant = variant, requested_backend = backend_name)
        rows[[length(rows) + 1L]] <- safe_case(
          function() run_dense_case(backend_name, group$op, size_label, variant),
          cell
        )
      }
    }
  }

  if (identical(group$suite, "sparse")) {
    if (group$op %in% sparse_iterative_ops()) {
      for (size_label in names(sparse_sizes())) {
        for (density in sparse_iterative_densities()) {
          for (variant in c("cold", "warm")) {
            params <- sparse_iterative_params(group$op, size_label)
            cell <- list(
              suite = "sparse",
              op = group$op,
              size_label = size_label,
              variant = variant,
              requested_backend = backend_name,
              rhs_width = params$k,
              density = density,
              nrow = sparse_sizes()[[size_label]]$nrow,
              ncol = sparse_sizes()[[size_label]]$ncol
            )
            rows[[length(rows) + 1L]] <- safe_case(
              function() run_sparse_iterative_case(backend_name, group$op, size_label, density, variant),
              cell
            )
          }
        }
      }
    } else {
      rhs_widths <- if (identical(group$op, "spmv")) 1L else spmm_rhs_widths()
      for (size_label in names(sparse_sizes())) {
        for (density in sparse_densities()) {
          for (rhs_width in rhs_widths) {
            for (variant in c("cold", "resident")) {
              cell <- list(
                suite = "sparse",
                op = group$op,
                size_label = size_label,
                variant = variant,
                requested_backend = backend_name,
                rhs_width = rhs_width,
                density = density,
                nrow = sparse_sizes()[[size_label]]$nrow,
                ncol = sparse_sizes()[[size_label]]$ncol
              )
              rows[[length(rows) + 1L]] <- safe_case(
                function() run_sparse_case(backend_name, group$op, size_label, density, rhs_width, variant),
                cell
              )
            }
          }
        }
      }
    }
  }

  do.call(rbind, rows)
}

key_columns <- c(
  "suite", "op", "size_label", "variant", "requested_backend",
  "dispatch_backend", "dispatch_path", "nrow", "ncol", "rhs_width",
  "density", "density_bucket"
)

baseline_min_key <- c("op", "size_label", "variant", "requested_backend")

add_cpu_reference <- function(results) {
  cpu <- results[
    results$status == "ok" &
      results$requested_backend == "cpu",
    c("suite", "op", "size_label", "variant", "nrow", "ncol", "rhs_width", "density", "density_bucket", "median_ms")
  ]
  names(cpu)[names(cpu) == "median_ms"] <- "cpu_reference_ms"

  merge(
    results,
    cpu,
    by = c("suite", "op", "size_label", "variant", "nrow", "ncol", "rhs_width", "density", "density_bucket"),
    all.x = TRUE
  )
}

normalize_baseline_schema <- function(baseline) {
  if ("size" %in% names(baseline) && !"size_label" %in% names(baseline)) {
    names(baseline)[names(baseline) == "size"] <- "size_label"
  }
  if ("backend" %in% names(baseline) && !"requested_backend" %in% names(baseline)) {
    names(baseline)[names(baseline) == "backend"] <- "requested_backend"
  }

  legacy_to_canonical <- c(
    "256x32"    = "small",
    "1024x128"  = "medium",
    "4096x128"  = "large",
    "4096x1024" = "xlarge"
  )
  if ("size_label" %in% names(baseline)) {
    hits <- baseline$size_label %in% names(legacy_to_canonical)
    if (any(hits)) {
      baseline$size_label[hits] <- unname(legacy_to_canonical[baseline$size_label[hits]])
    }
  }
  baseline
}

add_baseline_compare <- function(results, baseline_path) {
  results$baseline_ms <- NA_real_
  results$ratio_vs_baseline <- NA_real_

  if (!file.exists(baseline_path)) {
    return(results)
  }

  baseline <- tryCatch(read.csv(baseline_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(baseline) || !"median_ms" %in% names(baseline)) {
    return(results)
  }

  baseline <- normalize_baseline_schema(baseline)

  if ("status" %in% names(baseline)) {
    baseline <- baseline[baseline$status == "ok", , drop = FALSE]
  }

  if (!all(baseline_min_key %in% names(baseline))) {
    return(results)
  }

  join_key <- intersect(key_columns, names(baseline))
  if (!all(baseline_min_key %in% join_key)) {
    join_key <- baseline_min_key
  }

  baseline <- baseline[, c(join_key, "median_ms"), drop = FALSE]
  names(baseline)[names(baseline) == "median_ms"] <- "baseline_ms"
  baseline <- baseline[!duplicated(baseline[, join_key, drop = FALSE]), , drop = FALSE]

  results$baseline_ms <- NULL
  results$ratio_vs_baseline <- NULL
  merged <- merge(results, baseline, by = join_key, all.x = TRUE, sort = FALSE)
  merged$ratio_vs_baseline <- ifelse(
    is.na(merged$baseline_ms) | is.na(merged$median_ms),
    NA_real_,
    merged$median_ms / merged$baseline_ms
  )
  merged
}

summarize_results <- function(results) {
  display <- results[order(results$suite, results$op, results$size_label, results$variant, results$requested_backend), ]
  display$speedup_vs_cpu <- ifelse(
    is.na(display$cpu_reference_ms) | is.na(display$median_ms),
    NA_real_,
    display$cpu_reference_ms / display$median_ms
  )
  display$selected_backend <- display$dispatch_backend
  display$availability_state <- vapply(seq_len(nrow(display)), function(idx) {
    requested_backend <- display$requested_backend[[idx]]
    support_reason <- display$requested_support_reason[[idx]]
    status <- display$status[[idx]]

    if (identical(requested_backend, "cpu")) {
      return("cpu_baseline")
    }
    if (identical(status, "unavailable") || identical(support_reason, "backend unavailable")) {
      return("unavailable")
    }
    if (identical(support_reason, "backend not registered")) {
      return("unregistered")
    }
    "available"
  }, character(1))
  display$support_state <- vapply(seq_len(nrow(display)), function(idx) {
    requested_backend <- display$requested_backend[[idx]]
    support_reason <- display$requested_support_reason[[idx]]
    requested_supported <- display$requested_supported[[idx]]

    if (identical(requested_backend, "cpu")) {
      return("cpu_baseline")
    }
    if (identical(support_reason, "backend unavailable")) {
      return("backend_unavailable")
    }
    if (identical(support_reason, "backend not registered")) {
      return("backend_unregistered")
    }
    if (identical(support_reason, "precision incompatible")) {
      return("precision_incompatible")
    }
    if (identical(support_reason, "dispatch plan unavailable")) {
      return("unknown")
    }
    if (isTRUE(requested_supported)) {
      return("supported")
    }
    if (identical(support_reason, "calibration rejected")) {
      return("calibration_rejected")
    }
    if (identical(support_reason, "op unsupported") || identical(support_reason, "resident op unsupported")) {
      return("unsupported")
    }
    if (isFALSE(requested_supported)) {
      return("unsupported")
    }
    "unknown"
  }, character(1))
  display$execution_state <- ifelse(display$status == "ok", "executed", display$status)
  display$routing_state <- vapply(seq_len(nrow(display)), function(idx) {
    status <- display$status[[idx]]
    requested_backend <- display$requested_backend[[idx]]
    dispatch_backend <- display$dispatch_backend[[idx]]

    if (!identical(status, "ok")) {
      return("not_run")
    }
    if (identical(requested_backend, "cpu")) {
      return("cpu_baseline")
    }
    if (is.na(dispatch_backend)) {
      return("unknown_dispatch")
    }
    if (identical(dispatch_backend, requested_backend)) {
      return("accelerated")
    }
    if (identical(dispatch_backend, "cpu")) {
      return("cpu_fallback")
    }
    "rerouted"
  }, character(1))
  display$performance_state <- vapply(seq_len(nrow(display)), function(idx) {
    status <- display$status[[idx]]
    requested_backend <- display$requested_backend[[idx]]
    dispatch_backend <- display$dispatch_backend[[idx]]
    speedup_vs_cpu <- display$speedup_vs_cpu[[idx]]

    if (!identical(status, "ok")) {
      return("not_run")
    }
    if (identical(requested_backend, "cpu")) {
      return("cpu_baseline")
    }
    if (is.na(dispatch_backend) || !identical(dispatch_backend, requested_backend)) {
      return("not_accelerated")
    }
    if (is.na(speedup_vs_cpu)) {
      return("accelerated_unknown_speed")
    }
    if (speedup_vs_cpu > 1 + 1e-8) {
      return("accelerated_faster_than_cpu")
    }
    if (abs(speedup_vs_cpu - 1) <= 1e-8) {
      return("accelerated_at_cpu")
    }
    "accelerated_slower_than_cpu"
  }, character(1))
  display$dispatch_state <- display$routing_state
  display
}

benchmark_size_levels <- function() c("small", "medium", "large", "xlarge")

benchmark_ordered_size_label <- function(x) {
  factor(x, levels = benchmark_size_levels(), ordered = TRUE)
}

benchmark_report_regressions <- function(summary) {
  rows <- summary[
    summary$status == "ok" &
      summary$dispatch_state != "cpu_fallback" &
      !is.na(summary$ratio_vs_baseline) &
      summary$ratio_vs_baseline > 1.2,
    c(
      "suite", "op", "size_label", "variant", "requested_backend",
      "dispatch_backend", "median_ms", "baseline_ms", "ratio_vs_baseline"
    ),
    drop = FALSE
  ]
  if (nrow(rows) == 0L) {
    return(rows)
  }
  rows$regression_pct <- (rows$ratio_vs_baseline - 1) * 100
  rows[order(-rows$ratio_vs_baseline, rows$suite, rows$op, rows$size_label, rows$variant), , drop = FALSE]
}

benchmark_report_incidents <- function(summary) {
  summary[
    summary$status != "ok",
    c(
      "suite", "op", "size_label", "variant", "requested_backend",
      "status", "error_message", "requested_support_reason"
    ),
    drop = FALSE
  ]
}

benchmark_report_fallbacks <- function(summary) {
  summary[
    summary$status == "ok" &
      summary$requested_backend != "cpu" &
      summary$dispatch_state == "cpu_fallback",
    c(
      "suite", "op", "size_label", "variant", "requested_backend",
      "dispatch_backend", "median_ms", "cpu_reference_ms", "requested_support_reason"
    ),
    drop = FALSE
  ]
}

benchmark_report_warm_pairs <- function(summary) {
  cold <- summary[
    summary$status == "ok" & summary$variant == "cold",
    c(
      "suite", "op", "size_label", "requested_backend", "nrow", "ncol",
      "rhs_width", "density", "density_bucket", "median_ms", "sd_ms"
    ),
    drop = FALSE
  ]
  warm <- summary[
    summary$status == "ok" & summary$variant %in% c("warm", "resident"),
    c(
      "suite", "op", "size_label", "requested_backend", "nrow", "ncol",
      "rhs_width", "density", "density_bucket", "variant", "median_ms", "sd_ms"
    ),
    drop = FALSE
  ]
  if (nrow(cold) == 0L || nrow(warm) == 0L) {
    return(data.frame())
  }

  names(cold)[names(cold) == "median_ms"] <- "cold_ms"
  names(cold)[names(cold) == "sd_ms"] <- "cold_sd_ms"
  names(warm)[names(warm) == "variant"] <- "warm_variant"
  names(warm)[names(warm) == "median_ms"] <- "warm_ms"
  names(warm)[names(warm) == "sd_ms"] <- "warm_sd_ms"

  join_cols <- c("suite", "op", "size_label", "requested_backend", "nrow", "ncol", "rhs_width", "density", "density_bucket")
  pairs <- merge(cold, warm, by = join_cols, all = FALSE, sort = FALSE)
  if (nrow(pairs) == 0L) {
    return(pairs)
  }

  pairs$cold_to_warm_gain <- pairs$cold_ms / pairs$warm_ms
  pairs$warm_vs_cold_ratio <- pairs$warm_ms / pairs$cold_ms
  pairs$ms_saved <- pairs$cold_ms - pairs$warm_ms
  pairs[order(-pairs$cold_to_warm_gain, pairs$suite, pairs$op, pairs$size_label), , drop = FALSE]
}

benchmark_report_acceleration <- function(summary) {
  rows <- summary[
    summary$status == "ok" &
      summary$requested_backend != "cpu" &
      summary$dispatch_state == "accelerated",
    c(
      "suite", "op", "size_label", "variant", "requested_backend",
      "dispatch_backend", "median_ms", "cpu_reference_ms",
      "speedup_vs_cpu", "ratio_vs_baseline", "rel_err"
    ),
    drop = FALSE
  ]
  if (nrow(rows) == 0L) {
    return(rows)
  }
  rows[order(-rows$speedup_vs_cpu, rows$suite, rows$op, rows$size_label), , drop = FALSE]
}

benchmark_report_backend_overview <- function(summary, regressions) {
  backend_names <- unique(summary$requested_backend)
  rows <- lapply(backend_names, function(backend_name) {
    subset <- summary[summary$requested_backend == backend_name, , drop = FALSE]
    data.frame(
      requested_backend = backend_name,
      total_cells = nrow(subset),
      executed_cells = sum(subset$status == "ok"),
      accelerated_cells = sum(subset$dispatch_state == "accelerated", na.rm = TRUE),
      cpu_fallback_cells = sum(subset$dispatch_state == "cpu_fallback", na.rm = TRUE),
      incident_cells = sum(subset$status != "ok"),
      regression_cells = sum(regressions$requested_backend == backend_name),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$requested_backend), , drop = FALSE]
}

benchmark_report_suite_overview <- function(summary) {
  suites <- unique(summary$suite)
  rows <- lapply(suites, function(suite_name) {
    subset <- summary[summary$suite == suite_name, , drop = FALSE]
    data.frame(
      suite = suite_name,
      ops = length(unique(subset$op)),
      requested_backends = paste(sort(unique(subset$requested_backend)), collapse = ", "),
      total_cells = nrow(subset),
      executed_cells = sum(subset$status == "ok"),
      non_cpu_executed = sum(subset$status == "ok" & !is.na(subset$dispatch_backend) & subset$dispatch_backend != "cpu"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$suite), , drop = FALSE]
}

benchmark_report_op_coverage <- function(summary) {
  keys <- unique(summary[, c("suite", "op"), drop = FALSE])
  rows <- lapply(seq_len(nrow(keys)), function(idx) {
    subset <- summary[
      summary$suite == keys$suite[[idx]] &
        summary$op == keys$op[[idx]],
      ,
      drop = FALSE
    ]
    executed_dispatch <- sort(unique(subset$dispatch_backend[subset$status == "ok" & nzchar(subset$dispatch_backend)]))
    data.frame(
      suite = keys$suite[[idx]],
      op = keys$op[[idx]],
      requested_backends = paste(sort(unique(subset$requested_backend)), collapse = ", "),
      executed_backends = if (length(executed_dispatch) > 0L) paste(executed_dispatch, collapse = ", ") else "none",
      variants = paste(sort(unique(subset$variant)), collapse = ", "),
      cells = nrow(subset),
      executed_cells = sum(subset$status == "ok"),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$suite, out$op), , drop = FALSE]
}

benchmark_report_backend_status <- function() {
  status <- tryCatch(amatrix_backend_status(), error = function(e) NULL)
  if (is.null(status) || !is.data.frame(status) || nrow(status) == 0L) {
    return(data.frame())
  }

  keep <- intersect(
    c("name", "available", "health", "precision_modes", "residency_capable", "capabilities"),
    names(status)
  )
  status <- status[, keep, drop = FALSE]
  if ("capabilities" %in% names(status)) {
    status$capabilities <- vapply(
      status$capabilities,
      function(x) {
        x <- as.character(x)
        if (length(x) == 0L || is.na(x)) {
          return(NA_character_)
        }
        if (nchar(x) > 72L) {
          paste0(substr(x, 1L, 69L), "...")
        } else {
          x
        }
      },
      character(1)
    )
  }
  status[order(status$name), , drop = FALSE]
}

benchmark_report_routing_summary <- function(summary) {
  if (nrow(summary) == 0L) {
    return(data.frame())
  }

  keys <- unique(summary[, c("requested_backend", "routing_state", "execution_state"), drop = FALSE])
  counts <- integer(nrow(keys))
  for (idx in seq_len(nrow(keys))) {
    counts[[idx]] <- sum(
      summary$requested_backend == keys$requested_backend[[idx]] &
        summary$routing_state == keys$routing_state[[idx]] &
        summary$execution_state == keys$execution_state[[idx]]
    )
  }
  keys$cells <- counts
  keys[order(keys$requested_backend, keys$routing_state, keys$execution_state), , drop = FALSE]
}

benchmark_report_metrics <- function(summary, regressions, incidents, fallbacks, warm_pairs, baseline_path) {
  data.frame(
    metric = c(
      "total_cells",
      "executed_cells",
      "baseline_matched_cells",
      "regressions_gt_20pct",
      "incidents",
      "cpu_fallbacks",
      "warm_pairs",
      "requested_backends"
    ),
    value = c(
      nrow(summary),
      sum(summary$status == "ok"),
      sum(summary$status == "ok" & !is.na(summary$baseline_ms)),
      nrow(regressions),
      nrow(incidents),
      nrow(fallbacks),
      nrow(warm_pairs),
      length(unique(summary$requested_backend))
    ),
    detail = c(
      "Total benchmark cells in summary.csv",
      "Rows that completed successfully",
      if (file.exists(baseline_path)) "Executed rows with a baseline comparison" else "Baseline file not present",
      "Rows slower than baseline by more than 20%",
      "Rows with status != ok",
      "Successful rows that routed back to CPU",
      "Matched cold vs warm/resident pairs",
      paste(unique(summary$requested_backend), collapse = ", ")
    ),
    stringsAsFactors = FALSE
  )
}

benchmark_report_requires_cli_entry <- function(summary, metadata) {
  benchmark_metadata_on_apple_silicon(metadata) &&
    "arrayfire" %in% unique(summary$requested_backend)
}

benchmark_report_policy_notes <- function(summary, metadata) {
  rows <- list(
    data.frame(
      area = "Scope",
      policy = "Bundle-scoped results only",
      detail = "This report only covers the suites, operations, variants, and backends requested for this bundle. Missing algorithms or backends were not benchmarked here.",
      stringsAsFactors = FALSE
    ),
    data.frame(
      area = "Algorithmic scripts",
      policy = "Standalone scripts are separate from the canonical bundle",
      detail = paste(
        "Standalone exploratory scripts such as tools/benchmark-kmeans.R,",
        "tools/benchmark-cancor.R, and tools/benchmark-krr.R are not part of",
        "summary.csv unless their operations are promoted into tools/benchmark-regression.R."
      ),
      stringsAsFactors = FALSE
    )
  )

  if (benchmark_report_requires_cli_entry(summary, metadata)) {
    rows[[length(rows) + 1L]] <- data.frame(
      area = "Apple ArrayFire",
      policy = "Diagnostic-only on Apple Silicon",
      detail = paste(
        "For Apple Silicon bundles that include ArrayFire, prefer cpu/mlx/metal",
        "for canonical regression runs. ArrayFire requires explicit opt-in and",
        "the direct Rscript launcher; the shell wrapper is not the authoritative",
        "entry path for this case."
      ),
      stringsAsFactors = FALSE
    )
  } else if ("arrayfire" %in% unique(summary$requested_backend)) {
    rows[[length(rows) + 1L]] <- data.frame(
      area = "ArrayFire",
      policy = "Opt-in backend",
      detail = "ArrayFire only appears in the canonical harness when it is explicitly requested during backend discovery.",
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

benchmark_report_snippets <- function(summary, output_dir, baseline_path, metadata) {
  backend_arg <- paste(unique(summary$requested_backend), collapse = ",")
  suite_arg <- paste(unique(summary$suite), collapse = ",")
  summary_cols <- paste(names(summary), collapse = ", ")
  run_entry <- if (benchmark_report_requires_cli_entry(summary, metadata)) {
    "Rscript tools/benchmark-regression-cli.R"
  } else {
    "bash tools/benchmark-regression.sh"
  }
  arrayfire_debug <- if (benchmark_report_requires_cli_entry(summary, metadata)) {
    paste(
      "Rscript tools/debug-arrayfire-gpu.R --runtime=cpu",
      "Rscript tools/debug-arrayfire-gpu.R --runtime=opencl",
      sep = "\n"
    )
  } else {
    NULL
  }

  list(
    run_bundle = paste(
      run_entry,
      sprintf("--backends=%s", backend_arg),
      sprintf("--suites=%s", suite_arg),
      sprintf("--output-dir=%s", output_dir)
    ),
    render_html = paste(
      sprintf("cd %s", output_dir),
      "quarto render benchmark-report.qmd --to html --output benchmark-report.html",
      sep = "\n"
    ),
    render_pdf = paste(
      sprintf("cd %s", output_dir),
      "quarto render benchmark-report.qmd --to pdf --output benchmark-report.pdf",
      sep = "\n"
    ),
    read_bundle_r = paste(
      "report <- readRDS(\"report-data.rds\")",
      "summary <- read.csv(\"summary.csv\", stringsAsFactors = FALSE)",
      "str(report$tables$backend_overview)",
      "summary[c(\"suite\", \"op\", \"size_label\", \"variant\", \"median_ms\")]",
      sep = "\n"
    ),
    summary_schema = summary_cols,
    baseline_path = baseline_path,
    arrayfire_debug = arrayfire_debug
  )
}

benchmark_report_bundle <- function(summary, baseline_path, output_dir, metadata) {
  regressions <- benchmark_report_regressions(summary)
  incidents <- benchmark_report_incidents(summary)
  fallbacks <- benchmark_report_fallbacks(summary)
  warm_pairs <- benchmark_report_warm_pairs(summary)
  acceleration <- benchmark_report_acceleration(summary)
  backend_overview <- benchmark_report_backend_overview(summary, regressions)
  suite_overview <- benchmark_report_suite_overview(summary)
  op_coverage <- benchmark_report_op_coverage(summary)
  backend_status <- benchmark_report_backend_status()
  routing_summary <- benchmark_report_routing_summary(summary)
  metrics <- benchmark_report_metrics(summary, regressions, incidents, fallbacks, warm_pairs, baseline_path)
  policy_notes <- benchmark_report_policy_notes(summary, metadata)

  list(
    summary = summary,
    tables = list(
      metrics = metrics,
      policy_notes = policy_notes,
      backend_overview = backend_overview,
      suite_overview = suite_overview,
      op_coverage = op_coverage,
      backend_status = backend_status,
      regressions = regressions,
      warm_pairs = warm_pairs,
      acceleration = acceleration,
      incidents = incidents,
      fallbacks = fallbacks,
      routing_summary = routing_summary
    ),
    metadata = metadata,
    snippets = benchmark_report_snippets(summary, output_dir, baseline_path, metadata),
    baseline_path = baseline_path,
    baseline_present = file.exists(baseline_path),
    output_dir = output_dir
  )
}

benchmark_report_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15, colour = "#12263A"),
      plot.subtitle = ggplot2::element_text(colour = "#40566B"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold", colour = "#12263A"),
      axis.text = ggplot2::element_text(colour = "#23384D")
    )
}

benchmark_plot_colours <- function() {
  c(
    accelerated = "#00798C",
    cpu_baseline = "#4F6D7A",
    cpu_fallback = "#D1495B",
    not_run = "#B8C5D6",
    rerouted = "#EDAe49",
    available = "#3BB273",
    unavailable = "#8D99AE",
    dense = "#2563EB",
    sparse = "#F59E0B"
  )
}

benchmark_write_placeholder_plot <- function(path, title, subtitle) {
  grDevices::png(path, width = 1600, height = 900, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mar = c(0, 0, 0, 0), bg = "#F8FAFC")
  graphics::plot.new()
  graphics::rect(0, 0, 1, 1, col = "#F8FAFC", border = NA)
  graphics::text(0.5, 0.60, labels = title, cex = 2, font = 2, col = "#12263A")
  graphics::text(0.5, 0.46, labels = subtitle, cex = 1.2, col = "#40566B")
  invisible(path)
}

benchmark_write_plot <- function(plot, path, width = 1800, height = 1000, res = 160) {
  grDevices::png(path, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  invisible(path)
}

benchmark_write_report_plots <- function(bundle, output_dir) {
  plots_dir <- file.path(output_dir, "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  paths <- list(
    baseline = file.path(plots_dir, "baseline-drift.png"),
    warm = file.path(plots_dir, "warm-cold-gains.png"),
    routing = file.path(plots_dir, "routing-overview.png"),
    speedup = file.path(plots_dir, "speedup-vs-cpu.png")
  )

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    benchmark_write_placeholder_plot(paths$baseline, "Baseline drift unavailable", "Install ggplot2 to render benchmark charts.")
    benchmark_write_placeholder_plot(paths$warm, "Warm-state plot unavailable", "Install ggplot2 to render benchmark charts.")
    benchmark_write_placeholder_plot(paths$routing, "Routing plot unavailable", "Install ggplot2 to render benchmark charts.")
    benchmark_write_placeholder_plot(paths$speedup, "Speedup plot unavailable", "Install ggplot2 to render benchmark charts.")
    return(paths)
  }

  cols <- benchmark_plot_colours()
  regressions <- bundle$tables$regressions
  warm_pairs <- bundle$tables$warm_pairs
  routing <- bundle$tables$routing_summary
  acceleration <- bundle$tables$acceleration

  drift_rows <- bundle$summary[
    bundle$summary$status == "ok" & !is.na(bundle$summary$ratio_vs_baseline),
    c("suite", "op", "size_label", "variant", "requested_backend", "ratio_vs_baseline"),
    drop = FALSE
  ]
  if (nrow(drift_rows) > 0L) {
    drift_rows$drift_rank <- abs(log2(drift_rows$ratio_vs_baseline))
    drift_rows <- drift_rows[order(-drift_rows$drift_rank), , drop = FALSE]
    drift_rows <- utils::head(drift_rows, 14L)
    drift_rows$cell <- paste(
      drift_rows$suite,
      drift_rows$op,
      drift_rows$size_label,
      drift_rows$variant,
      drift_rows$requested_backend,
      sep = " / "
    )
    drift_rows$cell <- factor(drift_rows$cell, levels = unique(rev(drift_rows$cell)))
    drift_rows$direction <- ifelse(drift_rows$ratio_vs_baseline > 1, "slower", "faster")
    drift_plot <- ggplot2::ggplot(drift_rows, ggplot2::aes(x = ratio_vs_baseline, y = cell, fill = direction)) +
      ggplot2::geom_col(width = 0.72) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "#12263A") +
      ggplot2::scale_fill_manual(values = c(faster = "#3BB273", slower = "#D1495B")) +
      ggplot2::labs(
        title = "Current runtime versus baseline",
        subtitle = "Top cells by absolute drift; the dashed line is parity with baseline.csv.",
        x = "Current / baseline median runtime",
        y = NULL
      ) +
      benchmark_report_theme()
    benchmark_write_plot(drift_plot, paths$baseline)
  } else {
    benchmark_write_placeholder_plot(paths$baseline, "No baseline comparison", "This run did not match any rows in the baseline file.")
  }

  if (nrow(warm_pairs) > 0L) {
    warm_rows <- utils::head(warm_pairs, 14L)
    warm_rows$cell <- paste(
      warm_rows$suite,
      warm_rows$op,
      warm_rows$size_label,
      warm_rows$requested_backend,
      sep = " / "
    )
    warm_rows$cell <- factor(warm_rows$cell, levels = unique(rev(warm_rows$cell)))
    warm_plot <- ggplot2::ggplot(warm_rows, ggplot2::aes(x = cold_to_warm_gain, y = cell, fill = suite)) +
      ggplot2::geom_col(width = 0.72) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2fx", cold_to_warm_gain)),
        hjust = -0.1,
        size = 3.6,
        colour = "#12263A"
      ) +
      ggplot2::scale_fill_manual(values = c(dense = cols[["dense"]], sparse = cols[["sparse"]])) +
      ggplot2::labs(
        title = "Best warm-state wins",
        subtitle = "Cold-to-warm gain highlights the payoff from reuse and residency.",
        x = "Cold median / warm median",
        y = NULL
      ) +
      ggplot2::expand_limits(x = max(warm_rows$cold_to_warm_gain) * 1.15) +
      benchmark_report_theme()
    benchmark_write_plot(warm_plot, paths$warm)
  } else {
    benchmark_write_placeholder_plot(paths$warm, "No warm pairs", "This run did not produce matching cold and warm/resident cells.")
  }

  if (nrow(routing) > 0L) {
    routing$requested_backend <- factor(routing$requested_backend, levels = unique(routing$requested_backend))
    routing_plot <- ggplot2::ggplot(routing, ggplot2::aes(x = requested_backend, y = cells, fill = routing_state)) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::scale_fill_manual(values = cols[names(cols) %in% unique(routing$routing_state)], drop = FALSE) +
      ggplot2::labs(
        title = "Routing outcomes by requested backend",
        subtitle = "Counts are grouped by the requested backend and the actual routing state.",
        x = NULL,
        y = "Benchmark cells"
      ) +
      benchmark_report_theme()
    benchmark_write_plot(routing_plot, paths$routing)
  } else {
    benchmark_write_placeholder_plot(paths$routing, "No routing rows", "The summary did not contain any routing metadata.")
  }

  if (nrow(acceleration) > 0L) {
    speedup_rows <- acceleration
    speedup_rows$size_label <- benchmark_ordered_size_label(speedup_rows$size_label)
    line_groups <- interaction(speedup_rows$requested_backend, speedup_rows$op, speedup_rows$variant, drop = TRUE)
    speedup_plot <- ggplot2::ggplot(
      speedup_rows,
      ggplot2::aes(x = size_label, y = speedup_vs_cpu, colour = requested_backend, group = interaction(requested_backend, op, variant))
    ) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "#12263A") +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::facet_wrap(~ op, scales = "free_y") +
      ggplot2::labs(
        title = "Observed acceleration versus CPU",
        subtitle = "Only rows that truly executed on the requested non-CPU backend are shown.",
        x = "Size label",
        y = "CPU median / backend median"
      ) +
      benchmark_report_theme()
    if (any(table(line_groups) > 1L)) {
      speedup_plot <- speedup_plot + ggplot2::geom_line(alpha = 0.75)
    }
    benchmark_write_plot(speedup_plot, paths$speedup, width = 1800, height = 1200)
  } else {
    benchmark_write_placeholder_plot(paths$speedup, "No accelerated cells", "Every executed row stayed on CPU in this run.")
  }

  paths
}

benchmark_markdown_table <- function(df, digits = 3L, max_rows = 12L) {
  if (nrow(df) == 0L) {
    return("_None._")
  }

  df <- utils::head(df, max_rows)
  numeric_cols <- vapply(df, is.numeric, logical(1))
  for (name in names(df)[numeric_cols]) {
    df[[name]] <- round(df[[name]], digits = digits)
  }

  if (requireNamespace("knitr", quietly = TRUE)) {
    return(paste(capture.output(knitr::kable(df, format = "pipe")), collapse = "\n"))
  }

  paste0("```\n", paste(capture.output(print(df, row.names = FALSE)), collapse = "\n"), "\n```")
}

benchmark_write_markdown_summary <- function(bundle, output_dir, output_paths) {
  metrics <- bundle$tables$metrics
  metric_value <- function(name) metrics$value[match(name, metrics$metric)][[1L]]
  report_html <- file.path(output_dir, "benchmark-report.html")
  report_pdf <- file.path(output_dir, "benchmark-report.pdf")

  lines <- c(
    "# amatrix benchmark summary",
    "",
    sprintf("- Generated: `%s`", format(bundle$metadata$created_at, "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("- Host: `%s`", bundle$metadata$hostname),
    sprintf("- Baseline: `%s`", if (bundle$baseline_present) bundle$baseline_path else "not found"),
    sprintf("- Requested backends: `%s`", paste(unique(bundle$summary$requested_backend), collapse = ", ")),
    sprintf("- Executed cells: `%s / %s`", metric_value("executed_cells"), metric_value("total_cells")),
    sprintf("- Regressions >20%%: `%s`", metric_value("regressions_gt_20pct")),
    sprintf("- Incidents: `%s`", metric_value("incidents")),
    sprintf("- CPU fallbacks: `%s`", metric_value("cpu_fallbacks")),
    "",
    "## Reproduce and render",
    "",
    "```bash",
    bundle$snippets$run_bundle,
    "```",
    "",
    "```bash",
    bundle$snippets$render_pdf,
    "```",
    "",
    if (!is.null(bundle$snippets$arrayfire_debug)) "### ArrayFire diagnostics" else NULL,
    if (!is.null(bundle$snippets$arrayfire_debug)) "" else NULL,
    if (!is.null(bundle$snippets$arrayfire_debug)) "```bash" else NULL,
    if (!is.null(bundle$snippets$arrayfire_debug)) bundle$snippets$arrayfire_debug else NULL,
    if (!is.null(bundle$snippets$arrayfire_debug)) "```" else NULL,
    if (!is.null(bundle$snippets$arrayfire_debug)) "" else NULL,
    "## Policy notes",
    "",
    benchmark_markdown_table(bundle$tables$policy_notes, digits = 2L, max_rows = 12L),
    "",
    "## Backend overview",
    "",
    benchmark_markdown_table(bundle$tables$backend_overview, digits = 2L, max_rows = 12L),
    "",
    "## Top regressions",
    "",
    benchmark_markdown_table(bundle$tables$regressions, digits = 3L, max_rows = 12L),
    "",
    "## Warm-state gains",
    "",
    benchmark_markdown_table(bundle$tables$warm_pairs, digits = 3L, max_rows = 12L),
    "",
    "## Accelerated cells",
    "",
    benchmark_markdown_table(bundle$tables$acceleration, digits = 3L, max_rows = 12L),
    "",
    "## Incidents",
    "",
    benchmark_markdown_table(bundle$tables$incidents, digits = 3L, max_rows = 12L),
    "",
    "## Artifacts",
    "",
    "- Raw results: `raw-results.csv`",
    "- Decorated summary: `summary.csv`",
    "- Regressions: `regressions.csv`",
    "- Warm ratios: `warm-ratios.csv`",
    "- Routing summary: `routing-summary.csv`",
    "- Quarto source: `benchmark-report.qmd`",
    sprintf("- Rendered HTML: `%s`", if (file.exists(report_html)) basename(report_html) else "not rendered"),
    sprintf("- Rendered PDF: `%s`", if (file.exists(report_pdf)) basename(report_pdf) else "not rendered"),
    "- Charts: `plots/`"
  )

  path <- file.path(output_dir, "benchmark-summary.md")
  writeLines(lines, path)
  path
}

benchmark_render_quarto_report <- function(qmd_path, output_dir) {
  if (identical(Sys.getenv("AMATRIX_BENCHMARK_SKIP_INLINE_RENDER", unset = ""), "1")) {
    return(NA_character_)
  }

  quarto <- Sys.which("quarto")
  if (!nzchar(quarto)) {
    return(NA_character_)
  }

  html_path <- file.path(output_dir, "benchmark-report.html")
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(output_dir)
  capture_fun <- if (exists("benchmark_system2_capture", mode = "function")) {
    get("benchmark_system2_capture", mode = "function")
  } else {
    function(command, args) {
      warned_status <- NULL
      quoted_args <- vapply(args, shQuote, character(1), USE.NAMES = FALSE)
      output <- withCallingHandlers(
        system2(command, quoted_args, stdout = TRUE, stderr = TRUE),
        warning = function(w) {
          warned_status <<- attr(w, "status") %||% warned_status
          invokeRestart("muffleWarning")
        }
      )
      list(
        output = output,
        status = attr(output, "status") %||% warned_status %||% 0L
      )
    }
  }
  render <- capture_fun(
    quarto,
    c("render", basename(qmd_path), "--to", "html", "--output", basename(html_path))
  )
  if (render$status != 0L) {
    message("Quarto render failed: ", paste(render$output, collapse = "\n"))
    return(NA_character_)
  }

  if (!file.exists(html_path)) {
    return(NA_character_)
  }
  html_path
}

benchmark_regression_repo_root <- function() {
  script_path <- getOption("amatrix.benchmark_regression.script_path", "tools/benchmark-regression.R")
  env_root <- Sys.getenv("AMATRIX_BENCHMARK_REPO_ROOT", unset = "")
  candidates <- unique(c(
    env_root,
    dirname(dirname(script_path)),
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", "..")
  ))
  candidates <- candidates[nzchar(candidates)]
  for (candidate in candidates) {
    candidate <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
    if (file.exists(file.path(candidate, "tools", "benchmark-report-template.qmd"))) {
      return(candidate)
    }
  }
  normalizePath(dirname(dirname(script_path)), winslash = "/", mustWork = FALSE)
}

benchmark_prepare_report <- function(summary, output_dir, baseline_path, metadata) {
  bundle <- benchmark_report_bundle(summary, baseline_path, output_dir, metadata)
  plot_paths <- benchmark_write_report_plots(bundle, output_dir)

  bundle$plots <- lapply(plot_paths, basename)

  report_data_path <- file.path(output_dir, "report-data.rds")
  saveRDS(bundle, report_data_path)

  repo_root <- benchmark_regression_repo_root()
  qmd_template <- file.path(repo_root, "tools", "benchmark-report-template.qmd")
  css_template <- file.path(repo_root, "tools", "benchmark-report.css")
  tex_template <- file.path(repo_root, "tools", "benchmark-report-preamble.tex")
  report_qmd <- file.path(output_dir, "benchmark-report.qmd")
  report_css <- file.path(output_dir, "benchmark-report.css")
  report_tex <- file.path(output_dir, "benchmark-report-preamble.tex")
  if (file.exists(qmd_template)) {
    file.copy(qmd_template, report_qmd, overwrite = TRUE)
  }
  if (file.exists(css_template)) {
    file.copy(css_template, report_css, overwrite = TRUE)
  }
  if (file.exists(tex_template)) {
    file.copy(tex_template, report_tex, overwrite = TRUE)
  }

  paths <- list(
    regressions = file.path(output_dir, "regressions.csv"),
    warm_ratios = file.path(output_dir, "warm-ratios.csv"),
    routing_summary = file.path(output_dir, "routing-summary.csv"),
    backend_overview = file.path(output_dir, "backend-overview.csv"),
    report_data = report_data_path,
    report_qmd = report_qmd,
    report_css = report_css,
    report_tex = report_tex
  )

  write.csv(bundle$tables$regressions, paths$regressions, row.names = FALSE)
  write.csv(bundle$tables$warm_pairs, paths$warm_ratios, row.names = FALSE)
  write.csv(bundle$tables$routing_summary, paths$routing_summary, row.names = FALSE)
  write.csv(bundle$tables$backend_overview, paths$backend_overview, row.names = FALSE)

  report_html <- if (file.exists(report_qmd)) benchmark_render_quarto_report(report_qmd, output_dir) else NA_character_
  paths$report_html <- report_html
  paths$summary_md <- benchmark_write_markdown_summary(bundle, output_dir, c(paths, list(report_html = report_html)))
  paths$plots_dir <- file.path(output_dir, "plots")

  paths
}

write_outputs <- function(results, output_dir, baseline_path = "tools/baseline.csv") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  raw_path <- file.path(output_dir, "raw-results.csv")
  summary_path <- file.path(output_dir, "summary.csv")
  incidents_path <- file.path(output_dir, "incidents.csv")
  fallbacks_path <- file.path(output_dir, "fallbacks.csv")
  meta_path <- file.path(output_dir, "metadata.rds")

  summary <- summarize_results(results)
  incidents <- summary[summary$status != "ok", , drop = FALSE]
  fallbacks <- summary[
    summary$status == "ok" &
      summary$requested_backend != "cpu" &
      summary$dispatch_state == "cpu_fallback",
    ,
    drop = FALSE
  ]

  write.csv(results, raw_path, row.names = FALSE)
  write.csv(summary, summary_path, row.names = FALSE)
  write.csv(incidents, incidents_path, row.names = FALSE)
  write.csv(fallbacks, fallbacks_path, row.names = FALSE)
  metadata <- list(
    created_at = Sys.time(),
    hostname = Sys.info()[["nodename"]],
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = R.version$platform,
    sysname = Sys.info()[["sysname"]],
    arch = R.version$arch,
    timezone = Sys.timezone()
  )
  saveRDS(metadata, meta_path)

  report_paths <- benchmark_prepare_report(summary, output_dir, baseline_path, metadata)

  c(
    list(raw = raw_path, summary = summary_path, incidents = incidents_path, fallbacks = fallbacks_path, metadata = meta_path),
    report_paths
  )
}

print_run_summary <- function(results, baseline_path, output_paths, warm_only = FALSE) {
  summary <- summarize_results(results)
  display <- if (warm_only) {
    summary[summary$variant %in% c("warm", "resident"), , drop = FALSE]
  } else {
    summary
  }

  regressions <- display[
    display$status == "ok" &
      display$dispatch_state != "cpu_fallback" &
      !is.na(display$ratio_vs_baseline) &
      display$ratio_vs_baseline > 1.2,
    c("suite", "op", "size_label", "variant", "requested_backend", "median_ms", "baseline_ms", "ratio_vs_baseline")
  ]

  incidents <- display[display$status != "ok", c("suite", "op", "size_label", "variant", "requested_backend", "status", "error_message")]
  fallbacks <- display[
    display$status == "ok" &
      display$requested_backend != "cpu" &
      display$dispatch_state == "cpu_fallback",
    c("suite", "op", "size_label", "variant", "requested_backend", "dispatch_backend", "median_ms", "cpu_reference_ms")
  ]

  if (nrow(regressions) > 0L) {
    message("\nRegressions (>20% slower than baseline):")
    print(regressions, row.names = FALSE)
  } else if (file.exists(baseline_path)) {
    message("\nOK - no regressions vs baseline.")
  } else {
    message("\nNo baseline file found; skipped regression comparison.")
  }

  if (nrow(incidents) > 0L) {
    message("\nIncidents:")
    print(incidents, row.names = FALSE)
  } else {
    message("\nNo incidents recorded.")
  }

  backend_overview <- benchmark_report_backend_overview(display, regressions)
  message("\nBackend overview:")
  print(backend_overview, row.names = FALSE)

  warm_pairs <- benchmark_report_warm_pairs(display)
  if (nrow(warm_pairs) > 0L) {
    message("\nTop warm-state gains:")
    print(
      utils::head(
        warm_pairs[
          ,
          c("suite", "op", "size_label", "requested_backend", "warm_variant", "cold_ms", "warm_ms", "cold_to_warm_gain"),
          drop = FALSE
        ],
        8L
      ),
      row.names = FALSE
    )
  }

  accelerated <- benchmark_report_acceleration(display)
  if (nrow(accelerated) > 0L) {
    message("\nFastest accelerated rows:")
    print(
      utils::head(
        accelerated[
          ,
          c("suite", "op", "size_label", "variant", "requested_backend", "median_ms", "cpu_reference_ms", "speedup_vs_cpu"),
          drop = FALSE
        ],
        8L
      ),
      row.names = FALSE
    )
  } else {
    message("\nNo rows executed on a requested non-CPU backend in this run.")
  }

  if (nrow(fallbacks) > 0L) {
    message("\nCPU fallbacks:")
    print(utils::head(fallbacks, 8L), row.names = FALSE)
  }

  message("\nArtifacts:")
  message("  raw: ", output_paths$raw)
  message("  summary: ", output_paths$summary)
  message("  incidents: ", output_paths$incidents)
  message("  fallbacks: ", output_paths$fallbacks)
  message("  regressions: ", output_paths$regressions)
  message("  warm ratios: ", output_paths$warm_ratios)
  message("  routing summary: ", output_paths$routing_summary)
  message("  markdown summary: ", output_paths$summary_md)
  message("  quarto source: ", output_paths$report_qmd)
  if (!is.na(output_paths$report_html)) {
    message("  html report: ", output_paths$report_html)
  } else {
    message("  html report: not rendered (Quarto unavailable)")
  }
  message("  plots: ", output_paths$plots_dir)
  message("  metadata: ", output_paths$metadata)
}

run_worker <- function(args) {
  if (is.null(args$plan) || is.null(args$out) || is.null(args$group_id)) {
    stop("worker mode requires --plan, --group-id, and --out", call. = FALSE)
  }

  groups <- readRDS(args$plan)
  group <- Filter(function(x) identical(x$group_id, args$group_id), groups)[[1L]]
  prime_ok <- prime_requested_backend(group$requested_backend, include_arrayfire = args$include_arrayfire)
  benchmark_launch_debug(
    "worker group=", group$group_id,
    " ; requested_backend=", group$requested_backend,
    " ; prime_ok=", isTRUE(prime_ok),
    " ; diagnostics=", format_backend_diagnostics(group$requested_backend)
  )
  result <- run_group(group)
  saveRDS(result, args$out)
  invisible(NULL)
}

run_master <- function(args) {
  specs <- canonical_backend_specs(
    include_arrayfire = args$include_arrayfire,
    only = args$backends
  )
  specs <- filter_backend_specs(specs, args$backends)
  detected_backend_names <- vapply(specs, `[[`, character(1), "name")
  benchmark_launch_debug(
    "master detected_backends=", paste(detected_backend_names, collapse = ","),
    " ; requested=", paste(args$backends %||% character(), collapse = ","),
    " ; include_arrayfire=", args$include_arrayfire
  )
  backend_names <- if (is.null(args$backends)) {
    detected_backend_names
  } else {
    unique(args$backends[nzchar(args$backends)])
  }

  if (!is.null(args$backends)) {
    missing_backends <- setdiff(backend_names, detected_backend_names)
    if (length(missing_backends) > 0L) {
      message(
        "Requested backends unresolved during master discovery; worker validation will decide availability: ",
        paste(missing_backends, collapse = ", ")
      )
    }
  }

  worker_allowed <- vapply(
    backend_names,
    benchmark_worker_backend_allowed,
    logical(1),
    include_arrayfire = args$include_arrayfire
  )
  blocked_backends <- backend_names[!worker_allowed]
  if (length(blocked_backends) > 0L) {
    blocked_detail <- vapply(
      blocked_backends,
      benchmark_worker_backend_block_reason,
      character(1),
      include_arrayfire = args$include_arrayfire
    )
    message(
      "Skipping worker-ineligible backends in canonical harness: ",
      paste(sprintf("%s (%s)", blocked_backends, blocked_detail), collapse = ", ")
    )
    backend_names <- backend_names[worker_allowed]
  }

  if (length(backend_names) == 0L) {
    stop("No worker-eligible backends selected for canonical benchmark harness", call. = FALSE)
  }

  groups <- group_plan(backend_names, suites = args$suites)
  plan_path <- tempfile("benchmark-plan-", fileext = ".rds")
  saveRDS(groups, plan_path)
  on.exit(unlink(plan_path), add = TRUE)

  message(
    "Running canonical benchmark harness for backends: ",
    paste(backend_names, collapse = ", "),
    " | suites: ", paste(args$suites, collapse = ", ")
  )

  rows <- list()
  script_path <- normalizePath(
    getOption("amatrix.benchmark_regression.script_path", sys.frame(1)$ofile %||% "tools/benchmark-regression.R"),
    winslash = "/",
    mustWork = FALSE
  )

  for (idx in seq_along(groups)) {
    group <- groups[[idx]]
    out_path <- tempfile(sprintf("%s-", group$group_id), fileext = ".rds")
    message(sprintf("  [%02d/%02d] %s", idx, length(groups), group$group_id))
    launch <- benchmark_system2_capture(
      file.path(R.home("bin"), "Rscript"),
      benchmark_rscript_source_args(
        script_path,
        args = c("--worker", paste0("--plan=", plan_path), paste0("--group-id=", group$group_id), paste0("--out=", out_path)),
        main_call = "benchmark_regression_main(commandArgs(trailingOnly = TRUE))"
      )
    )
    log <- launch$output
    status <- launch$status

    if (status == 0L && file.exists(out_path)) {
      rows[[length(rows) + 1L]] <- readRDS(out_path)
    } else {
      failed <- expand_group_rows(group)
      failed$status <- "crash"
      failed$error_message <- sprintf(
        "worker exited with status %s%s",
        status,
        if (length(log) > 0L) paste0(": ", tail(log, 1L)) else ""
      )
      rows[[length(rows) + 1L]] <- failed
    }

    unlink(out_path)
  }

  results <- do.call(rbind, rows)
  results <- add_cpu_reference(results)
  results$speedup_vs_cpu <- ifelse(
    is.na(results$cpu_reference_ms) | is.na(results$median_ms),
    NA_real_,
    results$cpu_reference_ms / results$median_ms
  )
  results <- add_baseline_compare(results, args$baseline)
  results <- results[order(results$suite, results$op, results$size_label, results$variant, results$requested_backend), ]

  output_paths <- write_outputs(results, args$output_dir, baseline_path = args$baseline)

  history_path <- getOption(
    "amatrix.benchmark_history_path",
    file.path(dirname(args$baseline), "benchmark-history.csv")
  )
  tryCatch(
    append_benchmark_history(
      results[results$status == "ok", , drop = FALSE],
      history_path
    ),
    error = function(e) {
      message("benchmark-history append failed: ", conditionMessage(e))
    }
  )

  if (args$update || !file.exists(args$baseline)) {
    baseline_rows <- results[results$status == "ok", c(key_columns, "status", "median_ms")]
    write.csv(baseline_rows, args$baseline, row.names = FALSE)
    message("\nBaseline written to ", args$baseline)
  }

  print_run_summary(results, args$baseline, output_paths, warm_only = args$warm_only)
  invisible(results)
}

benchmark_regression_dispatch <- function(command_args = commandArgs(trailingOnly = TRUE), allow_relaunch = TRUE) {
  args <- parse_args(command_args)
  if (isTRUE(args$help)) {
    cat(benchmark_regression_usage(), sep = "\n")
    return(invisible(TRUE))
  }
  if (isTRUE(allow_relaunch)) {
    relaunch_safe_master_if_needed(args)
  }
  initialize_regression_benchmark_context()

  if (isTRUE(args$worker)) {
    run_worker(args)
  } else {
    run_master(args)
  }

  invisible(TRUE)
}

benchmark_regression_main <- function(command_args = commandArgs(trailingOnly = TRUE)) {
  benchmark_regression_dispatch(command_args = command_args, allow_relaunch = TRUE)
}

benchmark_regression_direct_file_info <- function() {
  raw_args <- commandArgs(trailingOnly = FALSE)
  direct_file_paths <- sub("^--file=", "", grep("^--file=", raw_args, value = TRUE))
  if (length(direct_file_paths) == 0L) {
    script_like_paths <- raw_args[basename(raw_args) == "benchmark-regression.R"]
    if (length(script_like_paths) > 0L) {
      direct_file_paths <- script_like_paths
    } else {
      existing_paths <- raw_args[file.exists(raw_args)]
      direct_file_paths <- existing_paths[basename(existing_paths) == "benchmark-regression.R"]
    }
  }

  if (length(direct_file_paths) == 0L) {
    return(NULL)
  }
  direct_file_paths <- direct_file_paths[basename(direct_file_paths) == "benchmark-regression.R"]
  if (length(direct_file_paths) == 0L) {
    return(NULL)
  }

  script_path <- normalizePath(direct_file_paths[[1L]], winslash = "/", mustWork = TRUE)
  list(
    script_path = script_path,
    repo_root = normalizePath(dirname(dirname(script_path)), winslash = "/", mustWork = TRUE)
  )
}

benchmark_regression_invoked_directly <- function() {
  !is.null(benchmark_regression_direct_file_info())
}

benchmark_regression_autorun_enabled <- function() {
  isTRUE(getOption("amatrix.benchmark_regression.autorun", TRUE))
}

benchmark_regression_direct_relaunch <- function(command_args = commandArgs(trailingOnly = TRUE)) {
  info <- benchmark_regression_direct_file_info()
  if (is.null(info)) {
    return(invisible(FALSE))
  }

  cli_info <- benchmark_regression_cli_script_path(info$script_path)
  relaunch_args <- c(cli_info$cli_path, command_args)
  warned_status <- NULL
  relaunch_output <- withCallingHandlers(
    system2(
      file.path(R.home("bin"), "Rscript"),
      vapply(relaunch_args, shQuote, character(1), USE.NAMES = FALSE),
      stdout = TRUE,
      stderr = TRUE
    ),
    warning = function(w) {
      warned_status <<- attr(w, "status") %||% warned_status
      invokeRestart("muffleWarning")
    }
  )
  relaunch_status <- attr(relaunch_output, "status") %||% warned_status %||% 0L

  if (length(relaunch_output) > 0L) {
    cat(paste(relaunch_output, collapse = "\n"), sep = "\n")
    if (!grepl("\n$", paste(relaunch_output, collapse = "\n"))) {
      cat("\n")
    }
  }

  quit(save = "no", status = relaunch_status)
}

if (benchmark_regression_autorun_enabled() && benchmark_regression_invoked_directly()) {
  direct_args <- commandArgs(trailingOnly = TRUE)
  parsed_args <- parse_args(direct_args)

  if (!isTRUE(parsed_args$safe_main) && !isTRUE(parsed_args$worker)) {
    benchmark_regression_direct_relaunch(direct_args)
  } else {
    benchmark_regression_main(direct_args)
  }
}
