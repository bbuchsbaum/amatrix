#!/usr/bin/env Rscript

source(file.path("tools", "benchmark-helpers.R"), local = TRUE)
load_benchmark_amatrix()

`%||%` <- function(x, y) if (is.null(x)) y else x

available_qr_backends <- function() {
  benchmark_backend_names(
    include_cpu = TRUE,
    include_mlx = FALSE,
    include_metal = FALSE,
    include_opencl = TRUE
  )
}

make_qr_operand <- function(x, backend) {
  precision <- if (identical(backend, "cpu")) "strict" else "fast"
  adgeMatrix(x, preferred_backend = backend, precision = precision)
}

benchmark_elapsed <- function(fn, reps = 5L) {
  timings <- numeric(reps)
  for (idx in seq_len(reps)) {
    gc()
    timings[[idx]] <- system.time(fn())[["elapsed"]]
  }
  median(timings)
}

dispatch_backend <- function(x, op, y = NULL) {
  amatrix_backend_plan(x, op, y = y)$chosen
}

experimental_qr_solve_requested <- function() {
  identical(Sys.getenv("AMATRIX_OPENCL_EXPERIMENTAL_QR_SOLVE", unset = ""), "1")
}

relative_error <- function(actual, expected) {
  denom <- max(1, max(abs(expected)))
  max(abs(actual - expected)) / denom
}

qr_case <- function(kind, n, p, rhs, seed) {
  set.seed(seed)
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y <- matrix(rnorm(n * rhs), nrow = n, ncol = rhs)

  if (identical(kind, "square")) {
    x <- x + diag(p) * 0.5
  }

  list(X = x, Y = y, kind = kind, n = n, p = p, rhs = rhs)
}

qr_cases <- function() {
  list(
    tall_2000x32 = qr_case("tall", n = 2000L, p = 32L, rhs = 4L, seed = 20260409L),
    tall_8000x64 = qr_case("tall", n = 8000L, p = 64L, rhs = 8L, seed = 20260410L),
    square_768 = qr_case("square", n = 768L, p = 768L, rhs = 8L, seed = 20260411L)
  )
}

bench_qr_case <- function(case_name, case, backend, reps = 5L) {
  x_arg <- make_qr_operand(case$X, backend)
  y_arg <- make_qr_operand(case$Y, backend)
  qr_fit <- qr(x_arg)
  qr_meta <- tryCatch(qr_info(qr_fit), error = function(e) NULL)
  qr_ref <- qr(case$X)
  coef_ref <- qr.coef(qr_ref, case$Y)
  fitted_ref <- qr.fitted(qr_ref, case$Y)
  resid_ref <- qr.resid(qr_ref, case$Y)
  qr_dispatch <- dispatch_backend(x_arg, "qr")

  rows <- list(
    data.frame(
      case = case_name,
      kind = case$kind,
      backend = backend,
      op = "factor",
      elapsed = benchmark_elapsed(function() qr(make_qr_operand(case$X, backend)), reps = reps),
      rel_error = 0,
      dispatch_backend = qr_dispatch,
      helper_path = if (is.null(qr_meta)) NA_character_ else qr_meta$helper_path,
      representation = if (is.null(qr_meta)) NA_character_ else qr_meta$representation,
      factor_source = if (is.null(qr_meta)) NA_character_ else qr_meta$compact_factor_source %||% qr_fit$state$factor_source %||% NA_character_,
      stringsAsFactors = FALSE
    ),
    data.frame(
      case = case_name,
      kind = case$kind,
      backend = backend,
      op = "coef_cached",
      elapsed = benchmark_elapsed(function() qr.coef(qr_fit, case$Y), reps = reps),
      rel_error = relative_error(as.matrix(qr.coef(qr_fit, case$Y)), coef_ref),
      dispatch_backend = qr_dispatch,
      helper_path = if (is.null(qr_meta)) NA_character_ else qr_meta$helper_path,
      representation = if (is.null(qr_meta)) NA_character_ else qr_meta$representation,
      factor_source = if (is.null(qr_meta)) NA_character_ else qr_meta$compact_factor_source %||% qr_fit$state$factor_source %||% NA_character_,
      stringsAsFactors = FALSE
    ),
    data.frame(
      case = case_name,
      kind = case$kind,
      backend = backend,
      op = "fitted_cached",
      elapsed = benchmark_elapsed(function() qr.fitted(qr_fit, case$Y), reps = reps),
      rel_error = relative_error(as.matrix(qr.fitted(qr_fit, case$Y)), fitted_ref),
      dispatch_backend = qr_dispatch,
      helper_path = if (is.null(qr_meta)) NA_character_ else qr_meta$helper_path,
      representation = if (is.null(qr_meta)) NA_character_ else qr_meta$representation,
      factor_source = if (is.null(qr_meta)) NA_character_ else qr_meta$compact_factor_source %||% qr_fit$state$factor_source %||% NA_character_,
      stringsAsFactors = FALSE
    ),
    data.frame(
      case = case_name,
      kind = case$kind,
      backend = backend,
      op = "resid_cached",
      elapsed = benchmark_elapsed(function() qr.resid(qr_fit, case$Y), reps = reps),
      rel_error = relative_error(as.matrix(qr.resid(qr_fit, case$Y)), resid_ref),
      dispatch_backend = qr_dispatch,
      helper_path = if (is.null(qr_meta)) NA_character_ else qr_meta$helper_path,
      representation = if (is.null(qr_meta)) NA_character_ else qr_meta$representation,
      factor_source = if (is.null(qr_meta)) NA_character_ else qr_meta$compact_factor_source %||% qr_fit$state$factor_source %||% NA_character_,
      stringsAsFactors = FALSE
    )
  )

  if (identical(case$kind, "square")) {
    solve_ref <- solve(case$X, case$Y)
    solve_dispatch <- dispatch_backend(x_arg, "solve", y = y_arg)
    rows[[length(rows) + 1L]] <- data.frame(
      case = case_name,
      kind = case$kind,
      backend = backend,
      op = "solve_rhs",
      elapsed = benchmark_elapsed(function() solve(x_arg, case$Y), reps = reps),
      rel_error = relative_error(as.matrix(solve(x_arg, case$Y)), solve_ref),
      dispatch_backend = solve_dispatch,
      helper_path = if (is.null(qr_meta)) NA_character_ else qr_meta$helper_path,
      representation = if (is.null(qr_meta)) NA_character_ else qr_meta$representation,
      factor_source = if (is.null(qr_meta)) NA_character_ else qr_meta$compact_factor_source %||% qr_fit$state$factor_source %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

benchmark_opencl_qr <- function(backends = available_qr_backends(), reps = 5L) {
  rows <- lapply(names(qr_cases()), function(case_name) {
    case <- qr_cases()[[case_name]]
    do.call(rbind, lapply(backends, function(backend) bench_qr_case(case_name, case, backend, reps = reps)))
  })
  do.call(rbind, rows)
}

if (sys.nframe() == 0L) {
  old <- options(amatrix.opencl.factor_gpu = TRUE)
  if (experimental_qr_solve_requested()) {
    old <- c(
      old,
      options(
        amatrix.opencl.experimental_qr_solve = TRUE,
        amatrix.opencl.solve_qr_min_dim = 1L,
        amatrix.opencl.solve_qr_min_rhs = 1L
      )
    )
  }
  on.exit(options(old), add = TRUE)
  print(benchmark_opencl_qr())
}
