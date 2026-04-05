#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(amatrix)
  library(amatrix.mlx)
})

benchmark_elapsed <- function(fn, reps = 5L, iterations = 10L, warmup = NULL) {
  if (is.function(warmup)) {
    warmup()
  }
  timings <- numeric(reps)
  for (idx in seq_len(reps)) {
    gc()
    timings[[idx]] <- system.time(
      for (iter in seq_len(iterations)) {
        fn()
      }
    )[["elapsed"]] / iterations
  }
  median(timings)
}

with_qr_options <- function(helper_mode, compact_method = NULL, block_rows = NULL, expr) {
  old <- options()
  on.exit(options(old), add = TRUE)
  opts <- list(amatrix.mlx.qr_helper_mode = helper_mode)
  if (!is.null(compact_method)) {
    opts$amatrix.mlx.qr_compact_method <- compact_method
  }
  if (!is.null(block_rows)) {
    opts$amatrix.mlx.qr_tsqr_block_rows <- as.integer(block_rows)
  }
  options(opts)
  force(expr)
}

benchmark_case <- function(n, p, rhs_cols, block_rows = NULL) {
  set.seed(20260404 + n + p + rhs_cols)
  iterations <- if (rhs_cols >= 128L) 3L else if (rhs_cols >= 32L) 5L else 8L
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y <- matrix(rnorm(n * rhs_cols), nrow = n, ncol = rhs_cols)
  x_mlx <- adgeMatrix(x, preferred_backend = "mlx", precision = "fast")

  base_fac <- base::qr(x)
  native_fac <- with_qr_options("native", expr = qr(x_mlx))
  compact_fac <- with_qr_options("compact", compact_method = "tsqr", block_rows = block_rows, expr = qr(x_mlx))

  rbind(
    data.frame(
      case = sprintf("%dx%d", n, p),
      rhs_cols = rhs_cols,
      block_rows = if (is.null(block_rows)) NA_integer_ else as.integer(block_rows),
      workload = "qr.coef",
      runtime = c("base_qr_cached", "mlx_native_resident", "mlx_compact_tsqr"),
      elapsed = c(
        benchmark_elapsed(function() base::qr.coef(base_fac, y), iterations = iterations),
        benchmark_elapsed(function() with_qr_options("native", expr = qr.coef(native_fac, y)), iterations = iterations),
        benchmark_elapsed(
          function() with_qr_options("compact", compact_method = "tsqr", block_rows = block_rows, expr = qr.coef(compact_fac, y)),
          iterations = iterations
        )
      ),
      stringsAsFactors = FALSE
    ),
    data.frame(
      case = sprintf("%dx%d", n, p),
      rhs_cols = rhs_cols,
      block_rows = if (is.null(block_rows)) NA_integer_ else as.integer(block_rows),
      workload = "am_many_lm_qr_hot",
      runtime = c("base_qr_cached", "mlx_native_resident", "mlx_compact_tsqr"),
      elapsed = c(
        benchmark_elapsed(function() base::qr.solve(base_fac, y), iterations = iterations),
        benchmark_elapsed(
          function() with_qr_options("native", expr = am_many_lm(x_mlx, y, include_residuals = FALSE, cache = TRUE, method = "qr")),
          warmup = function() with_qr_options("native", expr = am_many_lm(x_mlx, y, include_residuals = FALSE, cache = TRUE, method = "qr")),
          iterations = iterations
        ),
        benchmark_elapsed(
          function() with_qr_options("compact", compact_method = "tsqr", block_rows = block_rows, expr = am_many_lm(x_mlx, y, include_residuals = FALSE, cache = TRUE, method = "qr")),
          warmup = function() with_qr_options("compact", compact_method = "tsqr", block_rows = block_rows, expr = am_many_lm(x_mlx, y, include_residuals = FALSE, cache = TRUE, method = "qr")),
          iterations = iterations
        )
      ),
      stringsAsFactors = FALSE
    )
  )
}

if (!amatrix_mlx_is_available()) {
  options(amatrix.mlx.available = TRUE)
}

cases <- list(
  c(1024L, 128L, 8L, 256L),
  c(1024L, 128L, 32L, 256L),
  c(1024L, 128L, 128L, 256L),
  c(4096L, 128L, 32L, 512L),
  c(4096L, 128L, 128L, 512L)
)

rows <- do.call(
  rbind,
  lapply(cases, function(case) benchmark_case(case[[1]], case[[2]], case[[3]], block_rows = case[[4]]))
)
row.names(rows) <- NULL
rows$elapsed <- sprintf("%.6f", rows$elapsed)

cat("Notes:\n")
cat("- This harness isolates tall-skinny shared-X many-RHS QR.\n")
cat("- mlx_compact_tsqr forces the compact TSQR path.\n")
cat("- base_qr_cached measures the cached base QR helper path.\n")
cat("- am_many_lm_qr_hot is the flagship workload and excludes fitted/residual reconstruction.\n\n")
print(rows)
