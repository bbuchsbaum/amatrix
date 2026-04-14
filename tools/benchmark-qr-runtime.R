#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(amatrix)
  }
})

if (!requireNamespace("amatrix.mlx", quietly = TRUE)) {
  message("skipped: amatrix.mlx not installed")
  quit(save = "no", status = 0L)
}

benchmark_elapsed <- function(fn, reps = 5L, iterations = 20L, warmup = NULL) {
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

with_qr_mode <- function(mode, expr) {
  old <- getOption("amatrix.mlx.qr_helper_mode")
  options(amatrix.mlx.qr_helper_mode = mode)
  on.exit(options(amatrix.mlx.qr_helper_mode = old), add = TRUE)
  force(expr)
}

benchmark_case <- function(n, p, rhs_cols = 8L) {
  set.seed(20260404 + n + p)
  iterations <- if (n <= 256L) 50L else if (n <= 512L) 20L else 5L
  if (rhs_cols >= 64L) {
    iterations <- max(2L, iterations %/% 2L)
  }
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y <- matrix(rnorm(n * rhs_cols), nrow = n, ncol = rhs_cols)
  x_mlx <- adgeMatrix(x, preferred_backend = "mlx", precision = "fast")

  base_fac <- base::qr(x)
  native_fac <- with_qr_mode("native", qr(x_mlx))
  compact_fac <- with_qr_mode("compact", qr(x_mlx))

  list(
    qr_factor = data.frame(
      case = sprintf("%dx%d", n, p),
      rhs_cols = rhs_cols,
      workload = "qr.factor_cold",
      runtime = c("base_r", "mlx_resident_qr"),
      elapsed = c(
        benchmark_elapsed(function() base::qr(x), iterations = iterations),
        benchmark_elapsed(function() with_qr_mode("native", qr(x_mlx)), iterations = iterations)
      ),
      stringsAsFactors = FALSE
    ),
    qr_q = data.frame(
      case = sprintf("%dx%d", n, p),
      rhs_cols = rhs_cols,
      workload = "qr.Q_materialize",
      runtime = c("base_r", "mlx_native_resident"),
      elapsed = c(
        benchmark_elapsed(function() base::qr.Q(base_fac), iterations = iterations),
        benchmark_elapsed(function() with_qr_mode("native", qr.Q(native_fac)), iterations = iterations)
      ),
      stringsAsFactors = FALSE
    ),
    qr_coef = data.frame(
      case = sprintf("%dx%d", n, p),
      rhs_cols = rhs_cols,
      workload = "qr.coef",
      runtime = c("base_r", "mlx_native_resident", "mlx_compact"),
      elapsed = c(
        benchmark_elapsed(function() base::qr.coef(base_fac, y), iterations = iterations),
        benchmark_elapsed(function() with_qr_mode("native", qr.coef(native_fac, y)), iterations = iterations),
        benchmark_elapsed(function() with_qr_mode("compact", qr.coef(compact_fac, y)), iterations = iterations)
      ),
      stringsAsFactors = FALSE
    ),
    qr_qty = data.frame(
      case = sprintf("%dx%d", n, p),
      rhs_cols = rhs_cols,
      workload = "qr.qty",
      runtime = c("base_r", "mlx_native_resident", "mlx_compact"),
      elapsed = c(
        benchmark_elapsed(function() base::qr.qty(base_fac, y), iterations = iterations),
        benchmark_elapsed(function() with_qr_mode("native", qr.qty(native_fac, y)), iterations = iterations),
        benchmark_elapsed(function() with_qr_mode("compact", qr.qty(compact_fac, y)), iterations = iterations)
      ),
      stringsAsFactors = FALSE
    ),
    lm_fit = data.frame(
      case = sprintf("%dx%d", n, p),
      rhs_cols = rhs_cols,
      workload = "am_lm_fit_qr_hot",
      runtime = c("base_qr_cached", "mlx_native_resident", "mlx_compact"),
      elapsed = c(
        benchmark_elapsed(
          function() base::qr.solve(base_fac, y),
          iterations = iterations
        ),
        benchmark_elapsed(
          function() with_qr_mode("native", am_lm_fit(x_mlx, y, method = "qr", include_fitted = FALSE, include_residuals = FALSE, cache = TRUE)),
          warmup = function() with_qr_mode("native", am_lm_fit(x_mlx, y, method = "qr", include_fitted = FALSE, include_residuals = FALSE, cache = TRUE)),
          iterations = iterations
        ),
        benchmark_elapsed(
          function() with_qr_mode("compact", am_lm_fit(x_mlx, y, method = "qr", include_fitted = FALSE, include_residuals = FALSE, cache = TRUE)),
          warmup = function() with_qr_mode("compact", am_lm_fit(x_mlx, y, method = "qr", include_fitted = FALSE, include_residuals = FALSE, cache = TRUE)),
          iterations = iterations
        )
      ),
      stringsAsFactors = FALSE
    ),
    many_lm = data.frame(
      case = sprintf("%dx%d", n, p),
      rhs_cols = rhs_cols,
      workload = "many_lm_qr_hot",
      runtime = c("base_qr_cached", "mlx_native_resident", "mlx_compact"),
      elapsed = c(
        benchmark_elapsed(
          function() base::qr.solve(base_fac, y),
          iterations = iterations
        ),
        benchmark_elapsed(
          function() with_qr_mode("native", many_lm(x_mlx, y, include_residuals = FALSE, cache = TRUE, method = "qr")),
          warmup = function() with_qr_mode("native", many_lm(x_mlx, y, include_residuals = FALSE, cache = TRUE, method = "qr")),
          iterations = iterations
        ),
        benchmark_elapsed(
          function() with_qr_mode("compact", many_lm(x_mlx, y, include_residuals = FALSE, cache = TRUE, method = "qr")),
          warmup = function() with_qr_mode("compact", many_lm(x_mlx, y, include_residuals = FALSE, cache = TRUE, method = "qr")),
          iterations = iterations
        )
      ),
      stringsAsFactors = FALSE
    )
  )
}

# Activate the Metal GPU probe (safe in -e / interactive launch mode).
# In direct `Rscript file.R` mode, MLX is skipped rather than crashing.
Sys.setenv(AMATRIX_MLX_PROBE_GPU = "1")
if (!amatrix_mlx_is_available()) {
  message("MLX not available. To benchmark MLX, run via:\n",
          "  Rscript -e 'Sys.setenv(AMATRIX_MLX_PROBE_GPU=\"1\"); ",
          "setwd(\"", normalizePath("."), "\"); ",
          "source(\"tools/benchmark-qr-runtime.R\")'")
}

cases <- list(c(512L, 64L, 8L), c(1024L, 128L, 8L), c(1024L, 128L, 32L), c(1024L, 128L, 128L))
rows <- list()
for (dims in cases) {
  out <- benchmark_case(dims[[1]], dims[[2]], rhs_cols = dims[[3]])
  rows <- c(rows, out)
}

results <- do.call(rbind, rows)
row.names(results) <- NULL
results$elapsed <- sprintf("%.6f", results$elapsed)

cat("Notes:\n")
cat("- mlx_resident_qr keeps the MLX Q factor resident and returns only a resident key plus host R.\n")
cat("- mlx_native_resident uses backend-native helper ops against that resident Q (`Q^T Y`, `QY`, triangular solve) where possible.\n")
cat("- mlx_compact uses the compact MLX QR path; on tall-skinny matrices this currently means a TSQR-style blocked compact factor prototype, otherwise it falls back to the bridge-compact factor path.\n")
cat("- base_qr_cached measures the cached base QR helper path rather than a full cold re-factorization.\n\n")
print(results)
