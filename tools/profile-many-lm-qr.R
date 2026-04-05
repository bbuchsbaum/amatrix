#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(amatrix)
  library(amatrix.mlx)
})

profile_case <- function(n = 4096L, p = 128L, rhs_cols = 128L, block_rows = 512L) {
  set.seed(20260404 + n + p + rhs_cols)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Y <- matrix(rnorm(n * rhs_cols), nrow = n, ncol = rhs_cols)
  x <- adgeMatrix(X, preferred_backend = "mlx", precision = "fast")

  options(
    amatrix.mlx.available = TRUE,
    amatrix.mlx.qr_helper_mode = "compact",
    amatrix.mlx.qr_compact_method = "tsqr",
    amatrix.mlx.qr_tsqr_block_rows = block_rows
  )

  profile <- amatrix:::.amatrix_profile_many_lm_qr(x, Y, cache = TRUE)

  cat(sprintf("case=%dx%d rhs=%d block_rows=%d\n", n, p, rhs_cols, block_rows))
  cat(sprintf("cache_key=%s\n", profile$cache_key))
  cat(sprintf("cache_reused=%s\n", profile$cache_reused))
  cat(sprintf("qr_representation=%s\n", profile$qr_representation))
  cat(sprintf("qr_helper_path=%s\n", profile$qr_helper_path))
  cat(sprintf("qr_compact_factor_source=%s\n", profile$qr_compact_factor_source))
  cat(sprintf("timings(cache)=%.6f\n", profile$timings[["cache"]]))
  cat(sprintf("timings(solve)=%.6f\n", profile$timings[["solve"]]))
  cat(sprintf("timings(assemble)=%.6f\n", profile$timings[["assemble"]]))
}

profile_case()
