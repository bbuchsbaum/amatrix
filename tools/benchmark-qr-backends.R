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

if (!requireNamespace("amatrix.arrayfire", quietly = TRUE)) {
  message("skipped: amatrix.arrayfire not installed")
  quit(save = "no", status = 0L)
}

benchmark_elapsed <- function(fn, reps = 5L, iterations = 5L) {
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

explicit_qr_coef <- function(q, r, y, p) {
  qty <- crossprod(q, y)
  backsolve(r[seq_len(p), seq_len(p), drop = FALSE], qty[seq_len(p), , drop = FALSE])
}

relative_error <- function(x, y) {
  denom <- max(1, max(abs(y)))
  max(abs(x - y)) / denom
}

cases <- list(
  square_32 = c(32L, 32L),
  square_64 = c(64L, 64L),
  square_80 = c(80L, 80L)
)

set.seed(20260404)
rows <- list()

for (case_name in names(cases)) {
  dims <- cases[[case_name]]
  n <- dims[[1]]
  p <- dims[[2]]
  iterations <- if (n <= 32L) 1000L else if (n <= 64L) 500L else 250L
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y <- matrix(rnorm(n * 4L), nrow = n, ncol = 4L)
  base_fit <- base::qr(x)
  base_q <- qr.Q(base_fit)
  base_r <- qr.R(base_fit)
  base_coef <- base::qr.coef(base_fit, y)

  base_elapsed <- benchmark_elapsed(function() base::qr(x), iterations = iterations)
  mlx_elapsed <- benchmark_elapsed(function() amatrix.mlx:::amatrix_mlx_qr(x), iterations = iterations)

  mlx_fit <- amatrix.mlx:::amatrix_mlx_qr(x)
  mlx_coef <- explicit_qr_coef(mlx_fit$q, mlx_fit$r, y, p)

  rows[[length(rows) + 1L]] <- data.frame(
    case = case_name,
    backend = "base_r",
    elapsed = base_elapsed,
    rel_reconstruction_error = relative_error(base_q %*% base_r, x),
    rel_coef_error = relative_error(base_coef, base_coef),
    note = "base::qr",
    stringsAsFactors = FALSE
  )

  rows[[length(rows) + 1L]] <- data.frame(
    case = case_name,
    backend = "mlx",
    elapsed = mlx_elapsed,
    rel_reconstruction_error = relative_error(mlx_fit$q %*% mlx_fit$r, x),
    rel_coef_error = relative_error(mlx_coef, base_coef),
    note = "native MLX QR; current MLX runtime only supports QR on CPU stream",
    stringsAsFactors = FALSE
  )

  if (n >= 96L) {
    message(sprintf("skipped: known ArrayFire 96x96 segfault, see .omc/research/bench-audit/bugs.md (case=%s, n=%d)", case_name, n))
    rows[[length(rows) + 1L]] <- data.frame(
      case = case_name,
      backend = "arrayfire",
      elapsed = NA_real_,
      rel_reconstruction_error = NA_real_,
      rel_coef_error = NA_real_,
      note = "skipped: known ArrayFire 96x96 segfault, see .omc/research/bench-audit/bugs.md",
      stringsAsFactors = FALSE
    )
  } else {
    af_elapsed <- benchmark_elapsed(function() amatrix.arrayfire:::amatrix_arrayfire_qr(x), iterations = iterations)
    af_fit <- amatrix.arrayfire:::amatrix_arrayfire_qr(x)
    af_coef <- explicit_qr_coef(af_fit$q, af_fit$r, y, p)
    rows[[length(rows) + 1L]] <- data.frame(
      case = case_name,
      backend = "arrayfire",
      elapsed = af_elapsed,
      rel_reconstruction_error = relative_error(af_fit$q %*% af_fit$r, x),
      rel_coef_error = relative_error(af_coef, base_coef),
      note = "native ArrayFire QR",
      stringsAsFactors = FALSE
    )
  }
}

results <- do.call(rbind, rows)
row.names(results) <- NULL
results$elapsed <- sprintf("%.6f", results$elapsed)
results$rel_reconstruction_error <- sprintf("%.6e", results$rel_reconstruction_error)
results$rel_coef_error <- sprintf("%.6e", results$rel_coef_error)
cat("Notes:\n")
cat("- MLX QR is native, but the installed MLX runtime reports QR as unsupported on GPU; the MLX bridge uses a CPU stream.\n")
cat("- ArrayFire QR is stable up to 80x80 on this machine and segfaults at 96x96 and above.\n\n")
print(results)
