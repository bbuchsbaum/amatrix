#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(amatrix)
  }
  library(amatrix.mlx)
})

bench_expr <- function(fn, reps = 3L) {
  stopifnot(is.function(fn))
  system.time({
    for (i in seq_len(reps)) {
      invisible(fn())
    }
  })[["elapsed"]] / reps
}

chain_workload <- function(x, y, z) {
  ((x %*% y) * 2) + z
}

cross_workload <- function(x) {
  (crossprod(x) * 2) + diag(ncol(x))
}

bench_chain <- function(n, reps = 3L) {
  set.seed(1)
  x <- matrix(rnorm(n * n), n, n)
  y <- matrix(rnorm(n * n), n, n)
  z <- matrix(rnorm(n * n), n, n)

  cpu <- adgeMatrix(x, preferred_backend = "cpu")
  mlx <- adgeMatrix(x, preferred_backend = "mlx")
  z_mlx <- adgeMatrix(z, preferred_backend = "mlx")

  data.frame(
      workload = c("matmul_chain", "matmul_chain", "cross_chain", "cross_chain"),
      backend = c("cpu", "mlx", "cpu", "mlx"),
      n = n,
      elapsed = c(
      bench_expr(function() chain_workload(cpu, y, z), reps = reps),
      bench_expr(function() chain_workload(mlx, y, z_mlx), reps = reps),
      bench_expr(function() cross_workload(cpu), reps = reps),
      bench_expr(function() cross_workload(mlx), reps = reps)
    )
  )
}

results <- do.call(rbind, lapply(c(256L, 512L, 1024L), bench_chain))
print(results)
