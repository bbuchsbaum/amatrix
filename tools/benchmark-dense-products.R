#!/usr/bin/env Rscript

if (file.exists(file.path("tools", "benchmark-helpers.R"))) {
  source(file.path("tools", "benchmark-helpers.R"), local = FALSE)
}

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(amatrix)
  }
  library(bench)
})

if (!requireNamespace("amatrix.mlx", quietly = TRUE)) {
  message("skipped: amatrix.mlx not installed")
  quit(save = "no", status = 0L)
}

sizes <- c(256L, 512L, 1024L)
iterations <- 5L

options(amatrix.mlx.available = TRUE)

run_case <- function(n) {
  set.seed(n)
  x_host <- matrix(rnorm(n * n), nrow = n)
  y_host <- matrix(rnorm(n * n), nrow = n)

  x_cpu <- adgeMatrix(x_host, preferred_backend = "cpu")
  y_cpu <- adgeMatrix(y_host, preferred_backend = "cpu")
  x_mlx <- adgeMatrix(x_host, preferred_backend = "mlx")
  y_mlx <- adgeMatrix(y_host, preferred_backend = "mlx")

  x_mlx <- prime_backend(x_mlx, "mlx")
  y_mlx <- prime_backend(y_mlx, "mlx")

  bench::mark(
    cpu_matmul = { invisible(x_cpu %*% y_cpu) },
    mlx_matmul = { invisible(x_mlx %*% y_mlx) },
    cpu_crossprod = { invisible(crossprod(x_cpu)) },
    mlx_crossprod = { invisible(crossprod(x_mlx)) },
    cpu_tcrossprod = { invisible(tcrossprod(x_cpu)) },
    mlx_tcrossprod = { invisible(tcrossprod(x_mlx)) },
    iterations = iterations,
    check = FALSE,
    memory = FALSE,
    time_unit = "ms"
  )
}

results <- lapply(sizes, function(n) {
  out <- run_case(n)
  out$size <- n
  out
})

summary <- do.call(rbind, results)[, c("size", "expression", "median", "itr/sec")]
print(summary)
