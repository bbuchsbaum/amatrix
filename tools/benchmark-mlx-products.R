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

options(amatrix.mlx.available = TRUE)

sizes <- c(512L, 1024L, 1536L)
iterations <- 5L

run_case <- function(n) {
  set.seed(n)
  x_host <- matrix(rnorm(n * n), nrow = n)
  x_cpu <- adgeMatrix(x_host, preferred_backend = "cpu")
  x_mlx <- adgeMatrix(x_host, preferred_backend = "mlx")

  x_mlx <- prime_backend(x_mlx, "mlx")

  bench::mark(
    cpu_crossprod = invisible(crossprod(x_cpu)),
    mlx_crossprod = invisible(crossprod(x_mlx)),
    cpu_tcrossprod = invisible(tcrossprod(x_cpu)),
    mlx_tcrossprod = invisible(tcrossprod(x_mlx)),
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
