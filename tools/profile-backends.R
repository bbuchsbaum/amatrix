#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(amatrix)
})

load_optional_backend <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressWarnings(suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    ))
    TRUE
  } else {
    FALSE
  }
}

have_mlx <- load_optional_backend("amatrix.mlx")
have_arrayfire <- load_optional_backend("amatrix.arrayfire")

bench_expr <- function(fn, reps = 3L) {
  stopifnot(is.function(fn))
  system.time({
    for (i in seq_len(reps)) {
      invisible(fn())
    }
  })[["elapsed"]] / reps
}

row_sums_method <- methods::selectMethod("rowSums", "adgeMatrix")
col_sums_method <- methods::selectMethod("colSums", "adgeMatrix")

profile_product_case <- function(preferred_backend, op, x, y, reps = 3L) {
  chosen <- {
    obj <- adgeMatrix(x, preferred_backend = preferred_backend)
    if (identical(op, "matmul")) {
      amatrix_backend_plan(obj, op, y = y)$chosen
    } else {
      amatrix_backend_plan(obj, op)$chosen
    }
  }

  elapsed <- bench_expr(function() {
    obj <- adgeMatrix(x, preferred_backend = preferred_backend)
    switch(
      op,
      matmul = obj %*% y,
      crossprod = crossprod(obj),
      tcrossprod = tcrossprod(obj),
      stop("unsupported product op")
    )
  }, reps = reps)

  c(chosen = chosen, elapsed = elapsed)
}

profile_reduction_case <- function(preferred_backend, op, x, reps = 3L) {
  chosen <- {
    obj <- adgeMatrix(x, preferred_backend = preferred_backend)
    amatrix_backend_plan(obj, if (identical(op, "ewise_mul")) "ewise" else op)$chosen
  }

  elapsed <- bench_expr(function() {
    obj <- adgeMatrix(x, preferred_backend = preferred_backend)
    switch(
      op,
      ewise_mul = obj * 2,
      rowSums = row_sums_method(obj),
      colSums = col_sums_method(obj),
      stop("unsupported reduction op")
    )
  }, reps = reps)

  c(chosen = chosen, elapsed = elapsed)
}

profile_products <- function(n, reps = 3L) {
  set.seed(1)
  x <- matrix(rnorm(n * n), n, n)
  y <- matrix(rnorm(n * n), n, n)
  cpu_matmul <- profile_product_case("cpu", "matmul", x, y, reps = reps)
  cpu_crossprod <- profile_product_case("cpu", "crossprod", x, y, reps = reps)
  cpu_tcrossprod <- profile_product_case("cpu", "tcrossprod", x, y, reps = reps)
  mlx_matmul <- profile_product_case("mlx", "matmul", x, y, reps = reps)
  mlx_crossprod <- profile_product_case("mlx", "crossprod", x, y, reps = reps)
  mlx_tcrossprod <- profile_product_case("mlx", "tcrossprod", x, y, reps = reps)
  af_matmul <- profile_product_case("arrayfire", "matmul", x, y, reps = reps)
  af_crossprod <- profile_product_case("arrayfire", "crossprod", x, y, reps = reps)
  af_tcrossprod <- profile_product_case("arrayfire", "tcrossprod", x, y, reps = reps)

  data.frame(
    family = "products",
    n = n,
    backend = rep(c("cpu", "mlx", "arrayfire"), each = 3),
    op = rep(c("matmul", "crossprod", "tcrossprod"), 3),
    chosen = c(
      cpu_matmul[["chosen"]],
      cpu_crossprod[["chosen"]],
      cpu_tcrossprod[["chosen"]],
      mlx_matmul[["chosen"]],
      mlx_crossprod[["chosen"]],
      mlx_tcrossprod[["chosen"]],
      af_matmul[["chosen"]],
      af_crossprod[["chosen"]],
      af_tcrossprod[["chosen"]]
    ),
    elapsed = c(
      as.numeric(cpu_matmul[["elapsed"]]),
      as.numeric(cpu_crossprod[["elapsed"]]),
      as.numeric(cpu_tcrossprod[["elapsed"]]),
      as.numeric(mlx_matmul[["elapsed"]]),
      as.numeric(mlx_crossprod[["elapsed"]]),
      as.numeric(mlx_tcrossprod[["elapsed"]]),
      as.numeric(af_matmul[["elapsed"]]),
      as.numeric(af_crossprod[["elapsed"]]),
      as.numeric(af_tcrossprod[["elapsed"]])
    )
  )
}

profile_reductions <- function(n, reps = 3L) {
  set.seed(1)
  x <- matrix(rnorm(n * n), n, n)
  cpu_ewise <- profile_reduction_case("cpu", "ewise_mul", x, reps = reps)
  cpu_rowSums <- profile_reduction_case("cpu", "rowSums", x, reps = reps)
  cpu_colSums <- profile_reduction_case("cpu", "colSums", x, reps = reps)
  af_ewise <- profile_reduction_case("arrayfire", "ewise_mul", x, reps = reps)
  af_rowSums <- profile_reduction_case("arrayfire", "rowSums", x, reps = reps)
  af_colSums <- profile_reduction_case("arrayfire", "colSums", x, reps = reps)

  data.frame(
    family = "reductions",
    n = n,
    backend = rep(c("cpu", "arrayfire"), each = 3),
    op = rep(c("ewise_mul", "rowSums", "colSums"), 2),
    chosen = c(
      cpu_ewise[["chosen"]],
      cpu_rowSums[["chosen"]],
      cpu_colSums[["chosen"]],
      af_ewise[["chosen"]],
      af_rowSums[["chosen"]],
      af_colSums[["chosen"]]
    ),
    elapsed = c(
      as.numeric(cpu_ewise[["elapsed"]]),
      as.numeric(cpu_rowSums[["elapsed"]]),
      as.numeric(cpu_colSums[["elapsed"]]),
      as.numeric(af_ewise[["elapsed"]]),
      as.numeric(af_rowSums[["elapsed"]]),
      as.numeric(af_colSums[["elapsed"]])
    )
  )
}

product_sizes <- c(256L, 512L, 1024L)
reduction_sizes <- c(512L, 1024L, 2048L)

results <- rbind(
  do.call(rbind, lapply(product_sizes, profile_products)),
  do.call(rbind, lapply(reduction_sizes, profile_reductions))
)

print(results)
