# Seed-driven wrapper fuzzing.
#
# Focus:
# - finite dense/sparse cases agree with direct host references
# - non-finite contamination in X or Y is rejected consistently
# - empty / single-row / single-column / zero-column shapes stay well-defined

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(amatrix))

.wrapper_seed_count <- function(default = 8L) {
  raw <- suppressWarnings(as.integer(Sys.getenv("AMATRIX_STRESS_SEEDS", unset = as.character(default))))
  if (is.na(raw) || raw < 1L) {
    default
  } else {
    raw
  }
}

.wrapper_shape_pool <- list(
  c(0L, 3L),
  c(2L, 0L),
  c(1L, 3L),
  c(3L, 1L),
  c(2L, 2L),
  c(4L, 3L)
)

.wrapper_host_matrix <- function(seed, nr, nc) {
  set.seed(seed)
  if (nr == 0L || nc == 0L) {
    return(matrix(numeric(0), nrow = nr, ncol = nc))
  }
  matrix(rnorm(nr * nc), nrow = nr, ncol = nc)
}

.wrapper_host_sparse <- function(seed, nr, nc) {
  set.seed(seed)
  if (nr == 0L || nc == 0L) {
    return(Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(nr, nc)))
  }
  out <- Matrix::rsparsematrix(nr, nc, density = 0.35)
  methods::as(out, "dgCMatrix")
}

.wrapper_dist_ref <- function(X, Y = NULL, method = c("euclidean", "sqeuclidean", "cosine")) {
  method <- match.arg(method)
  Y_eff <- if (is.null(Y)) X else Y

  if (identical(method, "cosine")) {
    G <- tcrossprod(X, Y_eff)
    nx <- sqrt(rowSums(X^2))
    ny <- sqrt(rowSums(Y_eff^2))
    return(G / pmax(outer(nx, ny), .Machine$double.eps))
  }

  G <- tcrossprod(X, Y_eff)
  nx <- rowSums(X^2)
  ny <- rowSums(Y_eff^2)
  D_sq <- pmax(outer(nx, ny, "+") - 2 * G, 0)
  if (is.null(Y)) {
    diag(D_sq) <- 0
  }
  if (!is.null(Y) && identical(X, Y)) {
    diag(D_sq) <- 0
  }
  if (identical(method, "sqeuclidean")) {
    return(D_sq)
  }
  sqrt(D_sq)
}

.wrapper_kernel_ref <- function(X, Y = NULL, kernel = c("linear", "rbf", "cosine"), sigma = 1, zero_diag = FALSE) {
  kernel <- match.arg(kernel)
  Y_eff <- if (is.null(Y)) X else Y
  G <- tcrossprod(X, Y_eff)

  out <- switch(
    kernel,
    linear = G,
    cosine = {
      nx <- sqrt(rowSums(X^2))
      ny <- sqrt(rowSums(Y_eff^2))
      G / pmax(outer(nx, ny), .Machine$double.eps)
    },
    rbf = {
      nx <- rowSums(X^2)
      ny <- rowSums(Y_eff^2)
      D_sq <- pmax(outer(nx, ny, "+") - 2 * G, 0)
      if (is.null(Y)) {
        diag(D_sq) <- 0
      }
      exp(-D_sq / (2 * sigma^2))
    }
  )

  if (is.null(Y) && identical(kernel, "rbf")) {
    diag(out) <- 1
  }
  if (is.null(Y) && isTRUE(zero_diag)) {
    diag(out) <- 0
  }
  out
}

test_that("randomized finite wrapper cases match host references", {
  seed_count <- .wrapper_seed_count(default = 6L)

  for (seed in seq_len(seed_count)) {
    dims <- .wrapper_shape_pool[[((seed - 1L) %% length(.wrapper_shape_pool)) + 1L]]
    nr <- dims[[1L]]
    nc <- dims[[2L]]

    dense_host <- .wrapper_host_matrix(2026042100L + seed, nr, nc)
    sparse_host <- .wrapper_host_sparse(2026042300L + seed, nr, nc)

    cases <- list(
      list(name = "dense", host = dense_host, x = as_adgeMatrix(dense_host)),
      list(name = "sparse", host = sparse_host, x = as_adgCMatrix(sparse_host))
    )

    y_nr <- if (nr == 0L) 0L else max(1L, nr - 1L)
    y_host <- .wrapper_host_matrix(2026042400L + seed, y_nr, nc)

    for (case in cases) {
      info_prefix <- sprintf("seed=%d case=%s dim=%dx%d", seed, case$name, nr, nc)

      expect_equal(
        dist_matrix(case$x, method = "euclidean"),
        .wrapper_dist_ref(as.matrix(case$host), method = "euclidean"),
        tolerance = 1e-10,
        info = paste(info_prefix, "dist self")
      )

      expect_equal(
        dist_matrix(case$x, method = "sqeuclidean"),
        .wrapper_dist_ref(as.matrix(case$host), method = "sqeuclidean"),
        tolerance = 1e-10,
        info = paste(info_prefix, "sqdist self")
      )

      expect_equal(
        kernel_matrix(case$x, kernel = "linear"),
        .wrapper_kernel_ref(as.matrix(case$host), kernel = "linear"),
        tolerance = 1e-10,
        info = paste(info_prefix, "kernel linear")
      )

      expect_equal(
        kernel_matrix(case$x, kernel = "cosine"),
        .wrapper_kernel_ref(as.matrix(case$host), kernel = "cosine"),
        tolerance = 1e-10,
        info = paste(info_prefix, "kernel cosine")
      )

      expect_equal(
        kernel_matrix(case$x, kernel = "rbf", sigma = 0.7),
        .wrapper_kernel_ref(as.matrix(case$host), kernel = "rbf", sigma = 0.7),
        tolerance = 1e-10,
        info = paste(info_prefix, "kernel rbf")
      )

      expect_equal(
        kernel_matrix(case$x, kernel = "rbf", sigma = 0.7, zero_diag = TRUE),
        .wrapper_kernel_ref(as.matrix(case$host), kernel = "rbf", sigma = 0.7, zero_diag = TRUE),
        tolerance = 1e-10,
        info = paste(info_prefix, "kernel rbf zero_diag")
      )

      expect_equal(
        dist_matrix(case$x, y_host, method = "euclidean"),
        .wrapper_dist_ref(as.matrix(case$host), y_host, method = "euclidean"),
        tolerance = 1e-10,
        info = paste(info_prefix, "dist cross")
      )

      expect_equal(
        kernel_matrix(case$x, y_host, kernel = "linear"),
        .wrapper_kernel_ref(as.matrix(case$host), y_host, kernel = "linear"),
        tolerance = 1e-10,
        info = paste(info_prefix, "kernel cross")
      )
    }
  }
})

test_that("randomized non-finite contamination is rejected in X and Y", {
  seed_count <- .wrapper_seed_count(default = 6L)
  bad_values <- c(NA_real_, NaN, Inf, -Inf)

  for (seed in seq_len(seed_count)) {
    dims <- .wrapper_shape_pool[[((seed - 1L) %% length(.wrapper_shape_pool)) + 1L]]
    nr <- max(1L, dims[[1L]])
    nc <- max(1L, dims[[2L]])
    x <- .wrapper_host_matrix(2026042500L + seed, nr, nc)
    y <- .wrapper_host_matrix(2026042600L + seed, max(1L, nr - 1L), nc)

    for (bad in bad_values) {
      x_bad <- x
      x_bad[1L, 1L] <- bad
      y_bad <- y
      y_bad[1L, 1L] <- bad
      info <- sprintf("seed=%d bad=%s", seed, format(bad))

      expect_error(dist_matrix(x_bad), regexp = NULL, info = paste(info, "dist X"))
      expect_error(kernel_matrix(x_bad, kernel = "rbf", sigma = 1), regexp = NULL, info = paste(info, "kernel X"))
      expect_error(dist_matrix(x, y_bad), regexp = NULL, info = paste(info, "dist Y"))
      expect_error(kernel_matrix(x, y_bad, kernel = "linear"), regexp = NULL, info = paste(info, "kernel Y"))
    }
  }
})
