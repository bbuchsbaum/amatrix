suppressPackageStartupMessages(library(amatrix))

test_that("lu_factor returns amLU for plain matrix", {
  A   <- matrix(c(2, 1, 5, 3), 2, 2)
  fac <- lu_factor(A)
  expect_s4_class(fac, "amLU")
  expect_equal(dim(fac@A), c(2L, 2L))
})

test_that("lu_factor accepts adgeMatrix input", {
  set.seed(1)
  A   <- crossprod(matrix(rnorm(25), 5)) + diag(5)
  fac <- lu_factor(as_adgeMatrix(A))
  expect_s4_class(fac, "amLU")
  expect_identical(nrow(fac@A), 5L)
})

test_that("lu_factor rejects non-square matrix", {
  expect_error(lu_factor(matrix(1:6, 2, 3)), "square")
})

test_that("lu_solve matches base::solve for vector rhs", {
  set.seed(2)
  n <- 8L
  A <- crossprod(matrix(rnorm(n * n), n)) + diag(n)
  b <- rnorm(n)

  fac <- lu_factor(A)
  x   <- lu_solve(fac, b)
  ref <- base::solve(A, b)

  expect_equal(x, ref, tolerance = 1e-10)
  expect_type(x, "double")
  expect_length(x, n)
})

test_that("lu_solve handles matrix rhs", {
  set.seed(3)
  n <- 6L; k <- 3L
  A <- crossprod(matrix(rnorm(n * n), n)) + diag(n)
  B <- matrix(rnorm(n * k), n, k)

  fac <- lu_factor(A)
  X   <- lu_solve(fac, B)
  ref <- base::solve(A, B)

  expect_equal(X, ref, tolerance = 1e-10)
  expect_equal(dim(X), c(n, k))
})

test_that("lu_solve errors on wrong factor type", {
  expect_error(lu_solve(list(), rnorm(3)), "amLU")
})

test_that("lu_factor works on asymmetric (non-SPD) matrix", {
  set.seed(5)
  n <- 5L
  # General (non-symmetric) invertible matrix
  A <- matrix(rnorm(n * n), n, n) + diag(n) * n
  b <- rnorm(n)

  fac <- lu_factor(A)
  x   <- lu_solve(fac, b)
  ref <- base::solve(A, b)
  expect_equal(x, ref, tolerance = 1e-10)
})
