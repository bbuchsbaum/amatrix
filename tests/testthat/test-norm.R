suppressPackageStartupMessages(library(amatrix))

# ---- matrix norms -----------------------------------------------------------

test_that("norm Frobenius matches base::norm", {
  set.seed(1)
  A  <- matrix(rnorm(12), 3, 4)
  Xam <- as_adgeMatrix(A)
  expect_equal(norm(Xam, "F"), base::norm(A, "F"), tolerance = 1e-12)
})

test_that("norm 1-norm (max col sum) matches base::norm", {
  set.seed(2)
  A  <- matrix(rnorm(12), 3, 4)
  Xam <- as_adgeMatrix(A)
  expect_equal(norm(Xam, "1"), base::norm(A, "1"), tolerance = 1e-12)
})

test_that("norm infinity norm (max row sum) matches base::norm", {
  set.seed(3)
  A  <- matrix(rnorm(12), 3, 4)
  Xam <- as_adgeMatrix(A)
  expect_equal(norm(Xam, "I"), base::norm(A, "I"), tolerance = 1e-12)
})

test_that("norm max-entry norm matches base::norm", {
  set.seed(4)
  A  <- matrix(rnorm(12), 3, 4)
  Xam <- as_adgeMatrix(A)
  expect_equal(norm(Xam, "M"), base::norm(A, "M"), tolerance = 1e-12)
})

test_that("norm spectral norm (type='2') matches base::svd largest value", {
  set.seed(5)
  A  <- matrix(rnorm(15), 3, 5)
  Xam <- as_adgeMatrix(A)
  ref_sv1 <- base::svd(A, nu = 0, nv = 0)$d[[1L]]
  expect_equal(norm(Xam, "2"), ref_sv1, tolerance = 1e-6)
})

test_that("norm accepts plain matrix without wrapping", {
  A <- matrix(c(3, 4, 0, 0) * 1.0, 2, 2)
  expect_equal(norm(A, "F"), 5.0, tolerance = 1e-12)
})

# ---- vector norms -----------------------------------------------------------

test_that("norm Euclidean (2-norm) for vector", {
  x <- c(3.0, 4.0)
  expect_equal(norm(x, "2"), 5.0, tolerance = 1e-12)
})

test_that("norm 1-norm for vector", {
  x <- c(-3.0, 4.0, -5.0)
  expect_equal(norm(x, "1"), 12.0, tolerance = 1e-12)
})

test_that("norm infinity norm for vector", {
  x <- c(1.0, -7.0, 3.0)
  expect_equal(norm(x, "I"), 7.0, tolerance = 1e-12)
})
