# Regression repro for cov2cor(crossprod(adgeMatrix)).
# Seed: none; deterministic literal 3 x 3 matrix
# Shape: 3 x 3 dense input, 3 x 3 crossproduct
# Backend: cpu
# Precision mode: strict
# Dispatch path: crossprod(adgeMatrix) -> cov2cor()
# R/platform: R version 4.5.1 (2025-06-13) | aarch64-apple-darwin20
# Issue: amatrix-af1

test_that("cov2cor accepts crossprod results that preserve adgeMatrix", {
  X_host <- matrix(1:9 + 0.0, 3, 3)
  X <- adgeMatrix(X_host)

  gram <- crossprod(X)
  expect_s4_class(gram, "adgeMatrix")

  got <- cov2cor(gram)
  ref <- stats::cov2cor(crossprod(X_host))

  expect_true(is.matrix(got))
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("cov2cor accepts sparse amatrix covariance-like matrices", {
  S_host <- Matrix::Matrix(
    matrix(c(4, 1, 0, 1, 9, 2, 0, 2, 16), 3, 3),
    sparse = TRUE
  )
  S <- as_adgCMatrix(S_host)

  got <- cov2cor(S)
  ref <- stats::cov2cor(as.matrix(S_host))

  expect_true(is.matrix(got))
  expect_equal(got, ref, tolerance = 1e-12)
})
