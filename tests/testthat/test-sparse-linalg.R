test_that("chol works on sparse symmetric positive definite matrix", {
  skip_if_not_installed("Matrix")
  set.seed(42)
  n <- 20
  A_dn <- crossprod(matrix(rnorm(n*n), n, n)) + n * diag(n)
  A_sp <- new_adgCMatrix(Matrix::Matrix(A_dn, sparse=TRUE))
  R <- chol(A_sp)
  # Result should be upper triangular factor; verify A ≈ t(R) %*% R
  R_mat <- if (inherits(R, "adgCMatrix")) as.matrix(amatrix_materialize_host(R)) else as.matrix(R)
  expect_equal(t(R_mat) %*% R_mat, A_dn, tolerance = 1e-8)
})

test_that("solve works on sparse triangular system", {
  skip_if_not_installed("Matrix")
  set.seed(42)
  n <- 20
  b <- rnorm(n)
  A_dn <- crossprod(matrix(rnorm(n*n), n, n)) + n * diag(n)
  A_sp <- new_adgCMatrix(Matrix::Matrix(A_dn, sparse=TRUE))
  x <- solve(A_sp, b)
  expect_equal(as.numeric(A_dn %*% x), b, tolerance = 1e-8)
})
