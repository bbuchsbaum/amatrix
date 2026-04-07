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

test_that("qr works on sparse tall matrix", {
  skip_if_not_installed("Matrix")
  set.seed(42)
  n <- 40; p <- 20
  X_dn <- matrix(rnorm(n * p), n, p)
  X_sp <- new_adgCMatrix(Matrix::Matrix(X_dn, sparse = TRUE))
  qr_sp <- qr(X_sp)
  # QR result should be a valid decomposition; verify via solve
  b <- rnorm(n)
  coef_sp <- qr.coef(qr_sp, b)
  coef_ref <- qr.coef(qr(X_dn), b)
  expect_equal(as.numeric(coef_sp), coef_ref, tolerance = 1e-8)
})
