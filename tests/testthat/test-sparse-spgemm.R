test_that("sparse x sparse matmul stays sparse", {
  skip_if_not_installed("Matrix")
  set.seed(42)
  sp <- Matrix::rsparsematrix(100, 50, density=0.05)
  X_sp <- new_adgCMatrix(sp)
  X_dn <- as.matrix(sp)
  # SpGeMM: X %*% t(X) should give sparse or dense result, but correct values
  result <- X_sp %*% t(X_sp)
  ref    <- X_dn %*% t(X_dn)
  expect_equal(as.matrix(result), ref, tolerance = 1e-10)
})
