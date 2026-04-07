test_that("adgCMatrix arithmetic ops work correctly", {
  skip_if_not_installed("Matrix")
  set.seed(42)
  sp <- Matrix::rsparsematrix(100, 50, density = 0.05)
  X_sp <- new_adgCMatrix(sp)
  X_dn <- as.matrix(sp)

  # scalar ops
  expect_equal(as.matrix(X_sp + 1), X_dn + 1, tolerance = 1e-10)
  expect_equal(as.matrix(X_sp * 2), X_dn * 2, tolerance = 1e-10)
  expect_equal(as.matrix(X_sp / 2), X_dn / 2, tolerance = 1e-10)
  expect_equal(as.matrix(X_sp - 1), X_dn - 1, tolerance = 1e-10)
})

test_that("adgCMatrix self-subtraction gives zero", {
  skip_if_not_installed("Matrix")
  set.seed(42)
  sp <- Matrix::rsparsematrix(100, 50, density = 0.05)
  X_sp <- new_adgCMatrix(sp)
  result <- X_sp - X_sp
  expect_true(all(as.matrix(result) == 0))
})

test_that("adgCMatrix comparison ops work", {
  skip_if_not_installed("Matrix")
  set.seed(42)
  sp <- Matrix::rsparsematrix(100, 50, density = 0.05)
  X_sp <- new_adgCMatrix(sp)
  X_dn <- as.matrix(sp)

  expect_equal(as.matrix(X_sp > 0), X_dn > 0)
  expect_equal(as.matrix(X_sp == 0), X_dn == 0)
})

test_that("adgCMatrix ops with same-class rhs work", {
  skip_if_not_installed("Matrix")
  set.seed(42)
  sp <- Matrix::rsparsematrix(100, 50, density = 0.05)
  X_sp <- new_adgCMatrix(sp)
  X_dn <- as.matrix(sp)

  result <- X_sp * X_sp
  expect_equal(as.matrix(result), X_dn * X_dn, tolerance = 1e-10)
})
