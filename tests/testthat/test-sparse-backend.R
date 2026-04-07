test_that("adgCMatrix preserves preferred_backend slot", {
  skip_if_not_installed("Matrix")
  sp <- Matrix::rsparsematrix(50, 30, density=0.1)
  X <- new_adgCMatrix(sp)
  # Default preferred_backend should be accessible
  expect_true(!is.null(X@preferred_backend) || is.character(X@preferred_backend))
})

test_that("adgCMatrix preferred_backend defaults to cpu", {
  skip_if_not_installed("Matrix")
  sp <- Matrix::rsparsematrix(50, 30, density=0.1)
  X <- new_adgCMatrix(sp)
  expect_equal(X@preferred_backend, "cpu")
})

test_that("adgCMatrix constructor respects custom preferred_backend", {
  skip_if_not_installed("Matrix")
  sp <- Matrix::rsparsematrix(50, 30, density=0.1)
  X <- new_adgCMatrix(sp, preferred_backend = "cpu")
  expect_equal(X@preferred_backend, "cpu")
})
