suppressPackageStartupMessages(library(amatrix))

test_that("kron matches base::kronecker", {
  A <- matrix(1:6 * 1.0, 2, 3)
  B <- matrix(c(1, 0, 0, 2) * 1.0, 2, 2)
  r <- kron(as_adgeMatrix(A), as_adgeMatrix(B))
  expect_s4_class(r, "adgeMatrix")
  expect_equal(as.matrix(r), base::kronecker(A, B), tolerance = 1e-12)
})

test_that("kron dim is nrow(A)*nrow(B) x ncol(A)*ncol(B)", {
  A <- matrix(rnorm(6), 2, 3)
  B <- matrix(rnorm(12), 3, 4)
  r <- kron(as_adgeMatrix(A), as_adgeMatrix(B))
  expect_equal(dim(as.matrix(r)), c(6L, 12L))
})

test_that("kron with identity gives block-diagonal scaling", {
  A <- diag(c(2.0, 3.0))
  B <- diag(2L)
  r <- kron(as_adgeMatrix(A), as_adgeMatrix(B))
  ref <- base::kronecker(A, B)
  expect_equal(as.matrix(r), ref, tolerance = 1e-12)
})
