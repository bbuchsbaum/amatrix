# Tests for four API-gap fixes:
#   amatrix-jnd  kronecker() / %x% preserve the amatrix class
#   amatrix-x6a  KronMatrix gains a `[` method (materialize-on-subset)
#   amatrix-sxs  adgeMatrix gains a `[[` scalar-extraction method
#   amatrix-vbh  amatrix_release_resident() is implemented and exported
#
# All tests are CPU-only and set no GPU environment variables.

# ---------------------------------------------------------------------------
# amatrix-jnd: kronecker / %x% class preservation
# ---------------------------------------------------------------------------

test_that("kronecker(adge, adge) stays adgeMatrix and matches base::kronecker", {
  A <- matrix(c(1, 2, 3, 4) * 1.0, 2, 2)
  B <- matrix(c(5, 6, 7, 8, 9, 10) * 1.0, 2, 3)
  aA <- adgeMatrix(A)
  aB <- adgeMatrix(B)

  r <- kronecker(aA, aB)
  expect_s4_class(r, "adgeMatrix")
  expect_equal(as.matrix(r), base::kronecker(A, B), tolerance = 1e-12)
})

test_that("A %x% B preserves the amatrix class", {
  A <- matrix(c(1, 2, 3, 4) * 1.0, 2, 2)
  B <- diag(c(2.0, 3.0))
  aA <- adgeMatrix(A)
  aB <- adgeMatrix(B)

  r <- aA %x% aB
  expect_s4_class(r, "adgeMatrix")
  expect_equal(as.matrix(r), base::kronecker(A, B), tolerance = 1e-12)
})

test_that("kronecker(adgC, adgC) stays adgCMatrix and matches base::kronecker", {
  A <- matrix(c(1, 0, 0, 2) * 1.0, 2, 2)
  B <- matrix(c(0, 3, 4, 0) * 1.0, 2, 2)
  sA <- adgCMatrix(A)
  sB <- adgCMatrix(B)

  r <- kronecker(sA, sB)
  expect_s4_class(r, "adgCMatrix")
  expect_equal(as.matrix(r), base::kronecker(A, B), tolerance = 1e-12)

  r2 <- sA %x% sB
  expect_s4_class(r2, "adgCMatrix")
  expect_equal(as.matrix(r2), base::kronecker(A, B), tolerance = 1e-12)
})

test_that("kronecker with a base matrix operand stays amatrix", {
  A <- matrix(c(1, 2, 3, 4) * 1.0, 2, 2)
  B <- matrix(c(5, 6, 7, 8) * 1.0, 2, 2)
  aA <- adgeMatrix(A)

  r_right <- kronecker(aA, B)
  expect_s4_class(r_right, "adgeMatrix")
  expect_equal(as.matrix(r_right), base::kronecker(A, B), tolerance = 1e-12)

  r_left <- kronecker(A, aA)
  expect_s4_class(r_left, "adgeMatrix")
  expect_equal(as.matrix(r_left), base::kronecker(A, A), tolerance = 1e-12)
})

test_that("kronecker preserves backend metadata of the first amatrix operand", {
  aA <- adgeMatrix(matrix(1:4 * 1.0, 2, 2), preferred_backend = "cpu", precision = "fast")
  aB <- adgeMatrix(diag(2))

  r <- kronecker(aA, aB)
  expect_identical(r@preferred_backend, aA@preferred_backend)
  expect_identical(r@precision, aA@precision)
})

# ---------------------------------------------------------------------------
# amatrix-x6a: KronMatrix `[`
# ---------------------------------------------------------------------------

test_that("KronMatrix supports 2-D subsetting matching the dense product", {
  A <- matrix(c(1, 2, 3, 4) * 1.0, 2, 2)
  B <- matrix(c(5, 6, 7, 8, 9, 10) * 1.0, 2, 3)
  K <- kron_matrix(A, B)
  full <- base::kronecker(A, B)

  expect_equal(K[2, 3], full[2, 3])
  expect_equal(K[1, ], full[1, ])
  expect_equal(K[, 4], full[, 4])
  expect_equal(K[c(1, 3), c(2, 5)], full[c(1, 3), c(2, 5)])
})

test_that("KronMatrix subsetting honours drop and linear indexing", {
  A <- matrix(c(1, 2, 3, 4) * 1.0, 2, 2)
  B <- diag(2)
  K <- kron_matrix(A, B)
  full <- base::kronecker(A, B)

  expect_equal(K[1, , drop = FALSE], full[1, , drop = FALSE])
  expect_equal(dim(K[1, , drop = FALSE]), c(1L, 4L))
  # linear (single-subscript) indexing
  expect_equal(K[5], full[5])
  expect_equal(K[c(1, 6, 11)], full[c(1, 6, 11)])
  # whole-matrix K[]
  expect_equal(K[], full)
})

test_that("KronMatrix subsetting no longer errors as unsubsettable", {
  K <- kron_matrix(matrix(1:4 * 1.0, 2, 2), diag(2))
  expect_error(K[1, 1], NA)
})
