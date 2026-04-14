suppressPackageStartupMessages(library(amatrix))

# Tests covering S4 dispatch gotchas that were fixed.

set.seed(99)

# ---- eigen: symmetric arg is optional (mirrors base::eigen) -----------------

test_that("eigen(adgeMatrix) without symmetric auto-detects symmetry", {
  S <- matrix(c(3, 1, 1, 2) * 1.0, 2, 2)
  Am <- adgeMatrix(S)
  e_auto <- eigen(Am)               # no symmetric= argument
  e_sym  <- eigen(Am, symmetric = TRUE)
  e_base <- base::eigen(S)
  expect_equal(e_auto$values, e_base$values, tolerance = 1e-10)
  expect_equal(e_sym$values,  e_base$values, tolerance = 1e-10)
})

test_that("eigen(adgeMatrix, symmetric=FALSE) uses general path", {
  A  <- matrix(c(2, 1, 0, 3) * 1.0, 2, 2)  # upper triangular, not symmetric
  Am <- adgeMatrix(A)
  e_auto <- eigen(Am)               # auto-detects non-symmetric
  e_base <- base::eigen(A)
  expect_equal(sort(Re(e_auto$values)), sort(Re(e_base$values)), tolerance = 1e-10)
})

test_that("eigen(adgCMatrix) without symmetric auto-detects", {
  S   <- Matrix::Matrix(crossprod(matrix(rnorm(9), 3, 3)) + diag(3), sparse = TRUE)
  Asp <- adgCMatrix(S)
  e_auto <- eigen(Asp)
  e_base <- base::eigen(as.matrix(S), symmetric = TRUE)
  expect_equal(e_auto$values, e_base$values, tolerance = 1e-10)
})

# ---- solve: vector RHS returns plain numeric vector -------------------------

test_that("solve(adgeMatrix, numeric_vector) returns numeric not adgeMatrix", {
  A <- adgeMatrix(matrix(c(2, 1, 1, 3) * 1.0, 2, 2))
  b <- c(1.0, 2.0)
  r <- solve(A, b)
  expect_type(r, "double")           # plain numeric
  expect_null(dim(r))                # no dim → vector, not matrix
  expect_equal(r, base::solve(as.matrix(A), b), tolerance = 1e-10)
})

test_that("solve(adgeMatrix, matrix) returns adgeMatrix", {
  A <- adgeMatrix(matrix(c(2, 1, 1, 3) * 1.0, 2, 2))
  B <- matrix(c(1, 0, 0, 1) * 1.0, 2, 2)
  r <- solve(A, B)
  expect_s4_class(r, "adgeMatrix")
  expect_equal(as.matrix(r), base::solve(as.matrix(A), B), tolerance = 1e-10)
})

test_that("solve(adgeMatrix) returns adgeMatrix inverse", {
  A <- adgeMatrix(matrix(c(2, 1, 1, 3) * 1.0, 2, 2))
  r <- solve(A)
  expect_s4_class(r, "adgeMatrix")
  expect_equal(as.matrix(r), base::solve(as.matrix(A)), tolerance = 1e-10)
})

# ---- diag: extract mode vs create mode --------------------------------------

test_that("diag(adgeMatrix) extracts diagonal as numeric vector", {
  set.seed(1)
  M  <- crossprod(matrix(rnorm(16), 4, 4)) + diag(4)
  Am <- adgeMatrix(M)
  d  <- diag(Am)
  expect_type(d, "double")
  expect_length(d, 4L)
  expect_equal(d, base::diag(M), tolerance = 1e-12)
})

test_that("diag(adgCMatrix) extracts diagonal as numeric vector", {
  M   <- Matrix::Matrix(crossprod(matrix(rnorm(9), 3, 3)) + diag(3), sparse = TRUE)
  Asp <- adgCMatrix(M)
  d   <- diag(Asp)
  expect_type(d, "double")
  expect_equal(d, base::diag(as.matrix(M)), tolerance = 1e-12)
})

test_that("Math(adgeMatrix) preserves adgeMatrix class", {
  A <- adgeMatrix(matrix(c(1, 4, 9, 16) * 1.0, 2, 2))
  out <- sqrt(A)
  expect_s4_class(out, "adgeMatrix")
  expect_equal(as.matrix(out), sqrt(as.matrix(A)), tolerance = 1e-12)
})

test_that("Math(adgCMatrix) preserves adgCMatrix class when result stays sparse", {
  S <- adgCMatrix(as(Matrix::Diagonal(x = c(1, 4, 9)), "dgCMatrix"))
  out <- sqrt(S)
  expect_s4_class(out, "adgCMatrix")
  expect_equal(as.matrix(out), sqrt(as.matrix(S)), tolerance = 1e-12)
})

test_that("diag<-(adgeMatrix) preserves adgeMatrix class", {
  A <- adgeMatrix(matrix(0, 3, 3))
  diag(A) <- c(1, 2, 3)
  expect_s4_class(A, "adgeMatrix")
  expect_equal(as.matrix(A), diag(c(1, 2, 3)), tolerance = 1e-12)
})

test_that("diag<-(adgCMatrix) preserves adgCMatrix class", {
  S <- adgCMatrix(as(Matrix::Diagonal(3), "dgCMatrix"))
  diag(S) <- c(2, 4, 6)
  expect_s4_class(S, "adgCMatrix")
  expect_equal(as.matrix(S), diag(c(2, 4, 6)), tolerance = 1e-12)
})

# ---- Ops: adgeMatrix + adgeMatrix returns adgeMatrix (not dgeMatrix) --------

test_that("adgeMatrix + adgeMatrix returns adgeMatrix", {
  A <- adgeMatrix(matrix(1:6 * 1.0, 2, 3))
  r <- A + A
  expect_s4_class(r, "adgeMatrix")
  expect_equal(as.matrix(r), as.matrix(A) + as.matrix(A), tolerance = 1e-12)
})

test_that("adgeMatrix * adgeMatrix returns adgeMatrix", {
  A <- adgeMatrix(matrix(1:6 * 1.0, 2, 3))
  r <- A * A
  expect_s4_class(r, "adgeMatrix")
  expect_equal(as.matrix(r), as.matrix(A) * as.matrix(A), tolerance = 1e-12)
})

test_that("adgeMatrix - adgeMatrix returns adgeMatrix", {
  A <- adgeMatrix(matrix(1:6 * 1.0, 2, 3))
  r <- A - A
  expect_s4_class(r, "adgeMatrix")
  expect_true(all(as.matrix(r) == 0))
})

test_that("adgCMatrix + adgCMatrix returns adgCMatrix", {
  M  <- Matrix::Matrix(diag(3), sparse = TRUE)
  Sp <- adgCMatrix(M)
  r  <- Sp + Sp
  expect_s4_class(r, "adgCMatrix")
})
