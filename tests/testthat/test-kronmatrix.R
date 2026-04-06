suppressPackageStartupMessages(library(amatrix))

# helper: materialize reference via base::kronecker
kref <- function(A, B) base::kronecker(A, B)

set.seed(42)

# ---- construction & shape --------------------------------------------------

test_that("kron_matrix() returns a KronMatrix with correct dim", {
  A <- matrix(rnorm(6), 2, 3)
  B <- matrix(rnorm(12), 4, 3)
  K <- kron_matrix(A, B)
  expect_s4_class(K, "KronMatrix")
  expect_equal(dim(K), c(nrow(A) * nrow(B), ncol(A) * ncol(B)))
  expect_equal(nrow(K), nrow(A) * nrow(B))
  expect_equal(ncol(K), ncol(A) * ncol(B))
})

# ---- as.matrix (materialization) ------------------------------------------

test_that("as.matrix matches base::kronecker", {
  A <- matrix(rnorm(4), 2, 2)
  B <- matrix(rnorm(9), 3, 3)
  K <- kron_matrix(A, B)
  expect_equal(as.matrix(K), kref(A, B), tolerance = 1e-12)
})

# ---- transpose -------------------------------------------------------------

test_that("t(K) = A^T ⊗ B^T", {
  A <- matrix(rnorm(6), 2, 3)
  B <- matrix(rnorm(12), 4, 3)
  K  <- kron_matrix(A, B)
  Kt <- t(K)
  expect_equal(as.matrix(Kt), t(kref(A, B)), tolerance = 1e-12)
})

# ---- matrix-vector multiply ------------------------------------------------

test_that("K %*% y matches reference", {
  A <- matrix(rnorm(4), 2, 2)
  B <- matrix(rnorm(9), 3, 3)
  K <- kron_matrix(A, B)
  y <- rnorm(ncol(K))
  expect_equal(as.vector(K %*% y), as.vector(kref(A, B) %*% y), tolerance = 1e-12)
})

test_that("K %*% Y (matrix RHS) matches reference", {
  A <- matrix(rnorm(4), 2, 2)
  B <- matrix(rnorm(9), 3, 3)
  K <- kron_matrix(A, B)
  Y <- matrix(rnorm(ncol(K) * 4), ncol(K), 4)
  expect_equal(K %*% Y, kref(A, B) %*% Y, tolerance = 1e-12)
})

# ---- left multiply ---------------------------------------------------------

test_that("y %*% K matches reference", {
  A <- matrix(rnorm(4), 2, 2)
  B <- matrix(rnorm(9), 3, 3)
  K <- kron_matrix(A, B)
  y <- rnorm(nrow(K))
  expect_equal(as.vector(y %*% K), as.vector(y %*% kref(A, B)), tolerance = 1e-12)
})

test_that("X %*% K (matrix LHS) matches reference", {
  A <- matrix(rnorm(4), 2, 2)
  B <- matrix(rnorm(9), 3, 3)
  K <- kron_matrix(A, B)
  X <- matrix(rnorm(3 * nrow(K)), 3, nrow(K))
  expect_equal(X %*% K, X %*% kref(A, B), tolerance = 1e-12)
})

# ---- crossprod -------------------------------------------------------------

test_that("crossprod(K) = (A^T A) ⊗ (B^T B)", {
  A <- matrix(rnorm(6), 3, 2)
  B <- matrix(rnorm(12), 4, 3)
  K  <- kron_matrix(A, B)
  Cp <- crossprod(K)
  expect_s4_class(Cp, "KronMatrix")
  expect_equal(as.matrix(Cp), crossprod(kref(A, B)), tolerance = 1e-12)
})

test_that("crossprod(K, y) matches t(K) %*% y", {
  A <- matrix(rnorm(4), 2, 2)
  B <- matrix(rnorm(9), 3, 3)
  K <- kron_matrix(A, B)
  y <- rnorm(nrow(K))
  expect_equal(as.vector(crossprod(K, y)),
               as.vector(t(kref(A, B)) %*% y), tolerance = 1e-12)
})

# ---- solve -----------------------------------------------------------------

test_that("solve(K) gives A^{-1} ⊗ B^{-1}", {
  A <- matrix(c(2, 1, 0, 3), 2, 2)
  B <- matrix(c(4, 1, 2, 3), 2, 2)
  K    <- kron_matrix(A, B)
  Kinv <- solve(K)
  expect_s4_class(Kinv, "KronMatrix")
  expect_equal(as.matrix(Kinv), solve(kref(A, B)), tolerance = 1e-10)
})

test_that("solve(K, b) gives K^{-1} b", {
  A <- matrix(c(2, 1, 0, 3), 2, 2)
  B <- matrix(c(4, 1, 2, 3), 2, 2)
  K <- kron_matrix(A, B)
  b <- rnorm(nrow(K))
  expect_equal(solve(K, b), as.vector(solve(kref(A, B), b)), tolerance = 1e-10)
})

# ---- determinant -----------------------------------------------------------

test_that("determinant matches base on square factors", {
  A <- matrix(c(2, 1, 0, 3), 2, 2)
  B <- matrix(c(4, 1, 2, 3), 2, 2)
  K <- kron_matrix(A, B)
  d_kron <- determinant(K, logarithm = TRUE)
  d_ref  <- determinant(kref(A, B), logarithm = TRUE)
  expect_equal(as.double(d_kron$modulus), as.double(d_ref$modulus), tolerance = 1e-10)
  expect_equal(d_kron$sign, d_ref$sign)
})
