library(testthat)
library(amatrix)

# ── chol_logdet ────────────────────────────────────────────────────────────

test_that("chol_logdet matches log(det(K))", {
  set.seed(1)
  A <- matrix(rnorm(25), 5, 5)
  K <- A %*% t(A) + 5 * diag(5)
  K_am <- adgeMatrix(K)
  L <- chol_factor(K_am)

  expected <- as.double(determinant(K, logarithm = TRUE)$modulus)
  expect_equal(chol_logdet(L), expected, tolerance = 1e-10)
})

test_that("chol_logdet matches 2*sum(log(diag(chol(K))))", {
  set.seed(2)
  A <- matrix(rnorm(36), 6, 6)
  K <- A %*% t(A) + diag(6)
  K_am <- adgeMatrix(K)
  L <- chol_factor(K_am)

  expected <- 2 * sum(log(diag(chol(K))))
  expect_equal(chol_logdet(L), expected, tolerance = 1e-10)
})

test_that("chol_logdet rejects non-amChol", {
  expect_error(chol_logdet(matrix(1:4, 2, 2)), "amChol")
})

# ── chol_diag ─────────────────────────────────────────────────────────────

test_that("chol_diag returns diagonal of upper triangular factor", {
  set.seed(3)
  A <- matrix(rnorm(16), 4, 4)
  K <- A %*% t(A) + 4 * diag(4)
  K_am <- adgeMatrix(K)
  L <- chol_factor(K_am)

  expect_equal(chol_diag(L), diag(chol(K)), tolerance = 1e-10)
})

# ── quad_form ──────────────────────────────────────────────────────────────

test_that("quad_form(L, v) returns v' K^{-1} v (scalar)", {
  set.seed(4)
  A <- matrix(rnorm(25), 5, 5)
  K <- A %*% t(A) + 5 * diag(5)
  v <- rnorm(5)
  K_am <- adgeMatrix(K)
  L <- chol_factor(K_am)

  expected <- as.double(t(v) %*% solve(K, v))
  expect_equal(quad_form(L, v), expected, tolerance = 1e-8)
})

test_that("quad_form(L, V) returns V' K^{-1} V (matrix)", {
  set.seed(5)
  A <- matrix(rnorm(25), 5, 5)
  K <- A %*% t(A) + 5 * diag(5)
  V <- matrix(rnorm(10), 5, 2)
  K_am <- adgeMatrix(K)
  L <- chol_factor(K_am)

  expected <- t(V) %*% solve(K, V)
  expect_equal(quad_form(L, V), expected, tolerance = 1e-8)
})

test_that("quad_form result is positive for nonzero v", {
  set.seed(6)
  A <- matrix(rnorm(16), 4, 4)
  K <- A %*% t(A) + 4 * diag(4)
  v <- rnorm(4)
  L <- chol_factor(adgeMatrix(K))

  expect_gt(quad_form(L, v), 0)
})

# ── eigh ───────────────────────────────────────────────────────────────────

test_that("eigh returns named list with values and vectors", {
  set.seed(7)
  A <- matrix(rnorm(25), 5, 5)
  K <- A %*% t(A) + diag(5)
  res <- eigh(adgeMatrix(K))

  expect_named(res, c("values", "vectors"))
  expect_length(res$values, 5)
  expect_equal(dim(res$vectors), c(5, 5))
})

test_that("eigh eigenvalues are positive for SPD matrix", {
  set.seed(8)
  A <- matrix(rnorm(25), 5, 5)
  K <- A %*% t(A) + 5 * diag(5)
  res <- eigh(adgeMatrix(K))

  expect_true(all(res$values > 0))
})

test_that("eigh reconstruction K ≈ Q diag(lambda) Q'", {
  set.seed(9)
  A <- matrix(rnorm(25), 5, 5)
  K <- A %*% t(A) + 5 * diag(5)
  K_am <- adgeMatrix(K)
  res <- eigh(K_am)

  Q <- res$vectors
  lam <- res$values
  K_recon <- Q %*% diag(lam) %*% t(Q)
  expect_equal(K_recon, K, tolerance = 1e-8)
})

test_that("eigh eigenvectors are orthonormal: Q'Q ≈ I", {
  set.seed(10)
  A <- matrix(rnorm(36), 6, 6)
  K <- A %*% t(A) + diag(6)
  res <- eigh(adgeMatrix(K))

  QtQ <- t(res$vectors) %*% res$vectors
  expect_equal(QtQ, diag(6), tolerance = 1e-8)
})

# ── crossprod_weighted ─────────────────────────────────────────────────────

test_that("crossprod_weighted matches manual X' diag(w) X", {
  set.seed(11)
  X <- matrix(rnorm(30), 10, 3)
  w <- runif(10) + 0.1
  X_am <- adgeMatrix(X)

  expected <- t(X) %*% diag(w) %*% X
  result <- as.matrix(crossprod_weighted(X_am, w))
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("crossprod_weighted result is symmetric", {
  set.seed(12)
  X <- matrix(rnorm(40), 8, 5)
  w <- runif(8) + 0.1
  result <- as.matrix(crossprod_weighted(adgeMatrix(X), w))

  expect_equal(result, t(result), tolerance = 1e-12)
})

test_that("crossprod_weighted errors on length mismatch", {
  X_am <- adgeMatrix(matrix(1:12, 4, 3))
  expect_error(crossprod_weighted(X_am, rep(1, 5)), "nrow")
})

# ── xty_weighted ───────────────────────────────────────────────────────────

test_that("xty_weighted matches manual X' diag(w) y", {
  set.seed(13)
  X <- matrix(rnorm(30), 10, 3)
  y <- rnorm(10)
  w <- runif(10) + 0.1
  X_am <- adgeMatrix(X)

  expected <- t(X) %*% diag(w) %*% y
  result <- as.double(xty_weighted(X_am, w, y))
  expect_equal(result, as.double(expected), tolerance = 1e-10)
})

test_that("xty_weighted handles matrix y", {
  set.seed(14)
  X <- matrix(rnorm(30), 10, 3)
  Y <- matrix(rnorm(20), 10, 2)
  w <- runif(10) + 0.1
  X_am <- adgeMatrix(X)

  expected <- t(X) %*% diag(w) %*% Y
  result <- as.matrix(xty_weighted(X_am, w, Y))
  expect_equal(result, expected, tolerance = 1e-10)
})

# ── rowscale / colscale ─────────────────────────────────────────────────

test_that("rowscale matches diag(d) %*% X", {
  set.seed(15)
  X <- matrix(rnorm(30), 6, 5)
  d <- runif(6) + 0.5
  X_am <- adgeMatrix(X)

  expected <- diag(d) %*% X
  result <- as.matrix(rowscale(X_am, d))
  expect_equal(result, expected, tolerance = 1e-12)
})

test_that("colscale matches X %*% diag(d)", {
  set.seed(16)
  X <- matrix(rnorm(30), 6, 5)
  d <- runif(5) + 0.5
  X_am <- adgeMatrix(X)

  expected <- X %*% diag(d)
  result <- as.matrix(colscale(X_am, d))
  expect_equal(result, expected, tolerance = 1e-12)
})

test_that("rowscale errors on length mismatch", {
  X_am <- adgeMatrix(matrix(1:12, 4, 3))
  expect_error(rowscale(X_am, rep(1, 3)), "nrow")
})

test_that("colscale errors on length mismatch", {
  X_am <- adgeMatrix(matrix(1:12, 4, 3))
  expect_error(colscale(X_am, rep(1, 4)), "ncol")
})

test_that("rowscale result has same dims as input", {
  X_am <- adgeMatrix(matrix(rnorm(20), 4, 5))
  result <- rowscale(X_am, rep(2, 4))
  expect_equal(dim(result), c(4L, 5L))
})

test_that("rowscale followed by crossprod equals crossprod_weighted", {
  set.seed(17)
  X <- matrix(rnorm(30), 10, 3)
  w <- runif(10) + 0.1
  X_am <- adgeMatrix(X)

  via_rowscale  <- as.matrix(crossprod(rowscale(X_am, sqrt(w))))
  via_weighted  <- as.matrix(crossprod_weighted(X_am, w))
  expect_equal(via_rowscale, via_weighted, tolerance = 1e-10)
})

test_that("EMMA rotation pattern: rowscale(X_star, 1/d) then crossprod", {
  # Simulate the EMMA inner loop: X*' D^{-1} X* = crossprod(rowscale(X*, 1/sqrt(d)))
  set.seed(18)
  n <- 20; p <- 3
  X_star <- matrix(rnorm(n * p), n, p)
  lambda  <- sort(abs(rnorm(n))) + 0.1
  delta   <- 0.5
  d       <- lambda + delta

  X_am  <- adgeMatrix(X_star)
  XtDX  <- as.matrix(crossprod_weighted(X_am, 1 / d))
  XtDX2 <- as.matrix(crossprod(rowscale(X_am, 1 / sqrt(d))))
  expect_equal(XtDX, XtDX2, tolerance = 1e-10)
})
