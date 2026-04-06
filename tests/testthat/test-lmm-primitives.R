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

# ── rowmeans / colmeans ───────────────────────────────────────────────────────

test_that("rowmeans matches base rowMeans", {
  set.seed(19)
  X <- matrix(rnorm(30), 6, 5)
  expect_equal(rowmeans(adgeMatrix(X)), rowMeans(X), tolerance = 1e-12)
})

test_that("colmeans matches base colMeans", {
  set.seed(20)
  X <- matrix(rnorm(30), 6, 5)
  expect_equal(colmeans(adgeMatrix(X)), colMeans(X), tolerance = 1e-12)
})

test_that("rowMeans S4 dispatch works on adgeMatrix", {
  set.seed(21)
  X <- matrix(rnorm(20), 4, 5)
  expect_equal(rowMeans(adgeMatrix(X)), rowMeans(X), tolerance = 1e-12)
})

test_that("colMeans S4 dispatch works on adgeMatrix", {
  set.seed(22)
  X <- matrix(rnorm(20), 4, 5)
  expect_equal(colMeans(adgeMatrix(X)), colMeans(X), tolerance = 1e-12)
})

# ── trace ─────────────────────────────────────────────────────────────────────

test_that("trace matches sum(diag(A)) for square matrix", {
  set.seed(23)
  A <- matrix(rnorm(25), 5, 5)
  expect_equal(trace(adgeMatrix(A)), sum(diag(A)), tolerance = 1e-12)
})

test_that("trace of identity matrix equals n", {
  expect_equal(trace(adgeMatrix(diag(7))), 7, tolerance = 1e-12)
})

test_that("trace of SPD matrix equals sum of eigenvalues", {
  set.seed(24)
  A <- matrix(rnorm(25), 5, 5)
  K <- A %*% t(A) + diag(5)
  expect_equal(trace(adgeMatrix(K)), sum(eigen(K, only.values = TRUE)$values),
               tolerance = 1e-8)
})

# ── sym ───────────────────────────────────────────────────────────────────────

test_that("sym result is exactly symmetric", {
  set.seed(25)
  A <- matrix(rnorm(25), 5, 5)
  S <- sym(adgeMatrix(A))
  expect_equal(as.matrix(S), t(as.matrix(S)), tolerance = 1e-14)
})

test_that("sym matches manual (A + t(A)) / 2", {
  set.seed(26)
  A <- matrix(rnorm(16), 4, 4)
  expected <- (A + t(A)) / 2
  expect_equal(as.matrix(sym(adgeMatrix(A))), expected, tolerance = 1e-12)
})

test_that("sym is idempotent", {
  set.seed(27)
  A <- matrix(rnorm(16), 4, 4)
  A_am <- adgeMatrix(A)
  expect_equal(as.matrix(sym(sym(A_am))), as.matrix(sym(A_am)), tolerance = 1e-14)
})

# ── dot ───────────────────────────────────────────────────────────────────────

test_that("dot(x, y) matches sum(x * y) for vectors", {
  set.seed(28)
  x <- rnorm(10); y <- rnorm(10)
  expect_equal(dot(x, y), sum(x * y), tolerance = 1e-12)
})

test_that("dot(X, Y) matches sum(X * Y) for matrices (Frobenius inner product)", {
  set.seed(29)
  X <- matrix(rnorm(20), 4, 5)
  Y <- matrix(rnorm(20), 4, 5)
  expect_equal(dot(adgeMatrix(X), adgeMatrix(Y)), sum(X * Y), tolerance = 1e-12)
})

test_that("dot(x, x) equals squared norm", {
  set.seed(30)
  x <- rnorm(8)
  expect_equal(dot(x, x), sum(x^2), tolerance = 1e-12)
})

# ── crossprod_add_diag ────────────────────────────────────────────────────────

test_that("crossprod_add_diag(X, lambda) matches crossprod(X) + lambda*I", {
  set.seed(31)
  X <- matrix(rnorm(30), 6, 5)
  lambda <- 0.5
  expected <- t(X) %*% X + lambda * diag(5)
  result <- as.matrix(crossprod_add_diag(adgeMatrix(X), lambda))
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("crossprod_add_diag(X, d) matches crossprod(X) + diag(d)", {
  set.seed(32)
  X <- matrix(rnorm(30), 6, 5)
  d <- runif(5) + 0.1
  expected <- t(X) %*% X + diag(d)
  result <- as.matrix(crossprod_add_diag(adgeMatrix(X), d))
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("crossprod_add_diag result is symmetric", {
  set.seed(33)
  X <- matrix(rnorm(24), 8, 3)
  R <- as.matrix(crossprod_add_diag(adgeMatrix(X), 1.0))
  expect_equal(R, t(R), tolerance = 1e-12)
})

# ── mat_sqrt / mat_pow / mat_log ──────────────────────────────────────────────

test_that("mat_sqrt(K) %*% mat_sqrt(K) ≈ K", {
  set.seed(34)
  A <- matrix(rnorm(25), 5, 5)
  K <- crossprod(A) + 5 * diag(5)
  K_am <- adgeMatrix(K)
  S <- as.matrix(mat_sqrt(K_am))
  expect_equal(S %*% S, K, tolerance = 1e-8)
})

test_that("mat_pow(K, 2) ≈ K %*% K", {
  set.seed(35)
  A <- matrix(rnorm(16), 4, 4)
  K <- crossprod(A) + 4 * diag(4)
  K_am <- adgeMatrix(K)
  expect_equal(as.matrix(mat_pow(K_am, 2)), K %*% K, tolerance = 1e-8)
})

test_that("mat_pow(K, 0.5) matches mat_sqrt(K)", {
  set.seed(36)
  A <- matrix(rnorm(16), 4, 4)
  K <- crossprod(A) + 4 * diag(4)
  K_am <- adgeMatrix(K)
  expect_equal(as.matrix(mat_pow(K_am, 0.5)), as.matrix(mat_sqrt(K_am)),
               tolerance = 1e-8)
})

test_that("mat_log roundtrip: exp(mat_log(K)) ≈ K via eigenvalues", {
  set.seed(37)
  A <- matrix(rnorm(16), 4, 4)
  K <- crossprod(A) + 4 * diag(4)
  K_am <- adgeMatrix(K)
  logK <- mat_log(K_am)
  res <- eigh(logK)
  expect_equal(exp(res$values), eigh(K_am)$values, tolerance = 1e-6)
})

# ── solve_triangular ──────────────────────────────────────────────────────────

test_that("solve_triangular upper recovers x from R %*% x", {
  set.seed(38)
  A <- matrix(rnorm(25), 5, 5)
  K <- crossprod(A) + 5 * diag(5)
  R <- chol(K)
  x <- rnorm(5)
  expect_equal(solve_triangular(R, R %*% x), x, tolerance = 1e-10)
})

test_that("solve_triangular lower recovers x from L %*% x", {
  set.seed(39)
  A <- matrix(rnorm(25), 5, 5)
  K <- crossprod(A) + 5 * diag(5)
  L <- t(chol(K))
  x <- rnorm(5)
  expect_equal(solve_triangular(L, L %*% x, lower = TRUE), x, tolerance = 1e-10)
})

test_that("solve_triangular handles multiple RHS", {
  set.seed(40)
  A <- matrix(rnorm(25), 5, 5)
  K <- crossprod(A) + 5 * diag(5)
  R <- chol(K)
  X <- matrix(rnorm(15), 5, 3)
  expect_equal(solve_triangular(R, R %*% X), X, tolerance = 1e-10)
})

# ── trace_estim ───────────────────────────────────────────────────────────────

test_that("trace_estim converges to true trace with large k", {
  set.seed(41)
  A <- matrix(rnorm(25), 5, 5)
  K <- crossprod(A) + 5 * diag(5)
  K_am <- adgeMatrix(K)
  est <- trace_estim(K_am, k = 500, seed = 42)
  expect_equal(est, trace(K_am), tolerance = 1.0)  # stochastic, generous tol
})

test_that("trace_estim via solve_fn estimates tr(K^{-1})", {
  set.seed(42)
  A <- matrix(rnorm(16), 4, 4)
  K <- crossprod(A) + 4 * diag(4)
  L <- chol_factor(adgeMatrix(K))
  true_val <- sum(diag(solve(K)))
  est <- trace_estim(solve_fn = function(v) chol_solve(L, v), n = 4L,
                     k = 500L, seed = 1L)
  expect_equal(est, true_val, tolerance = 0.5)  # stochastic
})

test_that("trace_estim is reproducible with seed", {
  set.seed(43)
  K <- crossprod(matrix(rnorm(25), 5, 5)) + 5 * diag(5)
  K_am <- adgeMatrix(K)
  e1 <- trace_estim(K_am, k = 20, seed = 99)
  e2 <- trace_estim(K_am, k = 20, seed = 99)
  expect_equal(e1, e2)
})

# ---------------------------------------------------------------------------
# batch_chol / batch_solve / batch_crossprod
# ---------------------------------------------------------------------------

test_that("batch_chol returns a list of amChol objects", {
  set.seed(1)
  B <- 4L; n <- 5L
  mats <- lapply(seq_len(B), function(i) {
    A <- matrix(rnorm(n * n), n, n); crossprod(A) + diag(n)
  })
  Ls <- batch_chol(mats)
  expect_length(Ls, B)
  expect_true(all(vapply(Ls, inherits, logical(1), "amChol")))
})

test_that("batch_chol accepts a 3-D array", {
  set.seed(2)
  B <- 3L; n <- 4L
  arr <- array(0, c(n, n, B))
  for (b in seq_len(B)) { A <- matrix(rnorm(n*n), n, n); arr[,,b] <- crossprod(A) + diag(n) }
  Ls <- batch_chol(arr)
  expect_length(Ls, B)
})

test_that("batch_solve recovers RHS", {
  set.seed(3)
  B <- 5L; n <- 6L; k <- 2L
  mats <- lapply(seq_len(B), function(i) { A <- matrix(rnorm(n*n), n, n); crossprod(A) + diag(n) })
  Xs   <- lapply(seq_len(B), function(i) matrix(rnorm(n * k), n, k))
  Ls   <- batch_chol(mats)
  # compute RHS = A_b %*% x_b for each b
  rhs  <- Map(function(m, x) m %*% x, mats, Xs)
  sols <- batch_solve(Ls, rhs)
  for (b in seq_len(B)) {
    expect_equal(sols[[b]], Xs[[b]], tolerance = 1e-9)
  }
})

test_that("batch_solve works with 3-D B array", {
  set.seed(4)
  B <- 3L; n <- 4L; k <- 2L
  mats <- lapply(seq_len(B), function(i) { A <- matrix(rnorm(n*n), n, n); crossprod(A) + diag(n) })
  Ls   <- batch_chol(mats)
  rhs_arr <- array(rnorm(n * k * B), c(n, k, B))
  sols <- batch_solve(Ls, rhs_arr)
  expect_length(sols, B)
})

test_that("batch_crossprod returns list of p x p matrices", {
  set.seed(5)
  B <- 4L; n <- 8L; p <- 3L
  mats <- lapply(seq_len(B), function(i) matrix(rnorm(n * p), n, p))
  XtXs <- batch_crossprod(mats)
  expect_length(XtXs, B)
  for (b in seq_len(B)) {
    expect_equal(dim(XtXs[[b]]), c(p, p))
    expect_equal(XtXs[[b]], crossprod(mats[[b]]), tolerance = 1e-12)
  }
})
