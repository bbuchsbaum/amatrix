# test-lmm-building-blocks.R
# Tests for am_chol_logdet, am_chol_diag, am_quad_form (amatrix-rrv),
# am_eigh (amatrix-bg2), and weighted crossprod helpers (amatrix-sdl).

# ── Fixtures ──────────────────────────────────────────────────────────────────

.lmm_spd <- local({
  set.seed(99)
  p <- 6L
  Z  <- matrix(rnorm(p * p), p, p)
  S  <- crossprod(Z) + diag(p) * 3
  ch <- chol(S)
  list(S = S, ch = ch, p = p)
})

# ── am_chol_logdet ────────────────────────────────────────────────────────────

test_that("am_chol_logdet equals log(det(S))", {
  d  <- .lmm_spd
  fac <- am_chol_factor(adgeMatrix(d$S, preferred_backend = "cpu"))
  expect_equal(am_chol_logdet(fac),
               as.numeric(determinant(d$S, logarithm = TRUE)$modulus),
               tolerance = 1e-10)
})

test_that("am_chol_logdet errors on non-amChol", {
  expect_error(am_chol_logdet(list(factor = diag(3))), "amChol")
})

# ── am_chol_diag ─────────────────────────────────────────────────────────────

test_that("am_chol_diag returns diagonal of R", {
  d   <- .lmm_spd
  fac <- am_chol_factor(adgeMatrix(d$S, preferred_backend = "cpu"))
  expect_equal(am_chol_diag(fac), diag(d$ch), tolerance = 1e-10)
})

test_that("am_chol_diag values are positive (SPD input)", {
  d   <- .lmm_spd
  fac <- am_chol_factor(adgeMatrix(d$S, preferred_backend = "cpu"))
  expect_true(all(am_chol_diag(fac) > 0))
})

# ── am_quad_form ──────────────────────────────────────────────────────────────

test_that("am_quad_form vector matches t(v) %*% solve(S) %*% v", {
  d   <- .lmm_spd
  fac <- am_chol_factor(adgeMatrix(d$S, preferred_backend = "cpu"))
  set.seed(7)
  v   <- rnorm(d$p)
  ref <- as.numeric(t(v) %*% solve(d$S) %*% v)
  expect_equal(am_quad_form(fac, v), ref, tolerance = 1e-8)
})

test_that("am_quad_form matrix matches t(V) %*% solve(S) %*% V", {
  d   <- .lmm_spd
  fac <- am_chol_factor(adgeMatrix(d$S, preferred_backend = "cpu"))
  set.seed(8)
  V   <- matrix(rnorm(d$p * 3L), d$p, 3L)
  ref <- t(V) %*% solve(d$S) %*% V
  expect_equal(am_quad_form(fac, V), ref, tolerance = 1e-8)
})

# ── am_eigh ───────────────────────────────────────────────────────────────────

test_that("am_eigh eigenvalues match eigen() on SPD matrix", {
  d   <- .lmm_spd
  ev  <- am_eigh(adgeMatrix(d$S, preferred_backend = "cpu"))
  ref <- eigen(d$S, symmetric = TRUE)
  expect_equal(sort(ev$values, decreasing = TRUE), ref$values, tolerance = 1e-8)
})

test_that("am_eigh eigenvectors are orthonormal", {
  d   <- .lmm_spd
  ev  <- am_eigh(adgeMatrix(d$S, preferred_backend = "cpu"))
  UtU <- t(ev$vectors) %*% ev$vectors
  expect_equal(UtU, diag(d$p), tolerance = 1e-8)
})

test_that("am_eigh satisfies S V = V diag(lambda)", {
  d   <- .lmm_spd
  ev  <- am_eigh(adgeMatrix(d$S, preferred_backend = "cpu"))
  resid <- norm(d$S %*% ev$vectors -
                  ev$vectors %*% diag(ev$values, nrow = length(ev$values)), "F") /
           norm(d$S, "F")
  expect_lt(resid, 1e-8)
})

# ── am_crossprod_weighted ─────────────────────────────────────────────────────

test_that("am_crossprod_weighted matches t(X) %*% diag(w) %*% X", {
  set.seed(21)
  n <- 30L; p <- 5L
  X  <- matrix(rnorm(n * p), n, p)
  w  <- runif(n, 0.5, 2.0)
  aX <- adgeMatrix(X, preferred_backend = "cpu")

  ref <- t(X) %*% diag(w) %*% X
  res <- as.matrix(am_crossprod_weighted(aX, w))
  expect_equal(res, ref, tolerance = 1e-10)
})

test_that("am_crossprod_weighted errors when length(w) != nrow(X)", {
  set.seed(22)
  aX <- adgeMatrix(matrix(rnorm(20), 10, 2), preferred_backend = "cpu")
  expect_error(am_crossprod_weighted(aX, runif(5)), "nrow")
})

# ── am_tcrossprod_weighted ────────────────────────────────────────────────────

test_that("am_tcrossprod_weighted matches tcrossprod(X * sqrt(w))", {
  set.seed(23)
  n <- 20L; p <- 4L
  X  <- matrix(rnorm(n * p), n, p)
  w  <- runif(n, 0.5, 2.0)
  aX <- adgeMatrix(X, preferred_backend = "cpu")

  # Row i scaled by sqrt(w[i]): equivalent to diag(sqrt(w)) %*% X %*% t(X) %*% diag(sqrt(w))
  ref <- tcrossprod(X * sqrt(w))
  res <- as.matrix(am_tcrossprod_weighted(aX, w))
  expect_equal(res, ref, tolerance = 1e-10)
})

# ── am_xty_weighted ───────────────────────────────────────────────────────────

test_that("am_xty_weighted matches t(X) %*% diag(w) %*% y", {
  set.seed(24)
  n <- 25L; p <- 4L; k <- 3L
  X  <- matrix(rnorm(n * p), n, p)
  y  <- matrix(rnorm(n * k), n, k)
  w  <- runif(n, 0.5, 2.0)
  aX <- adgeMatrix(X, preferred_backend = "cpu")

  ref <- t(X) %*% diag(w) %*% y
  res <- as.matrix(am_xty_weighted(aX, w, y))
  expect_equal(res, ref, tolerance = 1e-10)
})

# ── am_wls_fit normal-equations path ─────────────────────────────────────────

test_that("am_wls_fit normal method matches lm.wfit", {
  set.seed(31)
  n <- 40L; p <- 4L
  X  <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  y  <- rnorm(n)
  w  <- runif(n, 0.5, 2.0)

  fit_ref <- lm.wfit(X, y, w = w)
  fit_am  <- am_wls_fit(X, y, weights = w, method = "normal", cache = FALSE)

  ref_coef <- as.numeric(fit_ref$coefficients)
  am_coef  <- as.numeric(coef(fit_am))
  expect_equal(am_coef, ref_coef, tolerance = 1e-8)
})
