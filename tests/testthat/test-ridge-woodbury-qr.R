suppressPackageStartupMessages(library(amatrix))

# ---- ridge_path -------------------------------------------------------------

test_that("ridge_path coefficients match direct formula", {
  set.seed(42)
  n <- 50L; p <- 8L
  X  <- matrix(rnorm(n * p), n, p)
  y  <- rnorm(n)
  lambdas <- c(0.1, 1.0, 10.0)

  rp <- ridge_path(X, y, lambdas = lambdas)

  expect_s3_class(rp, "ridge_path")
  expect_equal(rp$lambdas, lambdas)
  expect_equal(dim(rp$coef), c(p, 1L, length(lambdas)))

  # direct: beta = (X^T X + lambda I)^{-1} X^T y
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  for (i in seq_along(lambdas)) {
    ref <- solve(XtX + lambdas[[i]] * diag(p), Xty)
    expect_equal(rp$coef[, 1L, i], as.numeric(ref), tolerance = 1e-8,
                 label = paste0("ridge coef lambda=", lambdas[[i]]))
  }
})

test_that("ridge_path coef() extracts single lambda slice", {
  set.seed(1)
  X  <- matrix(rnorm(30 * 5), 30, 5)
  y  <- rnorm(30)
  rp <- ridge_path(X, y, lambdas = c(0.5, 2.0, 5.0))

  coef_all  <- coef(rp)
  coef_near <- coef(rp, lambda = 2.0)

  expect_equal(dim(coef_all),  c(5L, 1L, 3L))
  expect_equal(dim(coef_near), c(5L, 1L, 1L))
  expect_equal(coef_near[,, 1L], coef_all[,, 2L])
})

test_that("ridge_path predict matches manual X %*% coef", {
  set.seed(7)
  X    <- matrix(rnorm(40 * 4), 40, 4)
  y    <- rnorm(40)
  Xnew <- matrix(rnorm(10 * 4), 10, 4)
  rp   <- ridge_path(X, y, lambdas = c(1.0, 5.0))

  pred_single <- predict(rp, Xnew, lambda = 1.0)
  ref_coef    <- coef(rp, lambda = 1.0)[,, 1L]
  expect_equal(as.numeric(pred_single), as.numeric(Xnew %*% ref_coef),
               tolerance = 1e-10)
})

test_that("ridge_path multi-response Y works", {
  set.seed(3)
  X  <- matrix(rnorm(30 * 4), 30, 4)
  Y  <- matrix(rnorm(30 * 2), 30, 2)
  rp <- ridge_path(X, Y, lambdas = c(0.1, 1.0))

  expect_equal(dim(rp$coef), c(4L, 2L, 2L))
})

# ---- woodbury_solve ---------------------------------------------------------

test_that("woodbury_solve matches direct solve for rank-1 update", {
  set.seed(11)
  n <- 20L; k <- 1L
  A <- crossprod(matrix(rnorm(n * n), n)) + diag(n)   # SPD
  U <- matrix(rnorm(n * k), n, k)
  b <- matrix(rnorm(n * 2), n, 2)

  fac <- chol_factor(as_adgeMatrix(A))
  ws  <- woodbury_solve(fac, U, b)   # (A + U U^T)^{-1} b

  direct <- solve(A + U %*% t(U), b)
  expect_equal(ws, direct, tolerance = 1e-8)
})

test_that("woodbury_solve handles rank-k update and explicit V", {
  set.seed(13)
  n <- 15L; k <- 3L
  A <- crossprod(matrix(rnorm(n * n), n)) + 2 * diag(n)
  U <- matrix(rnorm(n * k), n, k)
  V <- matrix(rnorm(k * n), k, n)
  C_inv <- diag(k) * 0.5
  b <- rnorm(n)

  C <- solve(C_inv)
  fac <- chol_factor(as_adgeMatrix(A))
  ws  <- woodbury_solve(fac, U, matrix(b), V = V, C_inv = C_inv)
  direct <- solve(A + U %*% C %*% V, b)
  expect_equal(as.numeric(ws), direct, tolerance = 1e-8)
})

test_that("woodbury_solve auto-factors plain matrix A", {
  set.seed(17)
  n <- 10L
  A <- crossprod(matrix(rnorm(n * n), n)) + diag(n)
  U <- matrix(rnorm(n), n, 1L)
  b <- rnorm(n)

  ws     <- woodbury_solve(A, U, matrix(b))
  direct <- solve(A + U %*% t(U), b)
  expect_equal(as.numeric(ws), direct, tolerance = 1e-8)
})

# ---- woodbury_logdet --------------------------------------------------------

test_that("woodbury_logdet matches log(det(A + UU^T))", {
  set.seed(21)
  n <- 12L; k <- 2L
  A <- crossprod(matrix(rnorm(n * n), n)) + diag(n)
  U <- matrix(rnorm(n * k), n, k)

  fac <- chol_factor(as_adgeMatrix(A))
  wl  <- woodbury_logdet(fac, U)
  ref <- as.numeric(determinant(A + U %*% t(U), logarithm = TRUE)$modulus)
  expect_equal(wl, ref, tolerance = 1e-8)
})

test_that("woodbury_logdet works with explicit C_inv", {
  set.seed(23)
  n <- 10L; k <- 2L
  A     <- crossprod(matrix(rnorm(n * n), n)) + diag(n)
  U     <- matrix(rnorm(n * k), n, k)
  C_inv <- diag(k) * 2.0
  C     <- solve(C_inv)

  fac <- chol_factor(as_adgeMatrix(A))
  wl  <- woodbury_logdet(fac, U, C_inv = C_inv)
  ref <- as.numeric(determinant(A + U %*% C %*% t(U), logarithm = TRUE)$modulus)
  expect_equal(wl, ref, tolerance = 1e-8)
})

# ---- qr_downdate / lm_loo_cv ------------------------------------------------

test_that("qr_downdate.amQR gives same R as am_qr on subset", {
  set.seed(31)
  n <- 20L; p <- 4L
  X    <- matrix(rnorm(n * p), n, p)
  Xam  <- as_adgeMatrix(X)
  i    <- 7L
  qr_f <- am_qr(Xam)

  qr_d <- qr_downdate(qr_f, i, X = X)

  # Compare R matrices (upper triangular part)
  R_down <- qr.R(qr_d)
  R_ref  <- qr.R(am_qr(as_adgeMatrix(X[-i, , drop = FALSE])))

  # R is unique up to sign of columns
  for (j in seq_len(p)) {
    r1 <- R_down[, j]; r2 <- R_ref[, j]
    if (sign(r1[j]) != sign(r2[j])) r2 <- -r2
    expect_equal(r1, r2, tolerance = 1e-8,
                 label = paste0("R col ", j, " after downdate row ", i))
  }
})

test_that("lm_loo_cv residuals match manual LOO loop", {
  set.seed(37)
  n <- 15L; p <- 3L
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  loo <- lm_loo_cv(X, y)

  ref_resid <- numeric(n)
  for (i in seq_len(n)) {
    fit_i <- lm.fit(X[-i, , drop = FALSE], y[-i])
    ref_resid[[i]] <- y[[i]] - sum(X[i, ] * fit_i$coefficients)
  }

  expect_equal(loo$residuals, ref_resid, tolerance = 1e-8)
  expect_equal(loo$mse, mean(ref_resid^2), tolerance = 1e-10)
})

test_that("qr_downdate.default errors informatively", {
  qr_base <- qr(matrix(rnorm(12), 4, 3))
  expect_error(qr_downdate(qr_base, 1L), "amQR factor")
})
