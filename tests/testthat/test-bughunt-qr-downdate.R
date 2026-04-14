suppressPackageStartupMessages(library(amatrix))

# amatrix-a6n: qr_downdate.amQR crashes when X is adgeMatrix
# R/qr-downdate.R line 71: X[-row_idx, , drop=FALSE] fails for adgeMatrix
# because the S4 [ method requires an explicit j argument.
test_that("qr_downdate.amQR works when X is an adgeMatrix [amatrix-a6n]", {
  set.seed(42)
  n <- 20L; p <- 4L
  X_am <- as_adgeMatrix(matrix(rnorm(n * p), n, p))
  qr_f <- am_qr(X_am)

  # This currently errors: X[-row_idx, , drop=FALSE] on adgeMatrix
  qr_d <- qr_downdate(qr_f, row_idx = 3L, X = X_am)

  # Result should match downdate with plain matrix X
  X_plain <- as.matrix(X_am)
  qr_d_ref <- qr_downdate(qr_f, row_idx = 3L, X = X_plain)

  R_d   <- qr.R(qr_d)
  R_ref <- qr.R(qr_d_ref)

  # R is unique up to column sign flips
  for (j in seq_len(p)) {
    r1 <- R_d[, j]; r2 <- R_ref[, j]
    if (sign(r1[j]) != sign(r2[j])) r2 <- -r2
    expect_equal(r1, r2, tolerance = 1e-8,
                 label = paste0("R col ", j, " after adgeMatrix downdate"))
  }
})

# amatrix-a6n: lm_loo_cv crashes when X is adgeMatrix (same root cause)
test_that("lm_loo_cv works when X is adgeMatrix [amatrix-a6n]", {
  set.seed(37)
  n <- 15L; p <- 3L
  X      <- matrix(rnorm(n * p), n, p)
  X_am   <- as_adgeMatrix(X)
  y      <- rnorm(n)

  # Reference via plain matrix
  loo_ref <- lm_loo_cv(X, y)

  # This currently errors at qr_downdate(qr_full, i, X = X) where X is adgeMatrix
  loo_am  <- lm_loo_cv(X_am, y)

  expect_equal(loo_am$residuals, loo_ref$residuals, tolerance = 1e-8,
               label = "LOO residuals with adgeMatrix X")
  expect_equal(loo_am$mse, loo_ref$mse, tolerance = 1e-10,
               label = "LOO MSE with adgeMatrix X")
})
