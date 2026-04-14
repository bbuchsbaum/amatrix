# Regression repro metadata
# Seed: 20260414
# Dimensions: X = 6 x 3, Y = length 6 / 6 x 2
# Backend / precision / dispatch: cpu / strict / host-cold QR and normal paths
# R version / platform: R 4.5.1, aarch64-apple-darwin20 (macOS Sonoma 14.3)
# Issue: amatrix-aqm

suppressPackageStartupMessages(library(amatrix))

.model_nonfinite_fixture <- function() {
  set.seed(20260414L)
  list(
    X = matrix(rnorm(18L), nrow = 6L, ncol = 3L),
    y = rnorm(6L),
    Y = matrix(rnorm(12L), nrow = 6L, ncol = 2L)
  )
}

test_that("lm wrappers reject non-finite X with a classed bad-arg error [amatrix-aqm]", {
  dat <- .model_nonfinite_fixture()
  dat$X[2L, 3L] <- NaN

  expect_error(lm_fit(dat$X, dat$y, method = "qr"), class = "amatrix_bad_arg")
  expect_error(lm_fit(dat$X, dat$y, method = "normal"), class = "amatrix_bad_arg")
  expect_error(many_lm(dat$X, dat$Y, method = "qr"), class = "amatrix_bad_arg")
  expect_error(many_lm(dat$X, dat$Y, method = "normal"), class = "amatrix_bad_arg")
  expect_error(lm_loo_cv(dat$X, dat$y), class = "amatrix_bad_arg")
})

test_that("lm wrappers reject non-finite Y with a classed bad-arg error [amatrix-aqm]", {
  dat <- .model_nonfinite_fixture()
  dat$y[4L] <- Inf
  dat$Y[1L, 2L] <- NA_real_

  expect_error(lm_fit(dat$X, dat$y, method = "qr"), class = "amatrix_bad_arg")
  expect_error(lm_fit(dat$X, dat$y, method = "normal"), class = "amatrix_bad_arg")
  expect_error(many_lm(dat$X, dat$Y, method = "qr"), class = "amatrix_bad_arg")
  expect_error(many_lm(dat$X, dat$Y, method = "normal"), class = "amatrix_bad_arg")
  expect_error(lm_loo_cv(dat$X, dat$y), class = "amatrix_bad_arg")
})
