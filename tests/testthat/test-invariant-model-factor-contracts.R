# Invariant-driven model/factor contract checks.
#
# Focus:
# - many_lm agrees with column-wise lm_fit across normal and QR paths
# - ridge_fit agrees with the direct penalized normal-equation solution
# - lm_loo_cv matches explicit leave-one-out refits
# - chol factors preserve backend provenance and remain correct after resident release

suppressPackageStartupMessages(library(amatrix))

.model_contract_seed_count <- function(default = 4L) {
  raw <- suppressWarnings(as.integer(Sys.getenv("AMATRIX_STRESS_SEEDS", unset = as.character(default))))
  if (is.na(raw) || raw < 1L) {
    default
  } else {
    raw
  }
}

.model_contract_case <- function(seed) {
  set.seed(seed)
  n <- sample(8:12, 1L)
  p <- sample(3:5, 1L)
  q <- sample(2:4, 1L)

  list(
    X = matrix(rnorm(n * p), nrow = n, ncol = p),
    Y = matrix(rnorm(n * q), nrow = n, ncol = q)
  )
}

.model_contract_lm_ref <- function(X, Y, method) {
  do.call(cbind, lapply(seq_len(ncol(Y)), function(j) {
    as.matrix(coef(lm_fit(X, Y[, j, drop = FALSE], method = method)))
  }))
}

.model_contract_loo_ref <- function(X, y) {
  n <- nrow(X)
  residuals <- vapply(seq_len(n), function(i) {
    fit <- lm_fit(X[-i, , drop = FALSE], y[-i], method = "qr")
    coef_i <- as.numeric(coef(fit))
    y[[i]] - sum(X[i, ] * coef_i)
  }, numeric(1))

  list(residuals = residuals, mse = mean(residuals^2))
}

.invariant_mock_chol_backend <- function(counter) {
  backend <- make_recording_backend(
    counter,
    supported_ops = c("chol"),
    cold_supported_ops = character(),
    resident_supported_ops = c("chol"),
    precision_modes = "fast"
  )

  backend$supports <- function(op, x, y = NULL) {
    inherits(x, "adgeMatrix") &&
      identical(op, "chol") &&
      identical(x@preferred_backend, "mockchol") &&
      identical(x@precision, "fast")
  }

  backend$supports_resident <- function(op, x, y = NULL) {
    inherits(x, "adgeMatrix") &&
      identical(op, "chol") &&
      identical(x@preferred_backend, "mockchol") &&
      identical(x@precision, "fast")
  }

  backend$chol_resident <- function(lhs_key, out_key) {
    value <- chol(backend$resident_materialize(lhs_key))
    backend$resident_store(out_key, value)
    value
  }

  backend
}

test_that("many_lm matches column-wise lm_fit across solver paths", {
  for (seed in seq_len(.model_contract_seed_count(default = 3L))) {
    dat <- .model_contract_case(2026043000L + seed)

    for (method in c("normal", "qr")) {
      info <- sprintf("seed=%d method=%s", seed, method)
      fit_many <- many_lm(dat$X, dat$Y, method = method)
      fit_ref <- .model_contract_lm_ref(dat$X, dat$Y, method = method)

      expect_equal(
        as.matrix(coef(fit_many)),
        fit_ref,
        tolerance = 1e-10,
        info = info
      )
    }
  }
})

test_that("ridge_fit matches direct penalized normal-equation reference", {
  for (seed in seq_len(.model_contract_seed_count(default = 3L))) {
    dat <- .model_contract_case(2026043100L + seed)
    lambda <- c(0.25, 1.5, 3)[[(seed - 1L) %% 3L + 1L]]
    penalty <- diag(lambda, ncol(dat$X))
    ref <- solve(crossprod(dat$X) + penalty, crossprod(dat$X, dat$Y))

    fit <- ridge_fit(dat$X, dat$Y, lambda = lambda)
    expect_equal(
      as.matrix(coef(fit)),
      ref,
      tolerance = 1e-10,
      info = sprintf("seed=%d lambda=%s", seed, lambda)
    )
  }
})

test_that("lm_loo_cv matches explicit leave-one-out refits", {
  for (seed in seq_len(.model_contract_seed_count(default = 3L))) {
    dat <- .model_contract_case(2026043200L + seed)
    y <- dat$Y[, 1L]

    expect_equal(
      lm_loo_cv(dat$X, y),
      .model_contract_loo_ref(dat$X, y),
      tolerance = 1e-10,
      info = sprintf("seed=%d", seed)
    )
  }
})

test_that("chol_factor preserves backend provenance and survives resident release", {
  counter <- new.env(parent = emptyenv())
  backend <- .invariant_mock_chol_backend(counter)

  with_registered_backend("mockchol", backend, {
    set.seed(2026043301L)
    A <- crossprod(matrix(rnorm(60L), nrow = 12L, ncol = 5L)) + diag(5L)
    B <- matrix(rnorm(10L), nrow = 5L, ncol = 2L)
    X <- adgeMatrix(A, preferred_backend = "mockchol", precision = "fast")

    fac <- chol_factor(X)
    expect_identical(fac@backend, "mockchol")
    expect_true(inherits(fac@factor_obj, "adgeMatrix"))
    expect_identical(amatrix:::.amatrix_live_resident_backend(fac@factor_obj), "mockchol")
    expect_equal(chol_solve(fac, B), solve(A, B), tolerance = 1e-10)

    amatrix:::.amatrix_release_resident(fac@factor_obj)
    expect_null(amatrix:::.amatrix_live_resident_backend(fac@factor_obj))
    expect_equal(chol_solve(fac, B), solve(A, B), tolerance = 1e-10)
  })
})
