sim_logit_data <- function(n = 500L, p = 6L, seed = 42L) {
  set.seed(seed)
  X <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  colnames(X) <- c("(Intercept)", paste0("x", seq_len(p - 1L)))
  beta <- c(0.5, rnorm(p - 1L, sd = 0.5))
  y <- rbinom(n, 1L, plogis(as.numeric(X %*% beta)))
  list(X = X, y = y)
}

test_that("logit_fit matches stats::glm.fit on a base matrix", {
  d <- sim_logit_data()
  fit <- logit_fit(d$X, d$y)
  ref <- stats::glm.fit(d$X, d$y, family = stats::binomial())

  expect_true(fit$converged)
  expect_equal(unname(fit$coefficients), unname(coef(ref)), tolerance = 1e-6)
  expect_equal(fit$deviance, ref$deviance, tolerance = 1e-6)
})

test_that("logit_fit on a CPU adgeMatrix agrees with the base-matrix path", {
  d <- sim_logit_data(seed = 7L)
  fit_base <- logit_fit(d$X, d$y)
  X_am <- amatrix::adgeMatrix(d$X, preferred_backend = "cpu",
                              precision = "strict")
  fit_am <- logit_fit(X_am, d$y)

  expect_equal(fit_am$coefficients, fit_base$coefficients, tolerance = 1e-10)
  expect_equal(fit_am$deviance, fit_base$deviance, tolerance = 1e-10)
})

test_that("logit_fit on a fast-precision GPU adgeMatrix agrees to float32 tolerance", {
  skip_if_not_installed("amatrix.mlx")
  skip_if_not(
    isTRUE(tryCatch(amatrix.mlx::amatrix_mlx_is_available(),
                    error = function(e) FALSE)),
    "MLX backend not available on this machine"
  )

  d <- sim_logit_data(n = 2000L, p = 20L, seed = 11L)
  fit_base <- logit_fit(d$X, d$y)
  X_gpu <- amatrix::adgeMatrix(d$X, preferred_backend = "mlx",
                               precision = "fast")
  fit_gpu <- logit_fit(X_gpu, d$y)

  expect_true(fit_gpu$converged)
  expect_equal(fit_gpu$coefficients, fit_base$coefficients, tolerance = 1e-3)
})

test_that("logit_fit validates its inputs", {
  d <- sim_logit_data(n = 50L, p = 3L)
  expect_error(logit_fit(d$X, d$y[-1L]), "length")
  expect_error(logit_fit(d$X, d$y + 0.5), "binary")
})

test_that("logit_fit warns on perfect separation instead of failing silently", {
  set.seed(3)
  x <- c(rnorm(50, -3), rnorm(50, 3))
  X <- cbind(1, x)
  y <- as.numeric(x > 0)
  expect_warning(logit_fit(X, y, max_iter = 200L),
                 "probabilities numerically 0 or 1")
})

test_that("logit_fit fails loudly on a rank-deficient design", {
  set.seed(4)
  x <- rnorm(100)
  X <- cbind(1, x, x)  # duplicated column: X'WX exactly singular
  y <- rbinom(100, 1L, plogis(x))
  expect_error(logit_fit(X, y))
})

test_that("logit_fit handles an intercept-only model", {
  set.seed(1)
  y <- rbinom(200L, 1L, 0.3)
  X <- matrix(1, nrow = 200L, ncol = 1L)
  fit <- logit_fit(X, y)
  expect_equal(plogis(fit$coefficients[[1L]]), mean(y), tolerance = 1e-8)
})
