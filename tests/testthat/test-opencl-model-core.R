.opencl_model_spec <- function() {
  specs <- optional_backend_specs()
  specs[[match("opencl", vapply(specs, `[[`, character(1), "backend"))]]
}

.opencl_register_backend <- function(spec) {
  ns <- optional_backend_namespace(spec$package)
  get(spec$register_fun, envir = ns, inherits = FALSE)
}

.expect_opencl_fast_equal <- function(actual, expected, tolerance = 5e-6) {
  expect_equal(unname(as.matrix(actual)), unname(as.matrix(expected)), tolerance = tolerance)
}

test_that("OpenCL normal-equation lm_fit matches lm.fit", {
  spec <- .opencl_model_spec()
  skip_if_backend_package_missing(spec)

  register_backend <- .opencl_register_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    set.seed(20260408L)
    X_host <- matrix(rnorm(220 * 18), nrow = 220, ncol = 18)
    beta_host <- matrix(rnorm(18 * 3), nrow = 18, ncol = 3)
    Y_host <- X_host %*% beta_host + matrix(rnorm(220 * 3, sd = 1e-6), nrow = 220, ncol = 3)

    X_arg <- adgeMatrix(X_host, preferred_backend = "opencl", precision = "fast")
    Y_arg <- adgeMatrix(Y_host, preferred_backend = "opencl", precision = "fast")
    XtX <- am_crossprod(X_arg)
    XtY <- am_crossprod(X_arg, Y_arg)

    fit <- lm_fit(
      X_arg,
      Y_host,
      method = "normal",
      cache = FALSE,
      include_fitted = FALSE,
      include_residuals = FALSE
    )

    ref_coef <- do.call(cbind, lapply(seq_len(ncol(Y_host)), function(j) {
      lm.fit(X_host, Y_host[, j])$coefficients
    }))

    .expect_opencl_fast_equal(fit$coefficients, ref_coef)
    expect_identical(fit$backend, "opencl")
    expect_identical(fit$precision, "fast")
    expect_s4_class(fit$coefficients, "adgeMatrix")
  })
})

test_that("OpenCL normal-equation lm_fit seeds reusable Cholesky cache on XtX", {
  spec <- .opencl_model_spec()
  skip_if_backend_package_missing(spec)

  register_backend <- .opencl_register_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    set.seed(20260413L)
    X_host <- matrix(rnorm(240 * 20), nrow = 240, ncol = 20)
    Y_host <- matrix(rnorm(240 * 2), nrow = 240, ncol = 2)

    fit <- lm_fit(
      adgeMatrix(X_host, preferred_backend = "opencl", precision = "fast"),
      Y_host,
      method = "normal",
      cache = TRUE,
      include_fitted = FALSE,
      include_residuals = FALSE
    )

    expect_true(inherits(.amatrix_cache_get(paste0("chol:", fit$xtx@object_id)), "amChol"))
  })
})

test_that("OpenCL ridge_fit matches penalized normal equations", {
  spec <- .opencl_model_spec()
  skip_if_backend_package_missing(spec)

  register_backend <- .opencl_register_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    set.seed(20260409L)
    X_host <- cbind(1, matrix(rnorm(180 * 10), nrow = 180, ncol = 10))
    beta_host <- matrix(rnorm(ncol(X_host) * 2), nrow = ncol(X_host), ncol = 2)
    Y_host <- X_host %*% beta_host + matrix(rnorm(180 * 2, sd = 1e-6), nrow = 180, ncol = 2)
    lambda <- 0.75

    fit <- ridge_fit(
      adgeMatrix(X_host, preferred_backend = "opencl", precision = "fast"),
      Y_host,
      lambda = lambda,
      penalize_intercept = FALSE,
      include_fitted = FALSE,
      include_residuals = FALSE,
      cache = FALSE
    )

    penalty <- diag(c(0, rep(lambda, ncol(X_host) - 1L)))
    ref_coef <- solve(crossprod(X_host) + penalty, crossprod(X_host, Y_host))

    .expect_opencl_fast_equal(fit$coefficients, ref_coef)
    expect_identical(fit$backend, "opencl")
  })
})

test_that("OpenCL narrow ridge intermediates do not promote into resident solve", {
  spec <- .opencl_model_spec()
  skip_if_backend_package_missing(spec)

  register_backend <- .opencl_register_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    set.seed(20260412L)
    X_host <- matrix(rnorm(512 * 32), nrow = 512, ncol = 32)
    Y_host <- matrix(rnorm(512 * 4), nrow = 512, ncol = 4)

    X_arg <- adgeMatrix(X_host, preferred_backend = "opencl", precision = "fast")
    Y_arg <- adgeMatrix(Y_host, preferred_backend = "opencl", precision = "fast")

    XtX <- am_crossprod(X_arg)
    XtY <- am_crossprod(X_arg, Y_arg)
    penalty <- .amatrix_penalty_matrix(
      X_arg,
      lambda = 0.5,
      penalize_intercept = FALSE,
      has_intercept = FALSE
    )
    penalized <- ewise("+", XtX, penalty)

    expect_false(amatrix_residency_info(XtX)$live[[1]])
    expect_false(amatrix_residency_info(penalty)$live[[1]])
    expect_false(amatrix_residency_info(penalized)$live[[1]])
    fit <- ridge_fit(
      X_arg,
      Y_host,
      lambda = 0.5,
      include_fitted = FALSE,
      include_residuals = FALSE,
      cache = TRUE
    )

    ref_coef <- solve(crossprod(X_host) + diag(0.5, ncol(X_host)), crossprod(X_host, Y_host))
    .expect_opencl_fast_equal(fit$coefficients, ref_coef)
  })
})

test_that("OpenCL wls_fit normal path matches weighted normal equations", {
  spec <- .opencl_model_spec()
  skip_if_backend_package_missing(spec)

  register_backend <- .opencl_register_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    set.seed(20260410L)
    X_host <- cbind(1, matrix(rnorm(160 * 8), nrow = 160, ncol = 8))
    beta_host <- matrix(rnorm(ncol(X_host) * 2), nrow = ncol(X_host), ncol = 2)
    Y_host <- X_host %*% beta_host + matrix(rnorm(160 * 2, sd = 1e-6), nrow = 160, ncol = 2)
    weights <- runif(nrow(X_host), min = 0.2, max = 1.5)

    fit <- wls_fit(
      adgeMatrix(X_host, preferred_backend = "opencl", precision = "fast"),
      Y_host,
      weights = weights,
      method = "normal",
      cache = FALSE,
      include_fitted = FALSE,
      include_residuals = FALSE
    )

    Xw <- X_host * sqrt(weights)
    Yw <- Y_host * sqrt(weights)
    ref_coef <- solve(crossprod(Xw), crossprod(Xw, Yw))

    .expect_opencl_fast_equal(fit$coefficients, ref_coef)
    expect_identical(fit$backend, "opencl")
  })
})

test_that("OpenCL covariance uses the backend path and matches stats::cov", {
  spec <- .opencl_model_spec()
  skip_if_backend_package_missing(spec)

  register_backend <- .opencl_register_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    set.seed(20260411L)
    X_host <- matrix(rnorm(144 * 12), nrow = 144, ncol = 12)
    X_arg <- adgeMatrix(X_host, preferred_backend = "opencl", precision = "fast")

    expect_identical(amatrix_backend_plan(X_arg, "covariance")$chosen, "opencl")

    fit <- covariance(X_arg)
    .expect_opencl_fast_equal(fit, stats::cov(X_host))
  })
})

test_that("OpenCL tall-skinny QR uses resident Q helpers and matches base QR coefficients", {
  spec <- .opencl_model_spec()
  skip_if_backend_package_missing(spec)
  skip_if_not(
    isTRUE(amatrix_backend_status("opencl")$available),
    "opencl backend not available (no usable GPU device on this host)"
  )

  register_backend <- .opencl_register_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    old <- options(
      amatrix.opencl.factor_gpu = TRUE,
      amatrix.opencl.qr_min_n = 1L
    )
    on.exit(options(old), add = TRUE)

    set.seed(20260409L)
    X_host <- matrix(rnorm(2000 * 32), nrow = 2000, ncol = 32)
    Y_host <- matrix(rnorm(2000 * 4), nrow = 2000, ncol = 4)
    X_arg <- adgeMatrix(X_host, preferred_backend = "opencl", precision = "fast")

    qr_fit <- am_qr(X_arg)
    qr_ref <- qr(X_host)
    coef_fit <- qr.coef(qr_fit, Y_host)
    coef_ref <- qr.coef(qr_ref, Y_host)
    fitted_fit <- qr.fitted(qr_fit, Y_host)
    fitted_ref <- qr.fitted(qr_ref, Y_host)
    resid_fit <- qr.resid(qr_fit, Y_host)
    resid_ref <- qr.resid(qr_ref, Y_host)

    expect_identical(.amatrix_qr_kind(qr_fit), "explicit_qr")
    expect_identical(.amatrix_qr_helper_path(qr_fit), "native_resident_backend")
    expect_identical(.amatrix_qr_backend_ops(qr_fit), "opencl")
    expect_true(nzchar(.amatrix_qr_q_key(qr_fit)))
    .expect_opencl_fast_equal(coef_fit, coef_ref, tolerance = 5e-5)
    .expect_opencl_fast_equal(fitted_fit, fitted_ref, tolerance = 5e-5)
    .expect_opencl_fast_equal(resid_fit, resid_ref, tolerance = 5e-5)
  })
})

test_that("OpenCL explicit QR matches base qr.coef on rank-deficient inputs", {
  spec <- .opencl_model_spec()
  skip_if_backend_package_missing(spec)

  register_backend <- .opencl_register_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    old <- options(
      amatrix.opencl.factor_gpu = TRUE,
      amatrix.opencl.qr_min_n = 1L,
      amatrix.opencl.qr_max_p = 16L
    )
    on.exit(options(old), add = TRUE)

    set.seed(20260409L)
    X_base <- matrix(rnorm(240 * 8), nrow = 240, ncol = 8)
    X_host <- X_base
    X_host[, 8] <- X_host[, 1]
    Y_host <- matrix(rnorm(240 * 3), nrow = 240, ncol = 3)
    X_arg <- adgeMatrix(X_host, preferred_backend = "opencl", precision = "fast")

    qr_fit <- am_qr(X_arg)
    qr_ref <- qr(X_host)

    expect_equal(
      as.matrix(qr.coef(qr_fit, Y_host)),
      qr.coef(qr_ref, Y_host),
      tolerance = 5e-5
    )
  })
})

test_that("OpenCL explicit QR honors pivoting when host QR fallback pivots columns", {
  spec <- .opencl_model_spec()
  skip_if_backend_package_missing(spec)

  register_backend <- .opencl_register_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    old <- options(
      amatrix.opencl.factor_gpu = TRUE,
      amatrix.opencl.qr_min_n = 1L,
      amatrix.opencl.qr_max_p = 4L
    )
    on.exit(options(old), add = TRUE)

    set.seed(20260410L)
    X_host <- matrix(rnorm(240 * 8), nrow = 240, ncol = 8)
    X_host[, 1] <- 0
    Y_host <- matrix(rnorm(240 * 3), nrow = 240, ncol = 3)
    X_arg <- adgeMatrix(X_host, preferred_backend = "opencl", precision = "fast")

    qr_fit <- am_qr(X_arg)
    qr_ref <- qr(X_host)

    expect_false(identical(qr_ref$pivot, seq_len(ncol(X_host))))
    expect_true(isTRUE(qr_info(qr_fit)$pivoted))
    expect_equal(
      as.matrix(qr.coef(qr_fit, Y_host)),
      qr.coef(qr_ref, Y_host),
      tolerance = 5e-5
    )
  })
})

test_that("OpenCL planner avoids experimental QR solve by default for dense non-SPD systems", {
  spec <- .opencl_model_spec()
  skip_if_backend_package_missing(spec)

  register_backend <- .opencl_register_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    old <- options(
      amatrix.opencl.factor_gpu = TRUE,
      amatrix.opencl.solve_qr_min_dim = 1L,
      amatrix.opencl.solve_qr_min_rhs = 1L
    )
    on.exit(options(old), add = TRUE)

    set.seed(20260409L)
    A_host <- matrix(rnorm(96 * 96), nrow = 96, ncol = 96) + diag(96) * 0.5
    B_host <- matrix(rnorm(96 * 4), nrow = 96, ncol = 4)
    A_arg <- adgeMatrix(A_host, preferred_backend = "opencl", precision = "fast")
    B_arg <- adgeMatrix(B_host, preferred_backend = "opencl", precision = "fast")

    expect_identical(amatrix_backend_plan(A_arg, "solve", y = B_arg)$chosen, "cpu")
    .expect_opencl_fast_equal(solve(A_arg, B_arg), solve(A_host, B_host), tolerance = 5e-5)
  })
})

test_that("OpenCL QR solve helper matches base solve when experimental route is enabled", {
  spec <- .opencl_model_spec()
  skip_if_backend_package_missing(spec)
  skip_if_not(
    isTRUE(amatrix_backend_status("opencl")$available),
    "opencl backend not available (no usable GPU device on this host)"
  )

  register_backend <- .opencl_register_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    old <- options(
      amatrix.opencl.factor_gpu = TRUE,
      amatrix.opencl.experimental_qr_solve = TRUE,
      amatrix.opencl.solve_qr_min_dim = 1L,
      amatrix.opencl.solve_qr_min_rhs = 1L
    )
    on.exit(options(old), add = TRUE)

    ns <- optional_backend_namespace(spec$package)
    qr_solve_rhs <- get(".amatrix_opencl_qr_solve_rhs", envir = ns, inherits = FALSE)

    set.seed(20260409L)
    A_host <- matrix(rnorm(96 * 96), nrow = 96, ncol = 96) + diag(96) * 0.5
    B_host <- matrix(rnorm(96 * 4), nrow = 96, ncol = 4)
    A_arg <- adgeMatrix(A_host, preferred_backend = "opencl", precision = "fast")
    B_arg <- adgeMatrix(B_host, preferred_backend = "opencl", precision = "fast")

    expect_identical(amatrix_backend_plan(A_arg, "solve", y = B_arg)$chosen, "opencl")
    .expect_opencl_fast_equal(qr_solve_rhs(A_host, B_host), solve(A_host, B_host), tolerance = 5e-5)
    .expect_opencl_fast_equal(solve(A_arg, B_arg), solve(A_host, B_host), tolerance = 5e-5)
  })
})
