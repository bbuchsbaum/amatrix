# Regression: a backend whose available() probe THROWS must degrade to
# "unavailable" — status, policy, planning, and health must never propagate
# the error.
#
# Bug: with amatrix.mlx loaded from source without its compiled bridge,
# available() -> amatrix_mlx_native_available() -> .Call(...) errored with
# '"amatrix_mlx_native_available_bridge" not resolved from current namespace',
# and amatrix_backend_status() re-threw it (R/backend-registry.R available
# column), killing unrelated callers. Observed 2026-07-09 in
# test-backend-integration.R:122 under devtools::test() on aarch64-apple-darwin
# (R 4.5.1); same class as the red nightly-stress lane
# (mote bd-01KX33B04RZ2243Z21XN541TJF).
#
# Fix: .amatrix_backend_available_safe() wraps every available() consultation
# (status, default fast policy, resident + cold planning); the health probe
# records the error message instead of throwing.

.make_throwing_backend <- function() {
  list(
    capabilities    = function() "matmul",
    features        = function() "dense_f64",
    precision_modes = function() c("strict", "fast"),
    available       = function() stop("native bridge not resolved (simulated)"),
    supports        = function(op, x, y = NULL) identical(op, "matmul"),
    matmul          = function(x, y) x %*% y,
    crossprod       = function(x, y = NULL, ...) base::crossprod(x, y),
    tcrossprod      = function(x, y = NULL, ...) base::tcrossprod(x, y),
    ewise           = function(x, lhs, rhs = NULL, op, ...) {
      if (is.null(rhs)) do.call(op, list(lhs)) else do.call(op, list(lhs, rhs))
    },
    rowSums         = function(x, ...) base::rowSums(x, ...),
    colSums         = function(x, ...) base::colSums(x, ...)
  )
}

.with_throwing_backend <- function(name, code) {
  amatrix_register_backend(name, .make_throwing_backend(), overwrite = TRUE)
  on.exit({
    if (exists(name, envir = amatrix:::.amatrix_state$backends, inherits = FALSE)) {
      rm(list = name, envir = amatrix:::.amatrix_state$backends)
    }
    health <- amatrix:::.amatrix_state$backend_health
    if (!is.null(health) && exists(name, envir = health, inherits = FALSE)) {
      rm(list = name, envir = health)
    }
  }, add = TRUE)
  force(code)
}

test_that("amatrix_backend_status() reports a throwing backend as unavailable, without error", {
  .with_throwing_backend("zz_throwing", {
    status <- NULL
    expect_no_error(status <- amatrix_backend_status())
    row <- status[status$name == "zz_throwing", , drop = FALSE]
    expect_identical(nrow(row), 1L)
    expect_false(row$available)
  })
})

test_that("planning routes to cpu when the preferred backend's probe throws", {
  .with_throwing_backend("zz_throwing", {
    x <- adgeMatrix(matrix(as.numeric(1:4), nrow = 2),
                    preferred_backend = "zz_throwing", policy = "auto",
                    precision = "strict")
    plan <- NULL
    expect_no_error(plan <- amatrix_backend_plan(x, "matmul", y = diag(2)))
    expect_identical(plan$chosen, "cpu")
    entry <- Filter(function(e) identical(e$name, "zz_throwing"), plan$candidates)
    expect_identical(length(entry), 1L)
    expect_false(entry[[1]]$available)
  })
})

test_that("health probe records the probe error instead of throwing", {
  .with_throwing_backend("zz_throwing", {
    result <- NULL
    expect_no_error(result <- amatrix_backend_health_probe("zz_throwing"))
    expect_identical(result$status, "unhealthy")
    expect_true(grepl("available\\(\\) errored", result$reason))
    expect_true(grepl("simulated", result$reason))
  })
})
