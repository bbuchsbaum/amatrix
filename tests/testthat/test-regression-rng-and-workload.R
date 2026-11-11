# Regression tests for pre-CRAN review findings:
#   * loading the package and probing a backend must not touch the user's RNG
#   * calibration workloads must not overflow integer arithmetic
#   * a malformed on-disk calibration cache must not break dispatch

test_that(".onLoad does not perturb the user's RNG stream", {
  state <- amatrix:::.amatrix_state
  old_session_id <- state$session_id
  on.exit(state$session_id <- old_session_id, add = TRUE)

  set.seed(101)
  expected <- runif(3)

  set.seed(101)
  amatrix:::.onLoad("", "amatrix")
  expect_identical(runif(3), expected)

  # Session ids must still be unique across repeated initialisation.
  first <- state$session_id
  amatrix:::.onLoad("", "amatrix")
  expect_false(identical(state$session_id, first))
})

test_that("amatrix_backend_health_probe does not reset the user's RNG stream", {
  state <- amatrix:::.amatrix_state
  old_health <- state$backend_health
  on.exit(state$backend_health <- old_health, add = TRUE)

  set.seed(202)
  expected <- runif(3)

  set.seed(202)
  amatrix_backend_health_probe("cpu")
  expect_identical(runif(3), expected)
  expect_identical(amatrix:::.amatrix_backend_health_get("cpu")$status, "healthy")
})

test_that("dispatch workload does not overflow for large products", {
  x <- adgeMatrix(matrix(0, 2L, 2L))
  y <- adgeMatrix(matrix(0, 2L, 2L))
  # Pretend-large shapes without allocating: 2000 x 1200 %*% 1200 x 1000.
  x@Dim <- c(2000L, 1200L)
  y@Dim <- c(1200L, 1000L)
  w <- expect_no_warning(amatrix:::.amatrix_dispatch_workload(x, "matmul", y = y))
  expect_identical(w, 2000 * 1200 * 1000)

  expect_identical(
    expect_no_warning(amatrix:::.amatrix_benchmark_workload("gemm", 4096L, 4096L)),
    4096^3
  )
})

test_that("malformed calibration cache file is ignored rather than breaking dispatch", {
  state <- amatrix:::.amatrix_state
  old_calibration <- state$calibration
  cache_dir <- withr::local_tempdir()
  withr::local_envvar(R_USER_CACHE_DIR = cache_dir)
  on.exit(state$calibration <- old_calibration, add = TRUE)

  path <- amatrix:::.amatrix_calibration_path()
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(1:3, path)

  state$calibration <- NULL
  expect_no_error(amatrix:::.amatrix_load_calibration())
  expect_null(state$calibration)
})

test_that("weighted QR cache is not reused for different weights with equal moments", {
  set.seed(1)
  X <- matrix(rnorm(200), 50)
  y <- rnorm(50)
  w <- runif(50, 0.5, 2)
  XA <- adgeMatrix(X)

  f1 <- wls_fit(XA, y, w, method = "qr")
  # rev(w) has the same length, sum and sum of squares as w.
  f2 <- wls_fit(XA, y, rev(w), method = "qr")
  expect_equal(
    as.numeric(as.matrix(coef(f2))),
    as.numeric(stats::lm.wfit(X, y, rev(w))$coefficients),
    tolerance = 1e-10
  )
  # Identical weights still hit the cache.
  f3 <- wls_fit(XA, y, rev(w), method = "qr")
  expect_true(isTRUE(f3$cache_reused))
  expect_equal(as.matrix(coef(f3)), as.matrix(coef(f2)))
})

test_that("chol_solve rejects a right-hand side with the wrong number of rows", {
  S <- crossprod(matrix(sin(1:36), 6)) + diag(6)
  f <- chol_factor(adgeMatrix(S))
  expect_error(chol_solve(f, rep(1, 7)), class = "amatrix_bad_arg")
  expect_error(chol_solve(f, rep(1, 3)), class = "amatrix_bad_arg")
  expect_equal(as.numeric(chol_solve(f, rep(1, 6))), as.numeric(solve(S, rep(1, 6))))
})
