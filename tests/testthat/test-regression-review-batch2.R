# Regressions for the second batch of pre-CRAN review fixes.

test_that("lm_loo_cv matches brute-force leave-one-out", {
  set.seed(4)
  n <- 40
  X <- cbind(1, matrix(rnorm(n * 3), n))
  y <- rnorm(n)
  brute <- vapply(seq_len(n), function(i) {
    y[i] - sum(X[i, ] * qr.coef(qr(X[-i, ]), y[-i]))
  }, numeric(1))
  expect_equal(am_lm_loo_cv(X, y)$residuals, brute, tolerance = 1e-10)
  expect_equal(am_lm_loo_cv(adgeMatrix(X), y)$residuals, brute, tolerance = 1e-10)
  expect_equal(am_lm_loo_cv(X, y)$mse, mean(brute^2), tolerance = 1e-10)
  expect_error(am_lm_loo_cv(cbind(X, X[, 2]), y), class = "amatrix_bad_arg")
  expect_error(am_lm_loo_cv(X, cbind(y, y)), class = "amatrix_bad_arg")
})

test_that("lu_factor stores a factorization that lu_solve reuses", {
  set.seed(6)
  A <- matrix(rnorm(25), 5, dimnames = list(NULL, letters[1:5]))
  b <- rnorm(5)
  B <- matrix(rnorm(10), 5)
  f <- am_lu_factor(A)
  expect_s4_class(f@lu, "denseLU")
  expect_equal(am_lu_solve(f, b), solve(A, b), tolerance = 1e-12)
  expect_equal(am_lu_solve(f, B), solve(A, B), tolerance = 1e-12)
  expect_equal(am_lu_solve(am_lu_factor(adgeMatrix(A)), b), solve(A, b), tolerance = 1e-12)
  expect_error(am_lu_solve(am_lu_factor(matrix(1, 3, 3)), 1:3), "singular")
  expect_error(am_lu_solve(f, 1:3), class = "amatrix_bad_arg")
})

test_that("normal-equation fits report rank without a QR of X", {
  set.seed(7)
  X <- cbind(1, matrix(rnorm(150), 50))
  y <- rnorm(50)
  expect_identical(as.integer(am_lm_fit(X, y, method = "normal")$rank), 4L)
  expect_identical(as.integer(am_lm_fit(X, y, method = "qr")$rank), 4L)
  deficient <- adgeMatrix(crossprod(cbind(X, X[, 2])))
  expect_identical(amatrix:::.amatrix_xtx_rank(deficient), 4L)

  Xs <- Matrix::rsparsematrix(400, 12, 0.3)
  ys <- rnorm(400)
  fit <- am_lm_fit(adgCMatrix(Xs), ys)
  ref <- as.numeric(Matrix::solve(Matrix::crossprod(Xs), Matrix::crossprod(Xs, ys)))
  expect_equal(as.numeric(as.matrix(coef(fit))), ref, tolerance = 1e-8)
  expect_identical(as.integer(fit$rank), 12L)
})

test_that("zero weights do not count towards residual degrees of freedom", {
  set.seed(8)
  X <- cbind(1, matrix(rnorm(150), 50))
  y <- rnorm(50)
  w <- c(rep(0, 10), runif(40, 0.5, 2))
  ref <- stats::lm.wfit(X, y, w)
  for (method in c("normal", "qr")) {
    fit <- am_wls_fit(X, y, w, method = method)
    expect_equal(fit$df.residual, ref$df.residual, info = method)
    expect_equal(as.numeric(as.matrix(coef(fit))), as.numeric(ref$coefficients),
                 tolerance = 1e-8, info = method)
  }
  expect_error(am_wls_fit(X, y, c(Inf, rep(1, 49))), "finite")
})

test_that("ridge_path predict works for several responses and flags off-grid lambda", {
  set.seed(9)
  X <- matrix(rnorm(120), 30)
  Y <- matrix(rnorm(60), 30)
  rp <- am_ridge_path(X, Y, lambdas = c(0.5, 1, 2))
  pred <- predict(rp, X, lambda = 1)
  ref <- X %*% solve(crossprod(X) + diag(4), crossprod(X, Y))
  expect_equal(unname(as.matrix(pred)), unname(ref), tolerance = 1e-8)
  expect_warning(coef(rp, lambda = 3), "not on the fitted grid")
})

test_that("rsvd validates k and am_scatter_mean matches a per-cluster reference", {
  set.seed(10)
  A <- matrix(rnorm(60), 10)
  expect_error(am_rsvd(A, k = 7), "min\\(dim")

  X <- matrix(rnorm(400), 100)
  labels <- sample.int(6L, 100L, replace = TRUE)
  labels[labels == 4L] <- 5L   # leave cluster 4 empty
  ref <- t(vapply(1:6, function(k) {
    idx <- labels == k
    if (any(idx)) colMeans(X[idx, , drop = FALSE]) else rep(NA_real_, 4)
  }, numeric(4)))
  expect_equal(am_scatter_mean(X, labels, 6L), ref)
  expect_equal(am_scatter_mean(adgeMatrix(X), labels, 6L), ref)
  expect_error(am_scatter_mean(X, labels[-1], 6L), "length")
})

test_that("trace_estim(seed=) leaves the caller's RNG stream untouched", {
  A <- crossprod(matrix(sin(1:25), 5)) + 5 * diag(5)
  set.seed(303)
  expected <- runif(2)
  set.seed(303)
  first <- am_trace_estim(A, k = 20L, seed = 1L)
  expect_identical(runif(2), expected)
  expect_identical(am_trace_estim(A, k = 20L, seed = 1L), first)
})

test_that("a user error does not mark a healthy backend unhealthy", {
  counter <- new.env(parent = emptyenv())
  backend <- make_recording_backend(counter)
  # Host-only backend: exercise the cold dispatch path, not resident kernels.
  backend[grepl("resident", names(backend))] <- NULL
  backend$matmul <- function(x, y, ...) {
    if (ncol(x) != NROW(y)) stop("non-conformable arguments")
    as.matrix(x) %*% as.matrix(y)
  }
  with_registered_backend("review_user_error", backend, {
    amatrix:::.amatrix_backend_health_mark("review_user_error", "healthy", NULL)
    x <- adgeMatrix(matrix(1:6 + 0, 2), preferred_backend = "review_user_error")
    expect_error(x %*% matrix(1, 5, 5))
    expect_identical(
      amatrix:::.amatrix_backend_health_get("review_user_error")$status,
      "healthy"
    )
  })
})
