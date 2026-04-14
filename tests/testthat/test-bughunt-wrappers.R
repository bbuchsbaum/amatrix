library(testthat)

# Load package — tests run via devtools::test() which already loads the package,
# but be explicit for standalone runs.
# devtools::load_all(".")

# ── amatrix-3bx: dist_matrix no column dimension check ─────────────────────
# Bug: dist_matrix(X, Y) with ncol(X) != ncol(Y) should error; instead it
# silently passes mismatched matrices to tcrossprod, crashing or producing
# garbage dimensions.
test_that("amatrix-3bx: dist_matrix errors when ncol(X) != ncol(Y)", {
  X <- matrix(rnorm(12), nrow = 3, ncol = 4)
  Y <- matrix(rnorm(15), nrow = 3, ncol = 5)
  expect_error(dist_matrix(X, Y), regexp = NULL)
})

# ── amatrix-w20: vector input silently coerced to 1-col matrix ─────────────
# Bug: dist_matrix(v) where v is a plain numeric vector converts it to an
# ncol=1 matrix instead of raising an error. The user almost certainly meant
# a row-vector or made a mistake.
test_that("amatrix-w20: dist_matrix errors on plain vector input", {
  v <- rnorm(5)
  expect_error(dist_matrix(v), regexp = NULL)
})

test_that("amatrix-w20: kernel_matrix errors on plain vector input", {
  v <- rnorm(5)
  expect_error(kernel_matrix(v), regexp = NULL)
})

# ── amatrix-5ul: diagonal fixup skipped when Y is explicit ─────────────────
# Bug: dist_matrix(X) zeroes the diagonal, but dist_matrix(X, X) does not.
# For the Euclidean metric, self-distances must be exactly 0. Floating-point
# arithmetic can produce tiny nonzero values (~1e-15) that are not zeroed.
test_that("amatrix-5ul: dist_matrix(X,X) diagonal equals dist_matrix(X) diagonal", {
  set.seed(42)
  X <- matrix(rnorm(30), nrow = 6, ncol = 5)
  D1 <- dist_matrix(X)
  D2 <- dist_matrix(X, X)
  # dist_matrix(X) zeroes the diagonal; dist_matrix(X,X) should too
  expect_equal(diag(D1), rep(0.0, 6), tolerance = 0)
  expect_equal(diag(D2), rep(0.0, 6), tolerance = 0)
})

test_that("amatrix-5ul: dist_matrix(X,X) and dist_matrix(X) produce same matrix", {
  set.seed(7)
  X <- matrix(rnorm(20), nrow = 4, ncol = 5)
  D1 <- dist_matrix(X)
  D2 <- dist_matrix(X, X)
  expect_equal(D1, D2, tolerance = 1e-10)
})

# ── amatrix-inz: sigma=0 causes NaN in rbf / laplacian ────────────────────
# Bug: kernel_matrix with sigma=0 computes exp(-D_sq/(2*0^2)) = exp(-Inf) or
# exp(NaN). No guard — result is NaN. Should error with a clear message.
test_that("amatrix-inz: kernel_matrix rbf sigma=0 errors or returns finite values", {
  X <- matrix(rnorm(12), nrow = 3)
  expect_error(
    kernel_matrix(X, kernel = "rbf", sigma = 0),
    "sigma",
    label = "kernel_matrix rbf sigma=0 should fail fast"
  )
})

test_that("amatrix-inz: kernel_matrix laplacian sigma=0 errors or returns finite values", {
  X <- matrix(rnorm(12), nrow = 3)
  expect_error(
    kernel_matrix(X, kernel = "laplacian", sigma = 0),
    "sigma",
    label = "kernel_matrix laplacian sigma=0 should fail fast"
  )
})

# ── amatrix-n5x: negative sigma inverts laplacian kernel ──────────────────
# Bug: sigma=-1 in laplacian computes exp(-||x-y|| / -1) = exp(+||x-y||),
# producing values >> 1 and a non-PSD kernel. No validation of sigma sign.
test_that("amatrix-n5x: kernel_matrix laplacian with negative sigma errors or warns", {
  X <- matrix(rnorm(12), nrow = 3)
  expect_error(
    kernel_matrix(X, kernel = "laplacian", sigma = -1),
    "sigma",
    label = "laplacian kernel must reject negative sigma"
  )
})

test_that("amatrix-n5x: kernel_matrix rbf with negative sigma errors or warns", {
  X <- matrix(rnorm(12), nrow = 3)
  expect_error(
    kernel_matrix(X, kernel = "rbf", sigma = -1),
    "sigma",
    label = "rbf kernel must reject negative sigma"
  )
})

# ── amatrix-t7r: zero_diag=TRUE silently ignored in non-resident path ──────
# Bug: kernel_matrix(..., zero_diag=TRUE) is only honoured when a
# kernel_resident backend is available. The non-resident path (lines 2710-2711)
# calls .am_kernel_gpu() and drops zero_diag entirely. The diagonal of the
# returned matrix remains the kernel self-similarities (== 1 for normalised
# kernels) instead of 0.
test_that("amatrix-t7r: kernel_matrix zero_diag=TRUE zeroes diagonal (non-resident path)", {
  X <- matrix(rnorm(20), nrow = 4)
  # Force cpu non-resident path
  K <- kernel_matrix(X, kernel = "rbf", sigma = 1.0,
                     preferred_backend = NULL, zero_diag = TRUE)
  expect_equal(diag(K), rep(0.0, 4),
    tolerance = 1e-10,
    label = "diagonal should be zero when zero_diag=TRUE")
})

# ── amatrix-8li: NA in X propagates silently ──────────────────────────────
# Bug: dist_matrix / kernel_matrix accept X with NA values without error or
# warning. The NA propagates through rowSums and outer(), silently producing
# all-NA rows. Users get garbage results with no indication of the problem.
test_that("amatrix-8li: dist_matrix warns or errors on NA in X", {
  X <- matrix(rnorm(12), nrow = 3)
  X[2, 1] <- NA
  expect_error(
    dist_matrix(X),
    regexp = NULL,
    label = "dist_matrix should error/warn on NA input"
  )
})

test_that("amatrix-8li: kernel_matrix warns or errors on NA in X", {
  X <- matrix(rnorm(12), nrow = 3)
  X[1, 2] <- NA
  expect_error(
    kernel_matrix(X, kernel = "rbf", sigma = 1.0),
    regexp = NULL,
    label = "kernel_matrix should error/warn on NA input"
  )
})

# ── Regression: euclidean self-distance is exactly 0 ──────────────────────
# Basic sanity: every row's distance to itself is 0.
test_that("dist_matrix self-distance is zero on diagonal", {
  X <- matrix(rnorm(25), nrow = 5)
  D <- dist_matrix(X)
  expect_equal(diag(D), rep(0.0, 5), tolerance = 0)
})

# ── amatrix-f9k: rbf self-kernel diagonal not exactly 1 (float32 drift) ───
# Bug: kernel_matrix(X, kernel="rbf") diagonal should be exp(0)=1, but float32
# accumulation in rowSums(x^2) - 2*G can produce tiny nonzero D_sq for the
# diagonal, giving exp(-eps) < 1. The diag(...) <- 0 guard is only applied
# when is.null(y_host); without Y the diagonal should be forced to 1.
test_that("amatrix-f9k: kernel_matrix rbf diagonal is exactly 1 for self-kernel", {
  X <- matrix(rnorm(15), nrow = 3)
  K <- kernel_matrix(X, kernel = "rbf", sigma = 1.0)
  expect_equal(diag(K), rep(1.0, 3), tolerance = 1e-7)
})
