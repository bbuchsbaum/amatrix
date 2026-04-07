# test-segment-primitives.R
#
# Tests for am_segment_sum() and am_segment_mean() — amatrix-ylo
#
# Correctness tolerance:
#   1e-10  CPU backend (float64 throughout)
#   1e-4   GPU backends (float32)
#
# Edge cases tested:
#   - standard random labels
#   - unsorted labels
#   - empty cluster (label absent from labels)
#   - K=1 (all same group)
#   - all-same labels
#   - wide p (p > n)
#   - large K (K = n/2)
#   - non-contiguous resident (result of prior GPU op used as input)

.seg_cpu_tol  <- 1e-10
.seg_gpu_tol  <- 1e-4

# ── Reference implementations ──────────────────────────────────────────────

.ref_segment_sum <- function(X, labels, K) {
  out <- matrix(0, K, ncol(X))
  for (k in seq_len(K)) {
    idx <- which(labels == k)
    if (length(idx) > 0L) out[k, ] <- colSums(X[idx, , drop = FALSE])
  }
  out
}

.ref_segment_mean <- function(X, labels, K) {
  out <- matrix(NA_real_, K, ncol(X))
  for (k in seq_len(K)) {
    idx <- which(labels == k)
    if (length(idx) > 0L) out[k, ] <- colMeans(X[idx, , drop = FALSE])
  }
  out
}

# ── Helpers ────────────────────────────────────────────────────────────────

.check_segment_sum <- function(X, labels, K, tol = .seg_cpu_tol, info = "") {
  ref  <- .ref_segment_sum(X, labels, K)
  got  <- am_segment_sum(X, labels, K)
  if (inherits(got, "adgeMatrix")) got <- as.matrix(got)
  expect_equal(got, ref, tolerance = tol,
               label = paste("am_segment_sum", info))
}

.check_segment_mean <- function(X, labels, K, tol = .seg_cpu_tol, info = "") {
  ref <- .ref_segment_mean(X, labels, K)
  got <- am_segment_mean(X, labels, K)
  if (inherits(got, "adgeMatrix")) got <- as.matrix(got)
  # NA positions must match
  expect_equal(is.na(got), is.na(ref),
               label = paste("am_segment_mean NA mask", info))
  # Non-NA values must match
  expect_equal(got[!is.na(ref)], ref[!is.na(ref)], tolerance = tol,
               label = paste("am_segment_mean values", info))
}

# ── am_segment_sum — correctness ───────────────────────────────────────────

test_that("am_segment_sum: standard random labels, plain matrix input", {
  set.seed(1)
  X      <- matrix(rnorm(200 * 10), 200, 10)
  labels <- sample.int(5L, 200L, replace = TRUE)
  .check_segment_sum(X, labels, 5L)
})

test_that("am_segment_sum: unsorted labels", {
  set.seed(2)
  X      <- matrix(rnorm(100 * 8), 100, 8)
  labels <- c(rep(3L, 40L), rep(1L, 30L), rep(2L, 30L))
  .check_segment_sum(X, labels, 3L)
})

test_that("am_segment_sum: empty cluster (label missing)", {
  set.seed(3)
  X      <- matrix(rnorm(50 * 4), 50, 4)
  labels <- rep(c(1L, 3L), 25L)           # label 2 absent
  ref    <- .ref_segment_sum(X, labels, 3L)
  got    <- am_segment_sum(X, labels, 3L)
  if (inherits(got, "adgeMatrix")) got <- as.matrix(got)
  expect_equal(got[c(1L, 3L), ], ref[c(1L, 3L), ], tolerance = .seg_cpu_tol)
  expect_equal(got[2L, ], rep(0, 4))    # absent → row of zeros in sum
})

test_that("am_segment_sum: K=1 (all rows to one group)", {
  set.seed(4)
  X      <- matrix(rnorm(60 * 5), 60, 5)
  labels <- rep(1L, 60L)
  ref    <- matrix(colSums(X), 1L, 5L)
  got    <- am_segment_sum(X, labels, 1L)
  if (inherits(got, "adgeMatrix")) got <- as.matrix(got)
  expect_equal(got, ref, tolerance = .seg_cpu_tol)
})

test_that("am_segment_sum: wide matrix (p > n)", {
  set.seed(5)
  X      <- matrix(rnorm(20 * 100), 20, 100)
  labels <- sample.int(4L, 20L, replace = TRUE)
  .check_segment_sum(X, labels, 4L)
})

test_that("am_segment_sum: large K (K = n/2)", {
  set.seed(6)
  n      <- 40L; K <- 20L
  X      <- matrix(rnorm(n * 6), n, 6)
  labels <- sample.int(K, n, replace = TRUE)
  .check_segment_sum(X, labels, K)
})

# ── am_segment_mean — correctness ──────────────────────────────────────────

test_that("am_segment_mean: standard random labels, plain matrix input", {
  set.seed(10)
  X      <- matrix(rnorm(300 * 12), 300, 12)
  labels <- sample.int(8L, 300L, replace = TRUE)
  .check_segment_mean(X, labels, 8L)
})

test_that("am_segment_mean: unsorted labels", {
  set.seed(11)
  X      <- matrix(rnorm(80 * 6), 80, 6)
  labels <- c(rep(2L, 30L), rep(4L, 20L), rep(1L, 20L), rep(3L, 10L))
  .check_segment_mean(X, labels, 4L)
})

test_that("am_segment_mean: empty cluster → NA row", {
  set.seed(12)
  X      <- matrix(rnorm(60 * 5), 60, 5)
  labels <- rep(c(1L, 3L), 30L)           # label 2 absent
  ref    <- .ref_segment_mean(X, labels, 3L)
  got    <- am_segment_mean(X, labels, 3L)
  if (inherits(got, "adgeMatrix")) got <- as.matrix(got)
  expect_true(all(is.na(got[2L, ])),      info = "empty cluster row must be all-NA")
  expect_false(any(is.na(got[c(1L,3L),])),info = "non-empty cluster rows must not be NA")
  expect_equal(got[c(1L,3L), ], ref[c(1L,3L), ], tolerance = .seg_cpu_tol)
})

test_that("am_segment_mean: K=1", {
  set.seed(13)
  X      <- matrix(rnorm(80 * 4), 80, 4)
  labels <- rep(1L, 80L)
  ref    <- matrix(colMeans(X), 1L, 4L)
  got    <- am_segment_mean(X, labels, 1L)
  if (inherits(got, "adgeMatrix")) got <- as.matrix(got)
  expect_equal(got, ref, tolerance = .seg_cpu_tol)
})

test_that("am_segment_mean: all-same labels (every row in group 2)", {
  set.seed(14)
  X      <- matrix(rnorm(50 * 3), 50, 3)
  labels <- rep(2L, 50L)
  ref    <- .ref_segment_mean(X, labels, 3L)
  got    <- am_segment_mean(X, labels, 3L)
  if (inherits(got, "adgeMatrix")) got <- as.matrix(got)
  expect_equal(is.na(got), is.na(ref), label = "NA mask")
  expect_equal(got[!is.na(ref)], ref[!is.na(ref)], tolerance = .seg_cpu_tol)
})

test_that("am_segment_mean: wide matrix (p > n)", {
  set.seed(15)
  X      <- matrix(rnorm(30 * 200), 30, 200)
  labels <- sample.int(6L, 30L, replace = TRUE)
  .check_segment_mean(X, labels, 6L)
})

test_that("am_segment_mean: large K", {
  set.seed(16)
  n      <- 100L; K <- 50L
  X      <- matrix(rnorm(n * 8), n, 8)
  labels <- sample.int(K, n, replace = TRUE)
  .check_segment_mean(X, labels, K)
})

# ── adgeMatrix input (CPU backend — falls back to rowsum path) ─────────────

test_that("am_segment_sum: adgeMatrix input with cpu backend falls back correctly", {
  set.seed(20)
  X_mat  <- matrix(rnorm(150 * 7), 150, 7)
  X      <- adgeMatrix(X_mat, preferred_backend = "cpu")
  labels <- sample.int(6L, 150L, replace = TRUE)
  ref    <- .ref_segment_sum(X_mat, labels, 6L)
  got    <- am_segment_sum(X, labels, 6L)
  if (inherits(got, "adgeMatrix")) got <- as.matrix(got)
  expect_equal(got, ref, tolerance = .seg_cpu_tol)
})

test_that("am_segment_mean: adgeMatrix input with cpu backend falls back correctly", {
  set.seed(21)
  X_mat  <- matrix(rnorm(120 * 5), 120, 5)
  X      <- adgeMatrix(X_mat, preferred_backend = "cpu")
  labels <- sample.int(4L, 120L, replace = TRUE)
  .check_segment_mean(X_mat, labels, 4L)
})

# ── GPU backend tests (skip if not available) ───────────────────────────────

.seg_backends <- Filter(function(bk) {
  pkg <- paste0("amatrix.", bk)
  requireNamespace(pkg, quietly = TRUE) &&
    isTRUE(tryCatch(
      get(paste0(pkg, "_is_available"), envir = asNamespace(pkg))(),
      error = function(e) FALSE))
}, c("mlx", "arrayfire"))

for (.bk in .seg_backends) {
  local({
    bk <- .bk
    tol <- .seg_gpu_tol

    test_that(sprintf("am_segment_sum: %s backend returns resident adgeMatrix", bk), {
      set.seed(30)
      n <- 500L; p <- 20L; K <- 10L
      X_mat  <- matrix(rnorm(n * p), n, p)
      X      <- adgeMatrix(X_mat, preferred_backend = bk, precision = "fast")
      labels <- sample.int(K, n, replace = TRUE)
      result <- am_segment_sum(X, labels, K)
      expect_true(inherits(result, "adgeMatrix"),
                  info = sprintf("%s am_segment_sum must return adgeMatrix", bk))
      expect_equal(dim(result), c(K, p))
      ref <- .ref_segment_sum(X_mat, labels, K)
      expect_equal(as.matrix(result), ref, tolerance = tol,
                   label = sprintf("%s am_segment_sum values", bk))
    })

    test_that(sprintf("am_segment_mean: %s backend returns resident adgeMatrix", bk), {
      set.seed(31)
      n <- 500L; p <- 20L; K <- 10L
      X_mat  <- matrix(rnorm(n * p), n, p)
      X      <- adgeMatrix(X_mat, preferred_backend = bk, precision = "fast")
      labels <- sample.int(K, n, replace = TRUE)
      result <- am_segment_mean(X, labels, K)
      expect_true(inherits(result, "adgeMatrix"),
                  info = sprintf("%s am_segment_mean must return adgeMatrix", bk))
      expect_equal(dim(result), c(K, p))
      ref <- .ref_segment_mean(X_mat, labels, K)
      got <- as.matrix(result)
      expect_equal(is.na(got), is.na(ref))
      expect_equal(got[!is.na(ref)], ref[!is.na(ref)], tolerance = tol,
                   label = sprintf("%s am_segment_mean values", bk))
    })

    test_that(sprintf("am_segment_mean: %s empty cluster → NaN/NA propagation", bk), {
      set.seed(32)
      n <- 200L; p <- 8L; K <- 5L
      X_mat  <- matrix(rnorm(n * p), n, p)
      X      <- adgeMatrix(X_mat, preferred_backend = bk, precision = "fast")
      labels <- sample.int(K - 1L, n, replace = TRUE)    # label K absent
      result <- am_segment_mean(X, labels, K)
      got    <- as.matrix(result)
      expect_true(all(is.na(got[K, ])),
                  info = sprintf("%s: absent label K must yield NA row after materialization", bk))
    })

    test_that(sprintf("am_segment_mean: %s non-contiguous resident input", bk), {
      set.seed(33)
      n <- 300L; p <- 15L; K <- 8L
      X_mat  <- matrix(rnorm(n * p), n, p)
      X      <- adgeMatrix(X_mat, preferred_backend = bk, precision = "fast")
      # produce a non-contiguous resident by going through am_sweep
      X2     <- am_sweep(X, 1L, rowSums(X_mat), "-")   # X - rowSums broadcast
      labels <- sample.int(K, n, replace = TRUE)
      X2_mat <- as.matrix(X2)
      result <- am_segment_mean(X2, labels, K)
      ref    <- .ref_segment_mean(X2_mat, labels, K)
      got    <- as.matrix(result)
      expect_equal(got[!is.na(ref)], ref[!is.na(ref)], tolerance = tol,
                   label = sprintf("%s non-contiguous resident segment_mean", bk))
    })
  })
}
