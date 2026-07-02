# test-pairwise-argmin.R
#
# Regression/guard tests for pairwise_sqdist_argmin (bead amatrix-p24).
#
# amatrix-p24 alleged a CPU miscompute: that the squared-norm broadcast used
# `rep(c_norms, each = nrow)` where it should use `rep(c_norms, nrow)`. That
# allegation is FALSE, and its proposed fix would INTRODUCE a bug. In R, a
# length-(n*K) vector added to an n×K matrix fills column-major, so
# `rep(c_norms, each = n)` places c_norms[k] across every row of column k —
# which is exactly the correct per-centroid broadcast. The current
# `sweep(..., 2L, c_norms, "+")` is equivalent and also correct. Verified
# against an independent double-loop reference for non-square n≠k≠p inputs.
#
# The bead probe `pairwise_sqdist_argmin(X, C)` returning c(1,1,1) is the
# CORRECT answer: `Ct` is documented (and used by every caller) as p×K with
# *columns* as centroids, and the probe passed C non-transposed. For the
# square degenerate probe the columns of C are (0,5),(0,5), so c(1,1,1) is
# right. Passing t(C) per the contract yields the intuitive c(1,1,2).
#
# These tests pin the correct contract so the bad "fix" cannot be reapplied.

# Independent reference: nearest centroid by explicit double loop.
# `cents` here is k×p (rows are centroids), the natural human-facing layout.
.argmin_ref <- function(X, cents) {
  n <- nrow(X); k <- nrow(cents)
  out <- integer(n)
  for (i in seq_len(n)) {
    best <- Inf; bi <- 1L
    for (kk in seq_len(k)) {
      d <- sum((X[i, ] - cents[kk, ])^2)
      if (d < best) { best <- d; bi <- kk }
    }
    out[i] <- bi
  }
  out
}

test_that("bead probe: contract usage gives intuitive assignment [amatrix-p24]", {
  X <- matrix(c(0, 0, 1, 1, 5, 5), nrow = 3, byrow = TRUE)
  C <- matrix(c(0, 0, 5, 5), nrow = 2, byrow = TRUE)  # rows = centroids

  # Documented contract: Ct is p×K (columns are centroids) -> pass t(C).
  expect_equal(
    as.integer(pairwise_sqdist_argmin(X, t(C))),
    c(1L, 1L, 2L)
  )

  # Untransposed C computes distances to C's COLUMNS (0,5),(0,5); c(1,1,1)
  # is correct for those identical centroids. Pinned so nobody "fixes" it.
  expect_equal(
    as.integer(pairwise_sqdist_argmin(X, C)),
    c(1L, 1L, 1L)
  )
})

test_that("non-square n != k != p matches double-loop reference [amatrix-p24]", {
  set.seed(1)
  n <- 7L; p <- 4L; k <- 3L
  X     <- matrix(rnorm(n * p), n, p)
  cents <- matrix(rnorm(k * p), k, p)   # k×p, rows are centroids
  Ct    <- t(cents)                     # p×K per contract

  expect_equal(
    as.integer(pairwise_sqdist_argmin(X, Ct)),
    as.integer(.argmin_ref(X, cents))
  )
})

test_that("each= vs times= recycling: current code takes the correct branch [amatrix-p24]", {
  # Constructed so that column-major `rep(c_norms, each = n)` (correct) and
  # `rep(c_norms, times = n)` (the bead's proposed wrong fix) disagree.
  set.seed(1)
  n <- 7L; p <- 4L; k <- 3L
  X     <- matrix(rnorm(n * p), n, p)
  cents <- matrix(rnorm(k * p), k, p)
  Ct    <- t(cents)
  xn    <- rowSums(X^2)
  cn    <- colSums(Ct^2)
  cross <- X %*% Ct

  correct <- max.col(-(-2 * cross + xn + rep(cn, each = n)), ties.method = "first")
  wrong   <- max.col(-(-2 * cross + xn + rep(cn, times = n)), ties.method = "first")

  # Guard: this case must actually discriminate the two recyclings.
  expect_false(identical(as.integer(correct), as.integer(wrong)))

  # The exported function must produce the CORRECT (each=/sweep) answer.
  expect_equal(
    as.integer(pairwise_sqdist_argmin(X, Ct)),
    as.integer(correct)
  )
  expect_equal(as.integer(correct), as.integer(.argmin_ref(X, cents)))
})

test_that("precomputed x_norms / c_norms path agrees with reference [amatrix-p24]", {
  set.seed(2)
  n <- 6L; p <- 5L; k <- 4L
  X     <- matrix(rnorm(n * p), n, p)
  cents <- matrix(rnorm(k * p), k, p)
  Ct    <- t(cents)
  ref   <- .argmin_ref(X, cents)

  expect_equal(
    as.integer(pairwise_sqdist_argmin(X, Ct,
                                      x_norms = rowSums(X^2),
                                      c_norms = colSums(Ct^2))),
    as.integer(ref)
  )
})

test_that("mlx backend nearest-centroid matches CPU reference [amatrix-p24]", {
  skip_if_not_installed("amatrix.mlx")
  skip_if_not(
    isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE)),
    "mlx backend not available"
  )
  set.seed(3)
  n <- 12L; p <- 6L; k <- 4L
  X_host <- matrix(rnorm(n * p), n, p)
  cents  <- matrix(rnorm(k * p), k, p)
  Ct     <- t(cents)
  ref    <- .argmin_ref(X_host, cents)

  X_mlx  <- adgeMatrix(X_host, preferred_backend = "mlx")
  got    <- pairwise_sqdist_argmin(X_mlx, Ct)
  expect_equal(as.integer(got), as.integer(ref))
})
