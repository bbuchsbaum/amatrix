# Track 6 — determinism checks.
#
# Same input + seed + backend should produce bit-identical output across
# repetitions. Catches nondeterministic GPU reduction orders, thread-race
# bugs, and backend state leaks that would break downstream reproducibility.
#
# Scope: pure functional ops routed through amatrix. Iterative algorithms
# that depend on random sketches (rsvd, block_lanczos) must be run with a
# fixed seed under the test.
#
# Under test-residency-tripwire.R we already check that host-only ops do not
# fire residency transfers. Here we check that identical calls produce
# identical bits.

.td_repeat <- 10L  # trade off coverage vs test wall-time

test_that("dense matmul is bit-identical across repetitions", {
  set.seed(2026041400L)
  x_host <- matrix(rnorm(200L), nrow = 10L, ncol = 20L)
  y_host <- matrix(rnorm(200L), nrow = 20L, ncol = 10L)

  x <- adgeMatrix(x_host)
  y <- adgeMatrix(y_host)

  ref <- as.matrix(matmul(x, y))
  for (i in seq_len(.td_repeat)) {
    this <- as.matrix(matmul(x, y))
    expect_identical(this, ref, info = paste("matmul iteration", i))
  }
})

test_that("crossprod / tcrossprod are bit-identical across repetitions", {
  set.seed(2026041401L)
  x_host <- matrix(rnorm(200L), nrow = 10L, ncol = 20L)
  x <- adgeMatrix(x_host)

  ref_cp <- as.matrix(crossprod(x))
  ref_tcp <- as.matrix(tcrossprod(x))
  for (i in seq_len(.td_repeat)) {
    expect_identical(as.matrix(crossprod(x)), ref_cp,
                     info = paste("crossprod iteration", i))
    expect_identical(as.matrix(tcrossprod(x)), ref_tcp,
                     info = paste("tcrossprod iteration", i))
  }
})

test_that("reductions are bit-identical across repetitions", {
  set.seed(2026041402L)
  x_host <- matrix(rnorm(300L), nrow = 15L, ncol = 20L)
  x <- adgeMatrix(x_host)

  ref_row <- rowSums(x)
  ref_col <- colSums(x)
  ref_rm  <- rowmeans(x)
  ref_cm  <- colmeans(x)

  for (i in seq_len(.td_repeat)) {
    expect_identical(rowSums(x), ref_row, info = paste("rowSums", i))
    expect_identical(colSums(x), ref_col, info = paste("colSums", i))
    expect_identical(rowmeans(x), ref_rm, info = paste("rowmeans", i))
    expect_identical(colmeans(x), ref_cm, info = paste("colmeans", i))
  }
})

test_that("solve and chol_factor are bit-identical across repetitions", {
  set.seed(2026041403L)
  n <- 8L
  a <- crossprod(matrix(rnorm(n * n), n, n)) + diag(n) * 2
  x <- as_adgeMatrix(a)

  rhs <- rnorm(n)
  ref_solve <- solve(x, rhs)
  ref_chol <- as.matrix(chol_factor(x))

  for (i in seq_len(.td_repeat)) {
    expect_identical(solve(x, rhs), ref_solve,
                     info = paste("solve iteration", i))
    expect_identical(as.matrix(chol_factor(x)), ref_chol,
                     info = paste("chol_factor iteration", i))
  }
})

test_that("seeded iterative algorithms are deterministic at fixed seed", {
  # rsvd and block_lanczos use random sketches; determinism holds only under
  # a fixed seed. This test documents that contract and pins it as a
  # regression check.
  set.seed(2026041404L)
  x_host <- matrix(rnorm(300L), nrow = 30L, ncol = 10L)

  run_rsvd <- function() {
    set.seed(42L)
    rsvd(x_host, k = 3L)
  }
  ref <- run_rsvd()
  for (i in seq_len(.td_repeat)) {
    this <- run_rsvd()
    expect_equal(this$d, ref$d, tolerance = 0,
                 info = paste("rsvd singular values iteration", i))
    expect_equal(this$u, ref$u, tolerance = 0,
                 info = paste("rsvd left vectors iteration", i))
    expect_equal(this$v, ref$v, tolerance = 0,
                 info = paste("rsvd right vectors iteration", i))
  }
})
