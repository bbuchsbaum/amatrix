# Track 3 residency tripwire tests.
#
# The residency tripwire instruments every real GPU-to-host transfer in
# `R/residency.R::amatrix_materialize_dense()`. Tests enable it, exercise
# dispatch-sensitive ops, and assert that the counter is 0 afterwards.
#
# A non-zero counter means S4 dispatch fell through to a base R / Matrix
# generic and silently coerced a GPU-resident tensor back to host — which
# is both a correctness-visible performance bug and a stop-ship condition
# under planning_docs/quality-tracking.md §7 rule 6.
#
# On a CPU-only runner (no GPU backend loaded), the tripwire is trivially
# satisfied because no resident_materialize() call happens at all. The
# infrastructure tests below still verify that the counter / events / reset
# machinery is wired up correctly.

test_that("tripwire infrastructure: default disabled, counter tracks events", {
  withr::local_options(amatrix.residency.tripwire = NULL)
  withr::local_envvar(AMATRIX_RESIDENCY_TRIPWIRE = NA)
  amatrix:::.amatrix_tripwire_reset()

  expect_false(amatrix:::.amatrix_tripwire_enabled())
  expect_identical(amatrix:::.amatrix_tripwire_count(), 0L)
  expect_length(amatrix:::.amatrix_tripwire_events(), 0L)

  # With tripwire disabled, record() is a no-op.
  amatrix:::.amatrix_tripwire_record("test", "fake", "key")
  expect_identical(amatrix:::.amatrix_tripwire_count(), 0L)

  # Enable via option and verify record() fires.
  withr::local_options(amatrix.residency.tripwire = TRUE)
  expect_true(amatrix:::.amatrix_tripwire_enabled())
  amatrix:::.amatrix_tripwire_record("op_a", "cpu", "k1", c(2L, 3L))
  amatrix:::.amatrix_tripwire_record("op_b", "cpu", "k2", c(5L, 5L))
  expect_identical(amatrix:::.amatrix_tripwire_count(), 2L)

  events <- amatrix:::.amatrix_tripwire_events()
  expect_length(events, 2L)
  expect_identical(events[[1L]]$op, "op_a")
  expect_identical(events[[1L]]$key, "k1")
  expect_identical(events[[2L]]$op, "op_b")
  expect_identical(events[[2L]]$dim, c(5L, 5L))

  amatrix:::.amatrix_tripwire_reset()
  expect_identical(amatrix:::.amatrix_tripwire_count(), 0L)
  expect_length(amatrix:::.amatrix_tripwire_events(), 0L)
})

test_that("tripwire respects AMATRIX_RESIDENCY_TRIPWIRE env var", {
  withr::local_options(amatrix.residency.tripwire = NULL)
  withr::local_envvar(AMATRIX_RESIDENCY_TRIPWIRE = "1")
  expect_true(amatrix:::.amatrix_tripwire_enabled())

  withr::local_envvar(AMATRIX_RESIDENCY_TRIPWIRE = "0")
  expect_false(amatrix:::.amatrix_tripwire_enabled())
})

test_that("host-only ops never fire the tripwire", {
  withr::local_options(amatrix.residency.tripwire = TRUE)
  amatrix:::.amatrix_tripwire_reset()

  set.seed(2026041300L)
  x_host <- matrix(rnorm(20L), nrow = 5L, ncol = 4L)
  y_host <- matrix(rnorm(12L), nrow = 4L, ncol = 3L)
  z_host <- matrix(rnorm(15L), nrow = 5L, ncol = 3L)

  x <- adgeMatrix(x_host)
  y <- adgeMatrix(y_host)
  z <- adgeMatrix(z_host)

  # Exercise the core dispatch surface that Track 2 hardened.
  .quiet <- as.matrix(x %*% y)
  .quiet <- as.matrix(crossprod(x))
  .quiet <- as.matrix(tcrossprod(x))
  .quiet <- rowSums(x)
  .quiet <- colSums(x)
  .quiet <- as.matrix(crossprod(x, z))   # t(x) %*% z — aTransposeView path
  .quiet <- as.matrix(tcrossprod(x, adgeMatrix(t(y_host))))

  # Plain matrix LHS / numeric LHS — exercises dispatch-hardening.R paths.
  .quiet <- as.matrix(x_host %*% y)
  .quiet <- as.matrix(crossprod(x_host, z))
  .quiet <- as.matrix(tcrossprod(x_host, adgeMatrix(t(y_host))))

  expect_identical(
    amatrix:::.amatrix_tripwire_count(),
    0L,
    info = paste(
      "Host-only ops should never trigger GPU materialization.",
      "A non-zero count here means the tripwire fired on a code path",
      "that shouldn't have had any GPU residency."
    )
  )
})

test_that("tripwire counts repeated materializations of a resident-like object", {
  # A synthetic test that simulates what the tripwire should catch. We cannot
  # easily fake a resident backend in-process, so we verify the tripwire by
  # manually calling .amatrix_tripwire_record() with realistic event data.
  # This documents the counter's semantics without depending on a GPU.
  withr::local_options(amatrix.residency.tripwire = TRUE)
  amatrix:::.amatrix_tripwire_reset()

  for (i in seq_len(5L)) {
    amatrix:::.amatrix_tripwire_record(
      "materialize_dense.eager",
      "mlx",
      paste0("key-", i),
      c(4L, 4L)
    )
  }

  expect_identical(amatrix:::.amatrix_tripwire_count(), 5L)
  events <- amatrix:::.amatrix_tripwire_events()
  expect_length(events, 5L)
  expect_true(all(vapply(events, function(e) identical(e$backend, "mlx"), logical(1))))
})
