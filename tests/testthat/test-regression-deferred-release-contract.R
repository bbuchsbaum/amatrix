# Regression repro for amatrix-cth.
# Seed: deterministic literal matrix (no RNG required)
# Shape: 2 x 3 dense matrix
# Backend: cth_backend (recording backend)
# Precision mode: strict
# Dispatch path: deferred adgeMatrix residency release + product plan reuse
# R/platform: R version 4.5.1 (2025-06-13) | aarch64-apple-darwin20
# Issue: amatrix-cth
#
# A deferred adgeMatrix has no authoritative host copy (host @x is NaN
# placeholders); its only real data lives in the device buffer. The documented
# contract of amatrix_release_resident() promises "the object remains fully
# usable afterwards: its data is served from the host copy". For a deferred
# object that promise is broken unless release materializes the host copy first
# — otherwise host access throws (dead-deferred) or reads NaN. amatrix-cth
# observes this through a product plan that keeps computing after the operand
# is released.

test_that("amatrix-cth: releasing a deferred adgeMatrix keeps it usable via host copy", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "cth_backend",
    make_recording_backend(
      counter,
      supported_ops = c("matmul"),
      cold_supported_ops = c("matmul"),
      resident_supported_ops = c("matmul")
    ),
    {
      data <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2L, ncol = 3L)
      A <- adgeMatrix(data, preferred_backend = "cth_backend")
      h <- resident_handle(A, backend = "cth_backend")
      d <- amatrix:::as_adgeMatrix.resident_handle(h, defer_host = TRUE)

      # Setup sanity: d is genuinely deferred with no host copy yet.
      expect_true(isTRUE(d@finalizer_env$host_deferred))
      expect_null(d@finalizer_env$host_x)

      # Release device memory WITHOUT touching the host first.
      released <- amatrix_release_resident(d)
      expect_true(isTRUE(released))

      # Documented contract: the object remains fully usable afterwards.
      host <- as.matrix(amatrix_materialize_host(d))
      expect_equal(host, data, ignore_attr = TRUE)
    }
  )
})

test_that("amatrix-cth: a product plan keeps computing after its operand is released", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "cth_plan_backend",
    make_recording_backend(
      counter,
      supported_ops = c("matmul"),
      cold_supported_ops = c("matmul"),
      resident_supported_ops = c("matmul")
    ),
    {
      data <- matrix(c(2, 0, 1, 0, 3, 1), nrow = 2L, ncol = 3L)
      A <- adgeMatrix(data, preferred_backend = "cth_plan_backend")
      h <- resident_handle(A, backend = "cth_plan_backend")
      d <- amatrix:::as_adgeMatrix.resident_handle(h, defer_host = TRUE)

      plan <- amatrix_compile_product(d, op = "matmul", backend = "cth_plan_backend")
      rhs <- matrix(c(1, 0, 0, 1, 1, 1), nrow = 3L, ncol = 2L)

      before <- plan(rhs, materialize = "matrix")
      expect_equal(before, data %*% rhs, ignore_attr = TRUE)

      amatrix_release_resident(d)

      after <- plan(rhs, materialize = "matrix")
      expect_equal(after, data %*% rhs, ignore_attr = TRUE)
    }
  )
})
