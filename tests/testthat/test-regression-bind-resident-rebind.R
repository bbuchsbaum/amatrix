# Regression repro for amatrix-3ka.
# Seed: deterministic literal matrices (no RNG required)
# Shape: 2 x 2 dense matrices
# Backend: rebind_backend (recording backend)
# Precision mode: strict
# Dispatch path: .amatrix_bind_resident rebind (registry overwrite)
# R/platform: R version 4.5.1 (2025-06-13) | aarch64-apple-darwin20
# Issue: amatrix-3ka
#
# .amatrix_bind_resident overwrites an object's residency-registry entry. If the
# object already has a live entry pointing at key K1 and it is rebound to a new
# key K2, K1 must be released — otherwise the registry (and the finalizer) only
# tracks K2 and K1 leaks in device memory until process exit. A shared/aliased
# K1 (still referenced by a second object) must NOT be dropped.

test_that("amatrix-3ka: rebinding an object releases the prior resident key", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "rebind_backend",
    make_recording_backend(
      counter,
      supported_ops = c("matmul"),
      cold_supported_ops = c("matmul"),
      resident_supported_ops = c("matmul")
    ),
    {
      backend <- get(
        "rebind_backend",
        envir = amatrix:::.amatrix_state$backends,
        inherits = FALSE
      )

      x <- adgeMatrix(matrix(c(1, 2, 3, 4), nrow = 2L),
        preferred_backend = "rebind_backend"
      )

      k1 <- amatrix:::.amatrix_next_resident_key("rebind_backend")
      backend$resident_store(k1, matrix(c(1, 2, 3, 4), nrow = 2L))
      amatrix:::.amatrix_bind_resident(x, "rebind_backend", k1)

      k2 <- amatrix:::.amatrix_next_resident_key("rebind_backend")
      backend$resident_store(k2, matrix(c(5, 6, 7, 8), nrow = 2L))
      amatrix:::.amatrix_bind_resident(x, "rebind_backend", k2)

      # K1 must have been released on rebind; K2 is the live binding.
      expect_false(backend$resident_has(k1))
      expect_true(backend$resident_has(k2))
      expect_identical(amatrix:::.amatrix_resident_key(x, "rebind_backend"), k2)
    }
  )
})

test_that("amatrix-3ka: rebinding does not drop a key aliased by another object", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "rebind_alias_backend",
    make_recording_backend(
      counter,
      supported_ops = c("matmul"),
      cold_supported_ops = c("matmul"),
      resident_supported_ops = c("matmul")
    ),
    {
      backend <- get(
        "rebind_alias_backend",
        envir = amatrix:::.amatrix_state$backends,
        inherits = FALSE
      )

      x <- adgeMatrix(matrix(c(1, 2, 3, 4), nrow = 2L),
        preferred_backend = "rebind_alias_backend"
      )
      y <- adgeMatrix(matrix(c(1, 2, 3, 4), nrow = 2L),
        preferred_backend = "rebind_alias_backend"
      )

      shared <- amatrix:::.amatrix_next_resident_key("rebind_alias_backend")
      backend$resident_store(shared, matrix(c(1, 2, 3, 4), nrow = 2L))
      amatrix:::.amatrix_bind_resident(x, "rebind_alias_backend", shared)
      amatrix:::.amatrix_bind_resident(y, "rebind_alias_backend", shared)

      # Rebind x to a fresh key; the shared key is still used by y and must live.
      k2 <- amatrix:::.amatrix_next_resident_key("rebind_alias_backend")
      backend$resident_store(k2, matrix(c(5, 6, 7, 8), nrow = 2L))
      amatrix:::.amatrix_bind_resident(x, "rebind_alias_backend", k2)

      expect_true(backend$resident_has(shared))
      expect_true(backend$resident_has(k2))
      expect_identical(amatrix:::.amatrix_resident_key(y, "rebind_alias_backend"), shared)
    }
  )
})
