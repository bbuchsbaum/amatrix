# Regression repro for amatrix-qm2.
# Seed: deterministic literal matrices (no RNG required)
# Shape: 2 x 2 dense matrices
# Backend: qm2_backend_a / qm2_backend_b
# Precision mode: strict
# Dispatch path: mixed resident matmul promotes rhs to lhs backend
# R/platform: R version 4.5.1 (2025-06-13) | aarch64-apple-darwin20
# Issue: amatrix-qm2

test_that("amatrix-qm2: mixed-backend resident prep does not overwrite tracked binding", {
  counter_a <- new.env(parent = emptyenv())
  counter_b <- new.env(parent = emptyenv())
  backend_a <- make_recording_backend(
    counter_a,
    supported_ops = c("matmul"),
    cold_supported_ops = c("matmul"),
    resident_supported_ops = c("matmul")
  )
  backend_b <- make_recording_backend(
    counter_b,
    supported_ops = c("matmul"),
    cold_supported_ops = c("matmul"),
    resident_supported_ops = c("matmul")
  )

  resident_a <- environment(backend_a$resident_has)$resident
  resident_b <- environment(backend_b$resident_has)$resident

  with_registered_backend("qm2_backend_a", backend_a, {
    with_registered_backend("qm2_backend_b", backend_b, {
      A <- adgeMatrix(matrix(c(1, 2, 3, 4), 2, 2), preferred_backend = "qm2_backend_a")
      B <- adgeMatrix(matrix(c(5, 6, 7, 8), 2, 2), preferred_backend = "qm2_backend_b")

      A <- amatrix_bind_resident(A, "qm2_backend_a")
      B <- amatrix_bind_resident(B, "qm2_backend_b")

      b_key_before <- amatrix:::.amatrix_resident_key(B, backend = "qm2_backend_b")
      expect_true(isTRUE(backend_b$resident_has(b_key_before)))
      expect_identical(amatrix:::.amatrix_live_resident_backend(B), "qm2_backend_b")

      result <- A %*% B

      expect_s4_class(result, "adgeMatrix")
      expect_equal(
        as.matrix(result),
        as.matrix(A) %*% as.matrix(B),
        tolerance = 1e-10
      )

      expect_identical(
        amatrix:::.amatrix_live_resident_backend(B),
        "qm2_backend_b",
        label = "mixed-backend promotion must not retarget B's tracked residency"
      )
      expect_identical(
        amatrix:::.amatrix_resident_key(B, backend = "qm2_backend_b"),
        b_key_before
      )
      expect_true(isTRUE(backend_b$resident_has(b_key_before)))
      expect_false(
        exists(b_key_before, envir = resident_a, inherits = FALSE),
        label = "the rhs original key must not appear in the lhs backend store"
      )
      expect_true(
        length(ls(envir = resident_a, all.names = FALSE)) >= 2L,
        label = "lhs backend should still allocate temporary/result keys"
      )
      expect_true(
        exists(b_key_before, envir = resident_b, inherits = FALSE),
        label = "rhs original backend key must remain live after mixed op"
      )
    })
  })
})
