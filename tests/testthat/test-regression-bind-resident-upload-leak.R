# Regression repro for amatrix-aul (bind-resident upload remainder).
# Seed: deterministic literal matrix (no RNG required)
# Shape: 2 x 2 dense matrix
# Backend: aul_bind_backend (recording backend)
# Precision mode: strict
# Dispatch path: amatrix_bind_resident upload -> residency binding
# R/platform: R version 4.5.1 (2025-06-13) | aarch64-apple-darwin20
# Issue: amatrix-aul
#
# amatrix_bind_resident() allocates a resident key, uploads the host matrix to
# the device, and then records the residency binding. If the binding step
# throws after the upload, the device buffer has no registry entry and no
# finalizer, so it leaks until process exit. The fix registers an on.exit that
# releases the key unless the binding was recorded.

test_that("amatrix-aul: amatrix_bind_resident releases the key when binding fails", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "aul_bind_backend",
    make_recording_backend(
      counter,
      supported_ops = c("matmul"),
      cold_supported_ops = c("matmul"),
      resident_supported_ops = c("matmul")
    ),
    {
      backend <- get(
        "aul_bind_backend",
        envir = amatrix:::.amatrix_state$backends,
        inherits = FALSE
      )
      resident <- environment(backend$resident_has)$resident

      local_mocked_bindings(
        .amatrix_bind_resident = function(...) stop("bind boom", call. = FALSE),
        .package = "amatrix"
      )

      x <- matrix(c(1, 2, 3, 4), nrow = 2L)
      expect_error(amatrix_bind_resident(x, "aul_bind_backend"), "bind boom")

      # The uploaded buffer must be released; the resident env stays empty.
      expect_length(ls(envir = resident, all.names = TRUE), 0L)
      expect_true(isTRUE(counter$resident_drop >= 1L))
    }
  )
})
