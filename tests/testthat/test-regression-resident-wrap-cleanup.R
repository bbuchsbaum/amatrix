# Regression repro for amatrix-aul.
# Seed: deterministic literal matrix (no RNG required)
# Shape: 3 x 3 dense matrix
# Backend: aul_wrap_backend
# Precision mode: strict
# Dispatch path: resident crossprod -> resident wrap
# R/platform: R version 4.5.1 (2025-06-13) | aarch64-apple-darwin20
# Issue: amatrix-aul

test_that("amatrix-aul: resident output key is dropped when resident wrap fails", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "aul_wrap_backend",
    make_recording_backend(
      counter,
      supported_ops = c("crossprod"),
      cold_supported_ops = c("crossprod"),
      resident_supported_ops = c("crossprod")
    ),
    {
      backend <- get(
        "aul_wrap_backend",
        envir = amatrix:::.amatrix_state$backends,
        inherits = FALSE
      )
      resident <- environment(backend$resident_has)$resident

      x <- adgeMatrix(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3L),
        preferred_backend = "aul_wrap_backend"
      )
      x <- amatrix_bind_resident(x, "aul_wrap_backend")
      x_key <- amatrix:::.amatrix_resident_key(x, backend = "aul_wrap_backend")

      local_mocked_bindings(
        .amatrix_rewrap_like = function(...) {
          stop("resident wrap boom", call. = FALSE)
        },
        .package = "amatrix"
      )

      expect_error(am_crossprod(x), "resident wrap boom")
      expect_identical(sort(ls(envir = resident, all.names = FALSE)), x_key)
      expect_true(counter$crossprod_resident >= 1L)
      expect_true(counter$resident_drop >= 1L)
    }
  )
})
