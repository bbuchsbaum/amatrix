# Regression repro metadata
# Seed: deterministic literal matrix (no RNG required)
# Dimensions: 3 x 3 dense matrix
# Backend / precision / dispatch: planning_cleanup_backend / strict / resident dispatch fallback
# R version / platform: R 4.5.1, aarch64-apple-darwin20 (macOS Sonoma 14.3)
# Issue: amatrix-aul

suppressPackageStartupMessages(library(amatrix))

test_that("resident dispatch drops out_key even if lhs cleanup throws [amatrix-aul]", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "planning_cleanup_backend",
    make_recording_backend(
      counter,
      supported_ops = c("chol"),
      cold_supported_ops = c("chol"),
      resident_supported_ops = c("chol")
    ),
    {
      backend <- get(
        "planning_cleanup_backend",
        envir = amatrix:::.amatrix_state$backends,
        inherits = FALSE
      )
      resident <- environment(backend$resident_has)$resident
      backend$chol_resident <- function(lhs_key, out_key) {
        assign(out_key, diag(3L), envir = resident)
        stop("resident chol boom", call. = FALSE)
      }
      assign("planning_cleanup_backend", backend, envir = amatrix:::.amatrix_state$backends)

      x <- adgeMatrix(diag(3L), preferred_backend = "planning_cleanup_backend")
      x <- amatrix_bind_resident(x, "planning_cleanup_backend")
      x_key <- amatrix:::.amatrix_resident_key(x, backend = "planning_cleanup_backend")

      local_mocked_bindings(
        .amatrix_cleanup_temp_resident_safe = function(...) {
          stop("cleanup boom", call. = FALSE)
        },
        .package = "amatrix"
      )

      expect_equal(unname(as.matrix(amatrix_dispatch_op(
        x = x,
        op = "chol",
        method = "chol",
        fallback = function() "fallback"
      ))), diag(3L), tolerance = 0)

      expect_identical(sort(ls(envir = resident, all.names = FALSE)), x_key)
      expect_true(counter$resident_drop >= 1L)
    }
  )
})
