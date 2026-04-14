# Regression repro metadata
# Seed: none (deterministic mocked subspace error)
# Dimensions: dense 6 x 4
# Backend / precision / dispatch: cpu / strict / explicit subspace path with mocked error
# R version / platform: R 4.5.1, aarch64-apple-darwin20 (macOS Sonoma 14.3)
# Issue: amatrix-e4w

suppressPackageStartupMessages(library(amatrix))

test_that("svd_factor does not silently fall back when subspace path signals an error [amatrix-e4w]", {
  x <- adgeMatrix(matrix(rnorm(24L), nrow = 6L, ncol = 4L))

  local_mocked_bindings(
    am_svd = function(...) {
      stop("exact fallback should not run")
    },
    rsvd = function(...) {
      stop("rsvd fallback should not run")
    },
    .amatrix_subspace_svd = function(...) {
      stop(errorCondition(
        "forced subspace failure",
        class = "amatrix_subspace_error",
        call = NULL
      ))
    },
    .package = "amatrix"
  )

  expect_error(
    svd_factor(x, k = 2L, method = "subspace"),
    class = "amatrix_subspace_error"
  )
})
