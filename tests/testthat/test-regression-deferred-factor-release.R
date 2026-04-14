# Regression repro metadata
# Seed: 20260414
# Dimensions: dense SPD 4 x 4 with 4 x 2 RHS
# Backend / precision / dispatch: mock resident chol backend / fast / resident -> release -> host fallback
# R version / platform: R 4.5.1, aarch64-apple-darwin20 (macOS Sonoma 14.3)
# Issue: amatrix-ax8

suppressPackageStartupMessages(library(amatrix))

.mock_chol_release_backend <- function(counter) {
  backend <- make_recording_backend(
    counter,
    supported_ops = c("chol"),
    cold_supported_ops = character(),
    resident_supported_ops = c("chol"),
    precision_modes = "fast"
  )

  backend$supports <- function(op, x, y = NULL) {
    inherits(x, "adgeMatrix") &&
      identical(op, "chol") &&
      identical(x@preferred_backend, "mockcholrelease") &&
      identical(x@precision, "fast")
  }

  backend$supports_resident <- backend$supports

  backend$chol_resident <- function(lhs_key, out_key) {
    value <- chol(backend$resident_materialize(lhs_key))
    backend$resident_store(out_key, value)
    value
  }

  backend
}

test_that("resident chol factors keep a usable host copy after resident release [amatrix-ax8]", {
  counter <- new.env(parent = emptyenv())
  backend <- .mock_chol_release_backend(counter)

  with_registered_backend("mockcholrelease", backend, {
    set.seed(20260414L)
    A <- crossprod(matrix(rnorm(32L), nrow = 8L, ncol = 4L)) + diag(4L)
    B <- matrix(rnorm(8L), nrow = 4L, ncol = 2L)
    X <- adgeMatrix(A, preferred_backend = "mockcholrelease", precision = "fast")

    fac <- chol_factor(X)
    expect_identical(fac@backend, "mockcholrelease")
    expect_true(inherits(fac@factor_obj, "adgeMatrix"))
    expect_false(isTRUE(fac@factor_obj@finalizer_env$host_deferred))
    expect_identical(length(fac@factor), 0L)
    expect_identical(amatrix:::.amatrix_live_resident_backend(fac@factor_obj), "mockcholrelease")

    amatrix:::.amatrix_release_resident(fac@factor_obj)
    expect_null(amatrix:::.amatrix_live_resident_backend(fac@factor_obj))
    expect_equal(chol_solve(fac, B), solve(A, B), tolerance = 1e-10)
    expect_equal(as.matrix(fac), chol(A), tolerance = 1e-10)
  })
})
