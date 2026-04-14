# Regression repro for amatrix-v6y.
# Seed: deterministic literal matrix (no RNG required)
# Shape: 2 x 3 dense matrix
# Backend: regression_deferred_dimnames
# Precision mode: strict
# Dispatch path: internal deferred bind -> first host materialization
# R/platform: R version 4.5.1 (2025-06-13) | aarch64-apple-darwin20
# Issue: amatrix-v6y

test_that("amatrix-v6y: deferred host materialization preserves dimnames", {
  counter <- new.env(parent = emptyenv())
  backend <- make_recording_backend(
    counter,
    supported_ops = character(),
    cold_supported_ops = character(),
    resident_supported_ops = character()
  )

  resident_env <- environment(backend$resident_has)$resident
  backend$resident_materialize <- function(key) {
    if (is.null(counter$resident_materialize)) {
      counter$resident_materialize <- 0L
    }
    counter$resident_materialize <- counter$resident_materialize + 1L
    mat <- get(key, envir = resident_env, inherits = FALSE)
    dimnames(mat) <- list(NULL, NULL)
    mat
  }

  with_registered_backend("regression_deferred_dimnames", backend, {
    host <- matrix(
      c(1, 2, 3, 4, 5, 6),
      nrow = 2L,
      ncol = 3L,
      dimnames = list(c("row_a", "row_b"), c("col_x", "col_y", "col_z"))
    )

    x <- adgeMatrix(host, preferred_backend = "regression_deferred_dimnames")
    bound <- amatrix_bind_resident(x, "regression_deferred_dimnames")
    deferred <- amatrix:::new_adgeMatrix_deferred(
      dim = as.integer(dim(host)),
      dimnames = dimnames(host),
      preferred_backend = "regression_deferred_dimnames"
    )
    deferred <- amatrix:::.amatrix_bind_resident(
      deferred,
      "regression_deferred_dimnames",
      amatrix:::.amatrix_resident_key(bound)
    )

    host_result <- as.matrix(deferred)

    expect_equal(host_result, host, tolerance = 1e-10)
    expect_identical(dimnames(host_result), dimnames(host))
    expect_identical(counter$resident_materialize, 1L)
  })
})
