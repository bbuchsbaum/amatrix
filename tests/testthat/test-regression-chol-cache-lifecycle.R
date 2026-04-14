# Regression repro metadata
# Seed: 20260414
# Dimensions: dense SPD 5 x 5
# Backend / precision / dispatch: cpu / strict / host-cold chol cache
# R version / platform: R 4.5.1, aarch64-apple-darwin20 (macOS Sonoma 14.3)
# Issue: amatrix-tqm

suppressPackageStartupMessages(library(amatrix))

.chol_cache_keys <- function() {
  sort(ls(envir = amatrix:::.amatrix_state$model_cache, all.names = FALSE))
}

test_that("chol cache invalidates stale entries across repeated subassignment mutation [amatrix-tqm]", {
  old_max <- amatrix_cache_max_size()
  on.exit(amatrix_set_cache_max_size(old_max), add = TRUE)
  amatrix_set_cache_max_size(Inf)
  amatrix:::.amatrix_cache_clear()

  set.seed(20260414L)
  x <- as_adgeMatrix(crossprod(matrix(rnorm(25L), nrow = 5L, ncol = 5L)) + diag(5L))

  for (step in seq_len(6L)) {
    fac <- chol_factor(x)
    keys <- .chol_cache_keys()
    expect_identical(keys, paste0("chol:", x@object_id), info = paste("step", step))
    expect_identical(fac@source_id, x@object_id, info = paste("step", step))

    x[1L, 1L] <- x[1L, 1L] + step / 10
  }

  final_fac <- chol_factor(x)
  expect_identical(.chol_cache_keys(), paste0("chol:", x@object_id))
  expect_equal(as.matrix(final_fac), chol(as.matrix(x)), tolerance = 1e-10)
})

test_that("chol cache invalidates stale entries across diag replacement [amatrix-tqm]", {
  amatrix:::.amatrix_cache_clear()

  set.seed(202604141L)
  x <- as_adgeMatrix(crossprod(matrix(rnorm(16L), nrow = 4L, ncol = 4L)) + diag(4L))
  chol_factor(x)
  expect_identical(.chol_cache_keys(), paste0("chol:", x@object_id))

  diag(x) <- diag(x) + 0.5
  fac <- chol_factor(x)

  expect_identical(.chol_cache_keys(), paste0("chol:", x@object_id))
  expect_equal(as.matrix(fac), chol(as.matrix(x)), tolerance = 1e-10)
})
