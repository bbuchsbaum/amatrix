# Regression repro metadata
# Seed: none (deterministic deferred fixture)
# Dimensions: 2 x 3 deferred dense matrix
# Backend / precision / dispatch: cpu / strict / saveRDS-readRDS round-trip
# R version / platform: captured by CI sessionInfo() on failure
# Issues: amatrix-1i1

suppressPackageStartupMessages(library(amatrix))

test_that("dead deferred roundtrip is detected explicitly and fails consistently [amatrix-1i1]", {
  x <- amatrix:::new_adgeMatrix_deferred(dim = c(2L, 3L), preferred_backend = "cpu")
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)

  saveRDS(x, tmp)
  x2 <- readRDS(tmp)

  expect_true(amatrix:::.amatrix_is_dead_deferred(x2))

  expect_error(
    show(x2),
    "deferred adgeMatrix cannot survive saveRDS/readRDS without host materialization; GPU resident data is unavailable"
  )
  expect_error(
    as.matrix(x2),
    "deferred adgeMatrix cannot survive saveRDS/readRDS without host materialization; GPU resident data is unavailable"
  )
  expect_error(
    as.numeric(x2),
    "deferred adgeMatrix cannot survive saveRDS/readRDS without host materialization; GPU resident data is unavailable"
  )
  expect_error(
    as.vector(x2),
    "deferred adgeMatrix cannot survive saveRDS/readRDS without host materialization; GPU resident data is unavailable"
  )
})
