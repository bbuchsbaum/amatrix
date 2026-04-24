# Regression repro metadata
# Seed: none (deterministic deferred fixture)
# Dimensions: 2 x 3 deferred dense matrix
# Backend / precision / dispatch: cpu / strict / saveRDS-readRDS round-trip
# R version / platform: R 4.5.1, aarch64-apple-darwin20 (macOS Sonoma 14.3)
# Issue: amatrix-90k

suppressPackageStartupMessages(library(amatrix))

test_that("deferred adgeMatrix roundtrip has an explicit S3 as.matrix bridge and errors cleanly [amatrix-90k]", {
  x <- amatrix:::new_adgeMatrix_deferred(dim = c(2L, 3L), preferred_backend = "cpu")
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  x2 <- readRDS(tmp)

  expect_true(is.function(getS3method("as.matrix", "adgeMatrix", optional = TRUE)))
  err <- "deferred adgeMatrix cannot survive serialization \\(saveRDS/readRDS\\) without host materialization; GPU resident data is unavailable"
  expect_error(as.matrix(x2), err)
  expect_error(base::as.matrix(x2), err)
})
