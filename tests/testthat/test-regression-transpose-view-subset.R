# Regression repro metadata
# Seed: none (deterministic literal matrix)
# Dimensions: source 3 x 4, transpose view 4 x 3
# Backend / precision / dispatch: cpu / strict / aTransposeView subset path
# R version / platform: R 4.5.1, aarch64-apple-darwin20 (macOS Sonoma 14.3)
# Issue: amatrix-xnp

suppressPackageStartupMessages(library(amatrix))

test_that("aTransposeView supports scalar and matrix subsetting [amatrix-xnp]", {
  host <- matrix(1:12, nrow = 3L, ncol = 4L)
  tv <- t(adgeMatrix(host))
  ref <- t(host)

  expect_equal(tv[1L, 2L], ref[1L, 2L], tolerance = 0)
  expect_equal(tv[, 2L], ref[, 2L], tolerance = 0)

  block <- tv[1:3, 2:3, drop = FALSE]
  expect_s4_class(block, "adgeMatrix")
  expect_identical(block@preferred_backend, tv@preferred_backend)
  expect_identical(block@policy, tv@policy)
  expect_identical(block@precision, tv@precision)
  expect_equal(as.matrix(block), ref[1:3, 2:3, drop = FALSE], tolerance = 0)
})
