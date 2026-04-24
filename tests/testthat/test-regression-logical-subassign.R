# Regression repro metadata
# Seed: none (deterministic dense fixture)
# Dimensions: 2 x 2 dense matrix with one NA
# Backend / precision / dispatch: cpu / strict / one-subscript logical replacement
# R version / platform: captured by CI sessionInfo() on failure
# Issue: amatrix-e97

test_that("one-subscript logical matrix replacement only updates selected cells [amatrix-e97]", {
  host <- matrix(c(1, NA, 3, 4), nrow = 2)
  expected <- host
  expected[is.na(expected)] <- 0

  x <- adgeMatrix(host, preferred_backend = "cpu")
  x[is.na(x)] <- 0

  expect_s4_class(x, "adgeMatrix")
  expect_equal(as.matrix(x), expected, tolerance = 0)
})
