# Regression repro metadata
# Seed: none (deterministic fixture with NA values)
# Dimensions: 4 x 3 dense matrix
# Backend / precision / dispatch: cpu / strict / host rowsum path
# R version / platform: captured by CI sessionInfo() on failure
# Issue: amatrix-adu

test_that("rowsum.adgeMatrix honours na.rm=TRUE", {
  host <- matrix(
    c(1, 2, NA, 4,
      5, NA, 7, 8,
      9, 10, 11, NA),
    nrow = 4L,
    ncol = 3L
  )
  groups <- c("a", "a", "b", "b")
  x <- as_adgeMatrix(
    host,
    preferred_backend = "cpu",
    policy = "opencl",
    precision = "strict"
  )

  expect_equal(
    rowsum.adgeMatrix(x, groups, na.rm = TRUE),
    base::rowsum(host, groups, na.rm = TRUE),
    tolerance = 1e-12
  )
})
