# Regression repro metadata
# Seed: none (deterministic fixture with NA values)
# Dimensions: 3 x 4 dense matrix
# Backend / precision / dispatch: cpu / strict / cold reduction path
# R version / platform: captured by CI sessionInfo() on failure
# Issue: amatrix-fya

test_that("rowmeans(adgeMatrix, na.rm=TRUE) matches base rowMeans", {
  host <- matrix(
    c(1, NA, 3,
      4, 5, NA,
      7, 8, 9,
      10, NA, 12),
    nrow = 3L,
    ncol = 4L
  )
  x <- as_adgeMatrix(
    host,
    preferred_backend = "cpu",
    policy = "opencl",
    precision = "strict"
  )

  expect_equal(
    rowmeans(x, na.rm = TRUE),
    base::rowMeans(host, na.rm = TRUE),
    tolerance = 1e-12
  )
})
