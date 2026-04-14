# Regression repro metadata
# Seed: none (deterministic fixture)
# Dimensions: 3 x 4 sparse matrix
# Backend / precision / dispatch: cpu / strict / cold host fallback
# R version / platform: captured by CI sessionInfo() on failure
# Issue: amatrix-932

test_that("sparse row replacement with missing j preserves adgCMatrix class", {
  x <- as_adgCMatrix(
    Matrix::sparseMatrix(
      i = c(1L, 2L, 3L),
      j = c(1L, 2L, 4L),
      x = c(1, 2, 3),
      dims = c(3L, 4L)
    ),
    preferred_backend = "cpu",
    policy = "opencl",
    precision = "strict"
  )
  expected <- as.matrix(x)
  expected[2L, ] <- c(7, 0, 8, 0)

  x[2L, ] <- c(7, 0, 8, 0)
  expect_s4_class(x, "adgCMatrix")
  expect_identical(x@preferred_backend, "cpu")
  expect_identical(x@policy, "opencl")
  expect_equal(as.matrix(x), expected, tolerance = 1e-12)
})
