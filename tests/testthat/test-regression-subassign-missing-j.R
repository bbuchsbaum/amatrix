# Regression repro metadata
# Seed: none (deterministic fixture)
# Dimensions: 3 x 4 dense matrix
# Backend / precision / dispatch: cpu / strict / cold host fallback
# R version / platform: captured by CI sessionInfo() on failure
# Issue: amatrix-d8s

test_that("row replacement with missing j preserves matrix semantics", {
  x <- as_adgeMatrix(
    matrix(seq_len(12L), nrow = 3L, ncol = 4L),
    preferred_backend = "cpu",
    policy = "opencl",
    precision = "strict"
  )
  expected <- matrix(seq_len(12L), nrow = 3L, ncol = 4L)
  expected[1L, ] <- c(10, 11, 12, 13)

  expect_no_warning({
    x[1L, ] <- c(10, 11, 12, 13)
  })
  expect_s4_class(x, "adgeMatrix")
  expect_identical(x@preferred_backend, "cpu")
  expect_identical(x@policy, "opencl")
  expect_equal(as.matrix(x), expected, tolerance = 1e-12)
})
