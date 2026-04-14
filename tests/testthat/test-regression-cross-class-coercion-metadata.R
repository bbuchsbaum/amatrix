# Regression repro metadata
# Seed: none (deterministic fixture)
# Dimensions: dense 3 x 4, sparse 3 x 4
# Backend / precision / dispatch: cpu / strict+fast / coercion round-trip
# R version / platform: captured by CI sessionInfo() on failure
# Issue: amatrix-3su

test_that("sparse -> dgCMatrix -> adgeMatrix preserves metadata", {
  sp_host <- Matrix::sparseMatrix(
    i = c(1L, 2L, 3L),
    j = c(1L, 2L, 4L),
    x = c(1, 2, 3),
    dims = c(3L, 4L)
  )
  x <- as_adgCMatrix(
    sp_host,
    preferred_backend = "cpu",
    policy = "opencl",
    precision = "fast"
  )

  roundtrip <- methods::as(methods::as(x, "dgCMatrix"), "adgeMatrix")
  expect_identical(roundtrip@preferred_backend, x@preferred_backend)
  expect_identical(roundtrip@policy, x@policy)
  expect_identical(roundtrip@precision, x@precision)
  expect_equal(as.matrix(roundtrip), as.matrix(x), tolerance = 1e-12)
})

test_that("dense -> dgeMatrix -> adgCMatrix preserves metadata", {
  x <- as_adgeMatrix(
    matrix(seq_len(12L), nrow = 3L, ncol = 4L),
    preferred_backend = "cpu",
    policy = "opencl",
    precision = "strict"
  )

  roundtrip <- methods::as(methods::as(x, "dgeMatrix"), "adgCMatrix")
  expect_identical(roundtrip@preferred_backend, x@preferred_backend)
  expect_identical(roundtrip@policy, x@policy)
  expect_identical(roundtrip@precision, x@precision)
  expect_equal(as.matrix(roundtrip), as.matrix(x), tolerance = 1e-12)
})
