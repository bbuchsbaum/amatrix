# Regression repro metadata
# Seed: none (deterministic fixtures)
# Dimensions: dense 2 x 3, sparse 2 x 3
# Backend / precision / dispatch: cpu / strict / host subset path
# R version / platform: captured by CI sessionInfo() on failure
# Issue: amatrix-36t

test_that("dense subset with missing column or row index matches host semantics", {
  host <- matrix(1:6, nrow = 2L, ncol = 3L)
  x <- as_adgeMatrix(
    host,
    preferred_backend = "cpu",
    policy = "opencl",
    precision = "strict"
  )

  expect_equal(x[1L, ], host[1L, ], tolerance = 0)
  expect_equal(x[, 2L], host[, 2L], tolerance = 0)
  expect_equal(as.matrix(x[1L, , drop = FALSE]), host[1L, , drop = FALSE], tolerance = 0)
  expect_equal(as.matrix(x[, 2L, drop = FALSE]), host[, 2L, drop = FALSE], tolerance = 0)
})

test_that("sparse subset with missing column or row index matches host semantics", {
  host <- Matrix::sparseMatrix(
    i = c(1L, 2L),
    j = c(1L, 3L),
    x = c(1, 2),
    dims = c(2L, 3L)
  )
  x <- as_adgCMatrix(
    host,
    preferred_backend = "cpu",
    policy = "opencl",
    precision = "strict"
  )

  expect_equal(as.matrix(x[1L, ]), as.matrix(host[1L, ]), tolerance = 0)
  expect_equal(as.matrix(x[, 3L]), as.matrix(host[, 3L]), tolerance = 0)
  expect_equal(as.matrix(x[1L, , drop = FALSE]), as.matrix(host[1L, , drop = FALSE]), tolerance = 0)
  expect_equal(as.matrix(x[, 3L, drop = FALSE]), as.matrix(host[, 3L, drop = FALSE]), tolerance = 0)
})
