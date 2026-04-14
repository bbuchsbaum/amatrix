# Regression repro metadata
# Seed: none (deterministic fixtures)
# Dimensions: dense 2 x 3, sparse 2 x 3
# Backend / precision / dispatch: cpu / strict / subset dispatch path
# R version / platform: captured by CI sessionInfo() on failure
# Issue: amatrix-tva

test_that("default-drop dense subset keeps adgeMatrix when result stays matrix-like", {
  host <- matrix(1:6, nrow = 2L, ncol = 3L)
  x <- as_adgeMatrix(
    host,
    preferred_backend = "cpu",
    policy = "opencl",
    precision = "strict"
  )

  out <- x[, 1:2]
  expect_s4_class(out, "adgeMatrix")
  expect_identical(out@preferred_backend, x@preferred_backend)
  expect_identical(out@policy, x@policy)
  expect_identical(out@precision, x@precision)
  expect_equal(as.matrix(out), host[, 1:2], tolerance = 0)
})

test_that("default-drop sparse subset keeps adgCMatrix when result stays sparse matrix-like", {
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

  out <- x[, 1:2]
  expect_s4_class(out, "adgCMatrix")
  expect_identical(out@preferred_backend, x@preferred_backend)
  expect_identical(out@policy, x@policy)
  expect_identical(out@precision, x@precision)
  expect_equal(as.matrix(out), as.matrix(host[, 1:2]), tolerance = 0)
})
