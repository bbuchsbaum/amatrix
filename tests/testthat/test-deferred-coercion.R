test_that("deferred dense coercions fail loudly when resident data is missing", {
  x <- amatrix:::new_adgeMatrix_deferred(c(2L, 2L))

  expect_error(as.matrix(x), "deferred adgeMatrix lost its GPU resident data")
  expect_error(as.numeric(x), "deferred adgeMatrix lost its GPU resident data")
  expect_error(as.vector(x), "deferred adgeMatrix lost its GPU resident data")
})
