test_that("deferred dense coercions fail loudly when resident data is missing", {
  x <- amatrix:::new_adgeMatrix_deferred(c(2L, 2L))

  err <- "deferred adgeMatrix cannot survive serialization \\(saveRDS/readRDS\\) without host materialization; GPU resident data is unavailable"
  expect_error(as.matrix(x), err)
  expect_error(as.numeric(x), err)
  expect_error(as.vector(x), err)
})
