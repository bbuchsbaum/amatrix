# Regression repro metadata
# Seed: none (deterministic zero-dimension fixtures)
# Dimensions: 0 x 3 and 2 x 0 matrices
# Backend / precision / dispatch: auto metric wrapper path
# R version / platform: captured by CI sessionInfo() on failure
# Issue: amatrix-c6v

test_that("zero-dimension metric wrappers do not crash in auto mode", {
  skip_if_not_installed("callr")
  skip_if_not_installed("pkgload")

  result <- callr::r(
    function(path) {
      pkgload::load_all(path, quiet = TRUE)

      x0r <- matrix(numeric(0), nrow = 0L, ncol = 3L)
      y0r <- matrix(numeric(0), nrow = 0L, ncol = 3L)
      x0c <- matrix(numeric(0), nrow = 2L, ncol = 0L)

      list(
        dist_0r = dim(dist_matrix(x0r)),
        dist_0r_cross = dim(dist_matrix(x0r, y0r)),
        kernel_0r = dim(kernel_matrix(x0r)),
        kernel_0r_cross = dim(kernel_matrix(x0r, y0r)),
        dist_0c = dist_matrix(x0c),
        kernel_0c_linear = kernel_matrix(x0c, kernel = "linear"),
        kernel_0c_cosine = kernel_matrix(x0c, kernel = "cosine")
      )
    },
    args = list(getwd()),
    spinner = FALSE,
    show = FALSE
  )

  expect_identical(result$dist_0r, c(0L, 0L))
  expect_identical(result$dist_0r_cross, c(0L, 0L))
  expect_identical(result$kernel_0r, c(0L, 0L))
  expect_identical(result$kernel_0r_cross, c(0L, 0L))
  expect_equal(result$dist_0c, matrix(0, nrow = 2L, ncol = 2L), tolerance = 0)
  expect_equal(result$kernel_0c_linear, matrix(0, nrow = 2L, ncol = 2L), tolerance = 0)
  expect_equal(result$kernel_0c_cosine, matrix(0, nrow = 2L, ncol = 2L), tolerance = 0)
})
