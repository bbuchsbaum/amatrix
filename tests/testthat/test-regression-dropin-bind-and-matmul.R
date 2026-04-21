# Regression repro metadata
# Seed: none (deterministic dense fixtures)
# Dimensions: bind cases 3 x 2 with 3 x 1 / 1 x 2 addends; matmul case 3 x 2 by
#   plain length-2 vector
# Backend / precision / dispatch: cpu / strict / fresh in-session dense path
# R version / platform: captured by CI sessionInfo() on failure
# Issues: amatrix-0qt, amatrix-qic

test_that("mixed two-argument cbind/rbind preserve adgeMatrix wrapper [amatrix-0qt]", {
  x_host <- matrix(1:6 * 1.0, nrow = 3)
  x <- adgeMatrix(x_host, preferred_backend = "cpu")

  col_vec <- c(10, 11, 12)
  row_vec <- c(20, 21)
  col_mat <- matrix(col_vec, ncol = 1)
  row_mat <- matrix(row_vec, nrow = 1)

  cb_mat <- cbind(x, col_mat)
  cb_vec <- cbind(x, col_vec)
  rb_mat <- rbind(x, row_mat)
  rb_vec <- rbind(x, row_vec)

  expect_s4_class(cb_mat, "adgeMatrix")
  expect_s4_class(cb_vec, "adgeMatrix")
  expect_s4_class(rb_mat, "adgeMatrix")
  expect_s4_class(rb_vec, "adgeMatrix")

  expect_equal(as.matrix(cb_mat), cbind(x_host, col_mat), tolerance = 0)
  expect_equal(as.matrix(cb_vec), cbind(x_host, col_vec), tolerance = 0)
  expect_equal(as.matrix(rb_mat), rbind(x_host, row_mat), tolerance = 0)
  expect_equal(as.matrix(rb_vec), rbind(x_host, row_vec), tolerance = 0)
})

test_that("adgeMatrix %*% plain vector preserves matrix-like result [amatrix-qic]", {
  x_host <- matrix(1:6 * 1.0, nrow = 3)
  x <- adgeMatrix(x_host, preferred_backend = "cpu")
  v <- c(1, 2)

  result <- x %*% v
  reference <- x_host %*% v

  expect_s4_class(result, "adgeMatrix")
  expect_identical(dim(result), c(3L, 1L))
  expect_equal(as.matrix(result), reference, tolerance = 0)
})
