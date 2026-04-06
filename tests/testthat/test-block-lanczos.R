low_rank_matrix <- function(n = 80L, p = 50L, rank = 5L) {
  set.seed(42L + n + p + rank)
  u <- qr.Q(qr(matrix(rnorm(n * rank), nrow = n, ncol = rank)))
  v <- qr.Q(qr(matrix(rnorm(p * rank), nrow = p, ncol = rank)))
  d <- seq(from = rank, to = 1, length.out = rank)
  u %*% diag(d, nrow = rank, ncol = rank) %*% t(v)
}

svd_reconstruct <- function(res) {
  res$u %*% diag(res$d, nrow = length(res$d), ncol = length(res$d)) %*% t(res$v)
}

test_that("am_block_lanczos recovers a low-rank matrix accurately", {
  x <- low_rank_matrix()
  ref <- La.svd(x, nu = 5L, nv = 5L)

  set.seed(1001L)
  res <- am_block_lanczos(x, nv = 5L, nu = 5L, block_size = 4L, n_steps = 6L)

  expect_identical(dim(res$u), c(nrow(x), 5L))
  expect_identical(dim(res$v), c(ncol(x), 5L))
  expect_equal(res$d, ref$d[seq_len(5L)], tolerance = 1e-6)

  rel_err <- norm(svd_reconstruct(res) - x, type = "F") / max(norm(x, type = "F"), 1e-12)
  expect_lt(rel_err, 1e-6)
})

test_that("am_block_svd remains an alias for am_block_lanczos", {
  x <- low_rank_matrix(n = 60L, p = 40L, rank = 4L)

  set.seed(1002L)
  via_alias <- am_block_svd(x, k = 4L, block_size = 4L, n_steps = 6L)
  set.seed(1002L)
  via_named <- am_block_lanczos(x, nv = 4L, nu = 4L, block_size = 4L, n_steps = 6L)

  expect_equal(via_alias$d, via_named$d, tolerance = 1e-10)
  expect_equal(svd_reconstruct(via_alias), svd_reconstruct(via_named), tolerance = 1e-10)
  expect_identical(via_alias$iter, via_named$iter)
  expect_identical(via_alias$mprod, via_named$mprod)
})

test_that("am_irlba block implementation delegates to am_block_lanczos", {
  x <- low_rank_matrix(n = 72L, p = 48L, rank = 4L)

  set.seed(1003L)
  expect_warning(
    via_irlba <- am_irlba(
      x,
      nv = 4L,
      nu = 3L,
      implementation = "block",
      block_size = 4L,
      n_steps = 6L,
      tol = 1e-4
    ),
    "ignores irlba-specific arguments"
  )
  set.seed(1003L)
  via_named <- am_block_lanczos(x, nv = 4L, nu = 3L, block_size = 4L, n_steps = 6L)

  expect_equal(via_irlba$d, via_named$d, tolerance = 1e-10)
  expect_equal(via_irlba$u, via_named$u, tolerance = 1e-10)
  expect_equal(via_irlba$v, via_named$v, tolerance = 1e-10)
  expect_identical(dim(via_irlba$u), c(nrow(x), 3L))
  expect_identical(dim(via_irlba$v), c(ncol(x), 4L))
  expect_identical(via_irlba$iter, via_named$iter)
  expect_identical(via_irlba$mprod, via_named$mprod)
})
