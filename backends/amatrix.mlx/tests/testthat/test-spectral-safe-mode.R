test_that("mlx spectral safe mode falls back without native bridge calls", {
  old_opts <- options(amatrix.mlx.available = TRUE, amatrix.mlx.safe_spectral = TRUE)
  on.exit(options(old_opts), add = TRUE)

  set.seed(20260409L)
  x <- matrix(rnorm(96L * 64L), nrow = 96L, ncol = 64L)

  sv_fun <- getFromNamespace("amatrix_mlx_svd", "amatrix.mlx")
  rsvd_fun <- getFromNamespace("amatrix_mlx_rsvd", "amatrix.mlx")
  sv_fit <- sv_fun(x, nu = 8L, nv = 8L)
  rs_fit <- rsvd_fun(x, k = 8L, n_oversamples = 4L, n_iter = 1L)

  expect_named(sv_fit, c("d", "u", "v"))
  expect_equal(dim(sv_fit$u), c(nrow(x), 8L))
  expect_equal(dim(sv_fit$v), c(ncol(x), 8L))
  expect_length(sv_fit$d, min(nrow(x), ncol(x)))

  expect_named(rs_fit, c("u", "d", "v", "iter", "mprod"))
  expect_equal(dim(rs_fit$u), c(nrow(x), 8L))
  expect_equal(dim(rs_fit$v), c(ncol(x), 8L))
  expect_length(rs_fit$d, 8L)
  expect_true(all(is.finite(rs_fit$d)))
})
