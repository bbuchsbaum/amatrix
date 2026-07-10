test_that("mlx rsvd bridge and factor path smoke test on native backend", {
  skip_if_not(amatrix_mlx_native_available(), "mlx native backend not available")

  old <- getOption("amatrix.mlx.available")
  options(amatrix.mlx.available = TRUE)
  on.exit(options(amatrix.mlx.available = old), add = TRUE)

  set.seed(20260405L)
  x <- matrix(rnorm(192L * 128L), nrow = 192L, ncol = 128L)
  k <- 8L
  n_oversamples <- 6L
  n_iter <- 1L

  bridge_fit <- amatrix_mlx_rsvd(
    x,
    k = k,
    n_oversamples = n_oversamples,
    n_iter = n_iter
  )

  factor_fit <- amatrix::svd_factor(
    amatrix::adgeMatrix(x, preferred_backend = "mlx", precision = "fast"),
    k = k,
    method = "rsvd",
    n_oversamples = n_oversamples,
    n_iter = n_iter
  )

  ref <- base::svd(x, nu = k, nv = k)
  ref_d <- ref$d[seq_len(k)]
  bridge_rel_sv_err <- max(abs(bridge_fit$d - ref_d) / pmax(abs(ref_d), 1e-12))
  factor_rel_sv_err <- max(abs(factor_fit@d - ref_d) / pmax(abs(ref_d), 1e-12))

  # The bridge returns diagnostic fields (iter, mprod) beyond the core u/d/v
  # contract; require the contract fields, allow diagnostics.
  expect_contains(names(bridge_fit), c("u", "d", "v"))
  expect_equal(dim(bridge_fit$u), c(nrow(x), k))
  expect_length(bridge_fit$d, k)
  expect_equal(dim(bridge_fit$v), c(ncol(x), k))
  expect_true(all(is.finite(bridge_fit$u)))
  expect_true(all(is.finite(bridge_fit$d)))
  expect_true(all(is.finite(bridge_fit$v)))
  expect_true(all(diff(bridge_fit$d) <= 1e-8))
  expect_lte(bridge_rel_sv_err, 0.15)

  expect_s4_class(factor_fit, "amSVD")
  expect_identical(factor_fit@backend, "mlx")
  expect_equal(dim(factor_fit@u), c(nrow(x), k))
  expect_length(factor_fit@d, k)
  expect_equal(dim(factor_fit@v), c(ncol(x), k))
  expect_true(all(is.finite(factor_fit@u)))
  expect_true(all(is.finite(factor_fit@d)))
  expect_true(all(is.finite(factor_fit@v)))
  expect_lte(factor_rel_sv_err, 0.15)
})
