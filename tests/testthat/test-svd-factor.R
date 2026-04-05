make_dense <- function(n, p, seed = 42L) {
  set.seed(seed)
  X_host <- matrix(rnorm(n * p), nrow = n, ncol = p)
  list(host = X_host, am = adgeMatrix(X_host))
}

test_that("am_svd_factor singular values match base::svd", {
  dat <- make_dense(20L, 8L, seed = 1L)
  k <- 5L
  fac <- am_svd_factor(dat$am, k = k)
  base_svd <- base::svd(dat$host)

  expect_s4_class(fac, "amSVD")
  expect_equal(fac@d, base_svd$d[seq_len(k)], tolerance = 1e-10)
  expect_identical(fac@k, k)
  expect_identical(fac@source_id, dat$am@object_id)
  expect_equal(ncol(fac@u), k)
  expect_equal(ncol(fac@v), k)
  expect_equal(nrow(fac@u), 20L)
  expect_equal(nrow(fac@v), 8L)
})

test_that("am_svd_project matches t(U) %*% Y", {
  dat <- make_dense(15L, 6L, seed = 2L)
  k <- 4L
  fac <- am_svd_factor(dat$am, k = k)

  set.seed(99L)
  Y_single <- matrix(rnorm(15L), nrow = 15L, ncol = 1L)
  Y_multi <- matrix(rnorm(15L * 7L), nrow = 15L, ncol = 7L)

  Z_single <- am_svd_project(fac, Y_single)
  Z_multi <- am_svd_project(fac, Y_multi)

  expect_equal(dim(Z_single), c(k, 1L))
  expect_equal(dim(Z_multi), c(k, 7L))
  expect_equal(Z_single, crossprod(fac@u, Y_single), tolerance = 1e-12)
  expect_equal(Z_multi, crossprod(fac@u, Y_multi), tolerance = 1e-12)
})

test_that("am_svd_reconstruct matches V %*% diag(1/d) %*% Z", {
  dat <- make_dense(12L, 7L, seed = 3L)
  k <- 3L
  fac <- am_svd_factor(dat$am, k = k)

  set.seed(101L)
  Z <- matrix(rnorm(k * 5L), nrow = k, ncol = 5L)
  expected <- fac@v %*% diag(1 / fac@d) %*% Z

  expect_equal(am_svd_reconstruct(fac, Z), expected, tolerance = 1e-10)
})

test_that("am_pca_coef round-trip matches manual PCR formula", {
  dat <- make_dense(30L, 10L, seed = 4L)
  k <- 6L
  fac <- am_svd_factor(dat$am, k = k)

  set.seed(202L)
  Y <- matrix(rnorm(30L * 3L), nrow = 30L, ncol = 3L)

  # Manual PCR coefficient formula: V %*% diag(1/d) %*% t(U) %*% Y
  base_svd <- base::svd(dat$host, nu = k, nv = k)
  expected <- base_svd$v %*% diag(1 / base_svd$d[seq_len(k)]) %*% crossprod(base_svd$u, Y)

  expect_equal(am_pca_coef(fac, Y), expected, tolerance = 1e-8)
})

test_that("am_svd_factor cache reuse returns identical factor", {
  dat <- make_dense(18L, 9L, seed = 5L)
  k <- 4L
  fac1 <- am_svd_factor(dat$am, k = k)
  fac2 <- am_svd_factor(dat$am, k = k)

  # Identical same-object because of cache
  expect_identical(fac1, fac2)

  # Different k creates a different cache entry
  fac3 <- am_svd_factor(dat$am, k = 3L)
  expect_identical(fac3@k, 3L)
  expect_false(identical(fac1, fac3))
})

test_that("am_svd_project handles k=1, k=10, k=50 column Y matrices", {
  dat <- make_dense(40L, 12L, seed = 6L)
  k <- 5L
  fac <- am_svd_factor(dat$am, k = k)

  for (m in c(1L, 10L, 50L)) {
    set.seed(m)
    Y <- matrix(rnorm(40L * m), nrow = 40L, ncol = m)
    Z <- am_svd_project(fac, Y)
    expect_equal(dim(Z), c(k, m))
    expect_equal(Z, crossprod(fac@u, Y), tolerance = 1e-12)
  }
})
