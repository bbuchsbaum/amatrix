make_dense <- function(n, p, seed = 42L) {
  set.seed(seed)
  X_host <- matrix(rnorm(n * p), nrow = n, ncol = p)
  list(host = X_host, am = adgeMatrix(X_host))
}

register_mock_rsvd_backend <- function(name = "mockgpu", rsvd_fun = NULL) {
  if (exists(name, envir = amatrix:::.amatrix_state$backends, inherits = FALSE)) {
    stop(sprintf("backend '%s' already exists", name), call. = FALSE)
  }

  backend <- amatrix:::.amatrix_cpu_backend()
  cpu_capabilities <- backend$capabilities()
  cpu_features <- backend$features()

  backend$capabilities <- function() unique(c(cpu_capabilities, "rsvd"))
  backend$features <- function() unique(c(cpu_features, "rsvd"))
  backend$precision_modes <- function() "fast"
  backend$available <- function() TRUE
  backend$rsvd <- if (is.null(rsvd_fun)) {
    function(x, k, n_oversamples = 10L, n_iter = 4L) {
      exact <- base::svd(as.matrix(x), nu = k, nv = k)
      list(
        u = exact$u[, seq_len(k), drop = FALSE],
        d = exact$d[seq_len(k)],
        v = exact$v[, seq_len(k), drop = FALSE]
      )
    }
  } else {
    rsvd_fun
  }

  amatrix_register_backend(name, backend, overwrite = TRUE)
  invisible(name)
}

drop_mock_backend <- function(name = "mockgpu") {
  if (exists(name, envir = amatrix:::.amatrix_state$backends, inherits = FALSE)) {
    rm(list = name, envir = amatrix:::.amatrix_state$backends)
  }
  invisible(name)
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

test_that("am_svd_factor auto keeps the exact path in strict precision", {
  register_mock_rsvd_backend()
  on.exit(drop_mock_backend(), add = TRUE)

  dat <- make_dense(48L, 24L, seed = 7L)
  X <- adgeMatrix(dat$host, preferred_backend = "mockgpu", precision = "strict")

  calls <- new.env(parent = emptyenv())
  calls$exact <- 0L

  local_mocked_bindings(
    am_svd = function(x, nu = min(dim(x)), nv = min(dim(x)), LINPACK = FALSE, ...) {
      calls$exact <- calls$exact + 1L
      base::svd(as.matrix(x), nu = nu, nv = nv, ...)
    },
    am_rsvd = function(...) {
      stop("auto should not use am_rsvd for strict precision")
    },
    .package = "amatrix"
  )

  fac <- am_svd_factor(X, k = 5L)

  expect_identical(calls$exact, 1L)
  expect_equal(fac@d, base::svd(dat$host)$d[seq_len(5L)], tolerance = 1e-10)
})

test_that("am_rsvd standalone keeps CPU fallback in strict precision", {
  register_mock_rsvd_backend(
    rsvd_fun = function(...) {
      stop("strict precision should not dispatch to backend rsvd")
    }
  )
  on.exit(drop_mock_backend(), add = TRUE)

  dat <- make_dense(48L, 24L, seed = 70L)
  X <- adgeMatrix(dat$host, preferred_backend = "mockgpu", precision = "strict")

  set.seed(701L)
  fit <- am_rsvd(X, k = 4L, n_oversamples = 5L, n_iter = 1L)
  set.seed(701L)
  ref <- if (requireNamespace("irlba", quietly = TRUE)) {
    irlba::svdr(dat$host, k = 4L, extra = 5L, it = 1L)
  } else {
    base::svd(dat$host, nu = 4L, nv = 4L)
  }

  expect_equal(fit$d, ref$d[seq_len(4L)], tolerance = 1e-8)
  expect_equal(dim(fit$u), c(48L, 4L))
  expect_equal(dim(fit$v), c(24L, 4L))
})

test_that("am_rsvd standalone uses backend rsvd in fast precision", {
  calls <- new.env(parent = emptyenv())
  calls$rsvd <- 0L

  register_mock_rsvd_backend(
    rsvd_fun = function(x, k, n_oversamples = 10L, n_iter = 4L) {
      calls$rsvd <- calls$rsvd + 1L
      exact <- base::svd(as.matrix(x), nu = k, nv = k)
      list(
        u = exact$u[, seq_len(k), drop = FALSE],
        d = exact$d[seq_len(k)],
        v = exact$v[, seq_len(k), drop = FALSE]
      )
    }
  )
  on.exit(drop_mock_backend(), add = TRUE)

  dat <- make_dense(48L, 24L, seed = 71L)
  X <- adgeMatrix(dat$host, preferred_backend = "mockgpu", precision = "fast")

  fit <- am_rsvd(X, k = 4L, n_oversamples = 5L, n_iter = 1L)

  expect_identical(calls$rsvd, 1L)
  expect_equal(fit$d, base::svd(dat$host)$d[seq_len(4L)], tolerance = 1e-10)
})

test_that("am_svd_factor auto uses rsvd for fast low-rank GPU factors", {
  register_mock_rsvd_backend()
  on.exit(drop_mock_backend(), add = TRUE)
  old_threshold <- getOption("amatrix.svd_factor.rsvd_min_dim")
  options(amatrix.svd_factor.rsvd_min_dim = 256L)
  on.exit(options(amatrix.svd_factor.rsvd_min_dim = old_threshold), add = TRUE)

  dat <- make_dense(320L, 280L, seed = 8L)
  X <- adgeMatrix(dat$host, preferred_backend = "mockgpu", precision = "fast")

  calls <- new.env(parent = emptyenv())
  calls$rsvd <- 0L
  calls$n_oversamples <- NA_integer_
  calls$n_iter <- NA_integer_

  local_mocked_bindings(
    am_svd = function(...) {
      stop("auto should choose am_rsvd for this factorization")
    },
    am_rsvd = function(x, k, n_oversamples = 10L, n_iter = 2L, ...) {
      calls$rsvd <- calls$rsvd + 1L
      calls$n_oversamples <- as.integer(n_oversamples)
      calls$n_iter <- as.integer(n_iter)
      exact <- base::svd(as.matrix(x), nu = k, nv = k)
      list(
        u = exact$u[, seq_len(k), drop = FALSE],
        d = exact$d[seq_len(k)],
        v = exact$v[, seq_len(k), drop = FALSE]
      )
    },
    .package = "amatrix"
  )

  fac <- am_svd_factor(X, k = 6L, n_oversamples = 7L, n_iter = 1L)

  expect_identical(calls$rsvd, 1L)
  expect_identical(calls$n_oversamples, 7L)
  expect_identical(calls$n_iter, 1L)
  expect_identical(fac@backend, "mockgpu")
  expect_equal(fac@d, base::svd(dat$host)$d[seq_len(6L)], tolerance = 1e-10)
})

test_that("am_svd_factor auto cutoff keeps small factors exact and medium factors rsvd", {
  register_mock_rsvd_backend()
  on.exit(drop_mock_backend(), add = TRUE)

  plan_small <- amatrix:::.amatrix_svd_factor_plan(
    adgeMatrix(matrix(rnorm(400L * 320L), nrow = 400L, ncol = 320L),
               preferred_backend = "mockgpu",
               precision = "fast"),
    k = 20L,
    method = "auto",
    n_oversamples = 10L,
    n_iter = 2L
  )

  plan_medium <- amatrix:::.amatrix_svd_factor_plan(
    adgeMatrix(matrix(rnorm(500L * 400L), nrow = 500L, ncol = 400L),
               preferred_backend = "mockgpu",
               precision = "fast"),
    k = 20L,
    method = "auto",
    n_oversamples = 10L,
    n_iter = 2L
  )

  expect_identical(plan_small$method, "exact")
  expect_identical(plan_medium$method, "rsvd")
})

test_that("am_svd_factor cache key distinguishes exact and rsvd factors", {
  register_mock_rsvd_backend()
  on.exit(drop_mock_backend(), add = TRUE)

  dat <- make_dense(320L, 120L, seed = 9L)
  X <- adgeMatrix(dat$host, preferred_backend = "mockgpu", precision = "fast")

  local_mocked_bindings(
    am_svd = function(x, nu = min(dim(x)), nv = min(dim(x)), LINPACK = FALSE, ...) {
      list(
        u = matrix(0, nrow = nrow(x), ncol = nu),
        d = rep(11, nu),
        v = matrix(0, nrow = ncol(x), ncol = nv)
      )
    },
    am_rsvd = function(x, k, n_oversamples = 10L, n_iter = 2L, ...) {
      list(
        u = matrix(0, nrow = nrow(x), ncol = k),
        d = rep(7, k),
        v = matrix(0, nrow = ncol(x), ncol = k)
      )
    },
    .package = "amatrix"
  )

  fac_exact <- am_svd_factor(X, k = 4L, method = "exact")
  fac_rsvd <- am_svd_factor(X, k = 4L, method = "rsvd")

  expect_equal(fac_exact@d, rep(11, 4L))
  expect_equal(fac_rsvd@d, rep(7, 4L))
  expect_false(identical(fac_exact, fac_rsvd))
})
