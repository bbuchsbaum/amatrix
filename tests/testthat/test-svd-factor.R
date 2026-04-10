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

.opencl_svd_spec <- function() {
  specs <- optional_backend_specs()
  specs[[match("opencl", vapply(specs, `[[`, character(1), "backend"))]]
}

.register_optional_backend <- function(spec) {
  ns <- optional_backend_namespace(spec$package)
  get(spec$register_fun, envir = ns, inherits = FALSE)
}

test_that("svd_factor singular values match base::svd", {
  dat <- make_dense(20L, 8L, seed = 1L)
  k <- 5L
  fac <- svd_factor(dat$am, k = k)
  base_svd <- base::svd(dat$host)

  expect_s4_class(fac, "amSVD")
  expect_equal(fac@d, base_svd$d[seq_len(k)], tolerance = 1e-10)
  expect_identical(fac@k, k)
  expect_identical(fac@method, "exact")
  expect_identical(fac@engine, "exact_svd")
  expect_identical(fac@source_id, dat$am@object_id)
  expect_equal(ncol(fac@u), k)
  expect_equal(ncol(fac@v), k)
  expect_equal(nrow(fac@u), 20L)
  expect_equal(nrow(fac@v), 8L)
})

test_that("svd_factor subspace method recovers full-rank small problems", {
  dat <- make_dense(20L, 8L, seed = 11L)
  k <- 5L

  set.seed(1101L)
  fac <- svd_factor(dat$am, k = k, method = "subspace")
  base_svd <- base::svd(dat$host)

  expect_s4_class(fac, "amSVD")
  expect_identical(fac@method, "subspace")
  expect_identical(fac@engine, "gram")
  expect_equal(fac@d, base_svd$d[seq_len(k)], tolerance = 1e-8)
  expect_equal(base::crossprod(fac@u), diag(k), tolerance = 1e-8)
  expect_equal(base::crossprod(fac@v), diag(k), tolerance = 1e-8)
})

test_that("svd_factor subspace uses compiled resident plans for dense backends", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "svd_subspace_dense_plan",
    make_recording_backend(
      counter,
      supported_ops = character(),
      cold_supported_ops = character(),
      resident_supported_ops = c("matmul", "crossprod"),
      precision_modes = "fast"
    ),
    {
      dat <- make_dense(72L, 36L, seed = 13L)
      X <- adgeMatrix(
        dat$host,
        preferred_backend = "svd_subspace_dense_plan",
        precision = "fast"
      )

      set.seed(1301L)
      fac <- svd_factor(X, k = 6L, method = "subspace", n_oversamples = 6L, n_iter = 1L)
      ref <- base::svd(dat$host)$d[seq_len(6L)]
      rel_sv_err <- max(abs(fac@d - ref) / pmax(abs(ref), 1e-12))

      expect_identical(fac@backend, "svd_subspace_dense_plan")
      expect_lt(rel_sv_err, 0.05)
      expect_true(counter$resident_store >= 1L)
      expect_true(counter$matmul_resident >= 1L)
      expect_true(counter$crossprod_resident >= 1L)
      expect_true(counter$resident_drop >= 1L)
    }
  )
})

test_that("svd_factor subspace reuses sparse lhs resident plans", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "svd_subspace_sparse_plan",
    make_recording_backend(
      counter,
      supported_ops = character(),
      cold_supported_ops = character(),
      resident_supported_ops = character(),
      precision_modes = "fast",
      supports_sparse_ops = c("matmul", "crossprod"),
      supports_sparse_resident = TRUE
    ),
    {
      diag_vals <- seq(from = 48, to = 1, length.out = 48L)
      x_host <- Matrix::sparseMatrix(
        i = seq_along(diag_vals),
        j = seq_along(diag_vals),
        x = diag_vals,
        dims = c(64L, 48L)
      )
      X <- adgCMatrix(
        x_host,
        preferred_backend = "svd_subspace_sparse_plan",
        precision = "fast"
      )

      set.seed(1302L)
      fac <- svd_factor(X, k = 6L, method = "subspace", n_oversamples = 8L, n_iter = 2L)
      ref <- diag_vals[seq_len(6L)]
      rel_sv_err <- max(abs(fac@d - ref) / pmax(abs(ref), 1e-12))

      expect_identical(fac@backend, "svd_subspace_sparse_plan")
      expect_lt(rel_sv_err, 0.05)
      expect_true(counter$sparse_resident_store >= 1L)
      expect_true(counter$spmm_resident_key >= 1L)
      expect_true(counter$sparse_resident_drop >= 1L)
    }
  )
})

test_that("svd_factor subspace uses backend-native randomized SVD when available", {
  calls <- new.env(parent = emptyenv())
  calls$rsvd <- 0L

  register_mock_rsvd_backend(
    name = "mocksubspace",
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
  on.exit(drop_mock_backend("mocksubspace"), add = TRUE)

  dat <- make_dense(60L, 24L, seed = 12L)
  X <- adgeMatrix(dat$host, preferred_backend = "mocksubspace", precision = "fast")

  fac <- svd_factor(X, k = 8L, method = "subspace", n_oversamples = 6L, n_iter = 1L)

  expect_identical(calls$rsvd, 1L)
  expect_identical(fac@method, "subspace")
  expect_identical(fac@engine, "backend_rsvd")
  expect_identical(fac@backend, "mocksubspace")
  expect_equal(fac@d, base::svd(dat$host)$d[seq_len(8L)], tolerance = 1e-10)
})

test_that("svd_project matches t(U) %*% Y", {
  dat <- make_dense(15L, 6L, seed = 2L)
  k <- 4L
  fac <- svd_factor(dat$am, k = k)

  set.seed(99L)
  Y_single <- matrix(rnorm(15L), nrow = 15L, ncol = 1L)
  Y_multi <- matrix(rnorm(15L * 7L), nrow = 15L, ncol = 7L)

  Z_single <- svd_project(fac, Y_single)
  Z_multi <- svd_project(fac, Y_multi)

  expect_equal(dim(Z_single), c(k, 1L))
  expect_equal(dim(Z_multi), c(k, 7L))
  expect_equal(Z_single, base::crossprod(fac@u, Y_single), tolerance = 1e-12)
  expect_equal(Z_multi, base::crossprod(fac@u, Y_multi), tolerance = 1e-12)
})

test_that("svd_reconstruct matches V %*% diag(1/d) %*% Z", {
  dat <- make_dense(12L, 7L, seed = 3L)
  k <- 3L
  fac <- svd_factor(dat$am, k = k)

  set.seed(101L)
  Z <- matrix(rnorm(k * 5L), nrow = k, ncol = 5L)
  expected <- fac@v %*% diag(1 / fac@d) %*% Z

  expect_equal(svd_reconstruct(fac, Z), expected, tolerance = 1e-10)
})

test_that("pca_coef round-trip matches manual PCR formula", {
  dat <- make_dense(30L, 10L, seed = 4L)
  k <- 6L
  fac <- svd_factor(dat$am, k = k)

  set.seed(202L)
  Y <- matrix(rnorm(30L * 3L), nrow = 30L, ncol = 3L)

  # Manual PCR coefficient formula: V %*% diag(1/d) %*% t(U) %*% Y
  base_svd <- base::svd(dat$host, nu = k, nv = k)
  expected <- base_svd$v %*% base::diag(1 / base_svd$d[seq_len(k)]) %*% base::crossprod(base_svd$u, Y)

  expect_equal(pca_coef(fac, Y), expected, tolerance = 1e-8)
})

test_that("svd_factor cache reuse returns identical factor", {
  dat <- make_dense(18L, 9L, seed = 5L)
  k <- 4L
  fac1 <- svd_factor(dat$am, k = k)
  fac2 <- svd_factor(dat$am, k = k)

  # Identical same-object because of cache
  expect_identical(fac1, fac2)

  # Different k creates a different cache entry
  fac3 <- svd_factor(dat$am, k = 3L)
  expect_identical(fac3@k, 3L)
  expect_false(identical(fac1, fac3))
})

test_that("svd_project handles k=1, k=10, k=50 column Y matrices", {
  dat <- make_dense(40L, 12L, seed = 6L)
  k <- 5L
  fac <- svd_factor(dat$am, k = k)

  for (m in c(1L, 10L, 50L)) {
    set.seed(m)
    Y <- matrix(rnorm(40L * m), nrow = 40L, ncol = m)
    Z <- svd_project(fac, Y)
    expect_equal(dim(Z), c(k, m))
    expect_equal(Z, base::crossprod(fac@u, Y), tolerance = 1e-12)
  }
})

test_that("svd_factor auto keeps the exact path in strict precision", {
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
    rsvd = function(...) {
      stop("auto should not use rsvd for strict precision")
    },
    .package = "amatrix"
  )

  fac <- svd_factor(X, k = 5L)

  expect_identical(calls$exact, 1L)
  expect_identical(fac@method, "exact")
  expect_identical(fac@engine, "exact_svd")
  expect_identical(fac@backend, "cpu")
  expect_equal(fac@d, base::svd(dat$host)$d[seq_len(5L)], tolerance = 1e-10)
})

test_that("rsvd standalone keeps CPU fallback in strict precision", {
  register_mock_rsvd_backend(
    rsvd_fun = function(...) {
      stop("strict precision should not dispatch to backend rsvd")
    }
  )
  on.exit(drop_mock_backend(), add = TRUE)

  dat <- make_dense(48L, 24L, seed = 70L)
  X <- adgeMatrix(dat$host, preferred_backend = "mockgpu", precision = "strict")

  set.seed(701L)
  fit <- rsvd(X, k = 4L, n_oversamples = 5L, n_iter = 1L)
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

test_that("rsvd standalone uses backend rsvd in fast precision", {
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

  fit <- rsvd(X, k = 4L, n_oversamples = 5L, n_iter = 1L)

  expect_identical(calls$rsvd, 1L)
  expect_equal(fit$d, base::svd(dat$host)$d[seq_len(4L)], tolerance = 1e-10)
})

test_that("svd_factor auto uses rsvd for fast low-rank GPU factors", {
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
      stop("auto should choose rsvd for this factorization")
    },
    rsvd = function(x, k, n_oversamples = 10L, n_iter = 2L, ...) {
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

  fac <- svd_factor(X, k = 6L, n_oversamples = 7L, n_iter = 1L)

  expect_identical(calls$rsvd, 1L)
  expect_identical(calls$n_oversamples, 7L)
  expect_identical(calls$n_iter, 1L)
  expect_identical(fac@method, "rsvd")
  expect_identical(fac@engine, "backend_rsvd")
  expect_identical(fac@backend, "mockgpu")
  expect_equal(fac@d, base::svd(dat$host)$d[seq_len(6L)], tolerance = 1e-10)
})

test_that("svd_factor auto cutoff keeps small factors exact and medium factors rsvd", {
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

test_that("svd_factor auto chooses subspace for fast moderate-rank GPU factors", {
  register_mock_rsvd_backend()
  on.exit(drop_mock_backend(), add = TRUE)

  old_rsvd_min_dim <- getOption("amatrix.svd_factor.rsvd_min_dim")
  old_subspace_min_dim <- getOption("amatrix.svd_factor.subspace_min_dim")
  options(
    amatrix.svd_factor.rsvd_min_dim = 256L,
    amatrix.svd_factor.subspace_min_dim = 256L
  )
  on.exit(
    options(
      amatrix.svd_factor.rsvd_min_dim = old_rsvd_min_dim,
      amatrix.svd_factor.subspace_min_dim = old_subspace_min_dim
    ),
    add = TRUE
  )

  dat <- make_dense(320L, 280L, seed = 81L)
  X <- adgeMatrix(dat$host, preferred_backend = "mockgpu", precision = "fast")

  calls <- new.env(parent = emptyenv())
  calls$subspace <- 0L

  local_mocked_bindings(
    am_svd = function(...) {
      stop("auto should not choose exact for this factorization")
    },
    rsvd = function(...) {
      stop("auto should not choose rsvd for this factorization")
    },
    .amatrix_subspace_svd = function(X, k, n_oversamples = 10L, n_iter = 2L, target_backend = NULL, eps_rank = 1e-8) {
      calls$subspace <- calls$subspace + 1L
      exact <- base::svd(as.matrix(X), nu = k, nv = k)
      exact$rank_discovered <- k
      exact$core_solver <- "gram"
      exact$diag_history <- exact$d
      exact
    },
    .package = "amatrix"
  )

  fac <- svd_factor(X, k = 80L, n_oversamples = 6L, n_iter = 1L)

  expect_identical(calls$subspace, 1L)
  expect_identical(fac@method, "subspace")
  expect_identical(fac@engine, "gram")
  expect_identical(fac@backend, "mockgpu")
  expect_equal(fac@d, base::svd(dat$host)$d[seq_len(80L)], tolerance = 1e-10)
})

test_that("svd_factor auto chooses OpenCL subspace for fast moderate-rank factors", {
  spec <- .opencl_svd_spec()
  skip_if_backend_package_missing(spec)

  register_backend <- .register_optional_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    old_rsvd_min_dim <- getOption("amatrix.svd_factor.rsvd_min_dim")
    old_subspace_min_dim <- getOption("amatrix.svd_factor.subspace_min_dim")
    options(
      amatrix.svd_factor.rsvd_min_dim = 256L,
      amatrix.svd_factor.subspace_min_dim = 256L
    )
    on.exit(
      options(
        amatrix.svd_factor.rsvd_min_dim = old_rsvd_min_dim,
        amatrix.svd_factor.subspace_min_dim = old_subspace_min_dim
      ),
      add = TRUE
    )

    dat <- make_dense(420L, 300L, seed = 20260412L)
    X <- adgeMatrix(dat$host, preferred_backend = "opencl", precision = "fast")

    plan <- amatrix:::.amatrix_svd_factor_plan(
      X,
      k = 80L,
      method = "auto",
      n_oversamples = 6L,
      n_iter = 1L
    )
    expect_identical(plan$method, "subspace")
    expect_identical(plan$subspace_backend, "opencl")

    set.seed(20260412L)
    fac <- svd_factor(X, k = 80L, n_oversamples = 6L, n_iter = 1L)
    ref_d <- base::svd(dat$host, nu = 80L, nv = 80L)$d[seq_len(80L)]
    rel_sv <- abs(fac@d - ref_d) / pmax(abs(ref_d), 1e-12)

    expect_identical(fac@method, "subspace")
    expect_identical(fac@backend, "opencl")
    expect_identical(fac@precision, "fast")
    expect_true(fac@engine %in% c("gram", "qr", "svd_core"))
    expect_lt(max(rel_sv[seq_len(20L)]), 0.05)
  })
})

test_that("svd_factor sparse subspace keeps explicit OpenCL residency", {
  spec <- .opencl_svd_spec()
  skip_if_backend_package_missing(spec)
  skip_if_not_installed("Matrix")

  register_backend <- .register_optional_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    old_subspace_min_dim <- getOption("amatrix.svd_factor.subspace_min_dim")
    old_spmv_min_nnz <- getOption("amatrix.opencl.spmv_min_nnz")
    old_spmm_min_nnz <- getOption("amatrix.opencl.spmm_min_nnz")
    options(
      amatrix.svd_factor.subspace_min_dim = 128L,
      amatrix.opencl.spmv_min_nnz = 1L,
      amatrix.opencl.spmm_min_nnz = 1L
    )
    on.exit(
      options(
        amatrix.svd_factor.subspace_min_dim = old_subspace_min_dim,
        amatrix.opencl.spmv_min_nnz = old_spmv_min_nnz,
        amatrix.opencl.spmm_min_nnz = old_spmm_min_nnz
      ),
      add = TRUE
    )

    diag_vals <- seq(from = 256, to = 1, length.out = 256L)
    x_host <- Matrix::sparseMatrix(
      i = seq_along(diag_vals),
      j = seq_along(diag_vals),
      x = diag_vals,
      dims = c(384L, 256L)
    )
    X <- adgCMatrix(x_host, preferred_backend = "opencl", precision = "fast")

    plan <- amatrix:::.amatrix_svd_factor_plan(
      X,
      k = 32L,
      method = "subspace",
      n_oversamples = 8L,
      n_iter = 1L
    )
    expect_identical(plan$subspace_backend, "opencl")
    expect_identical(plan$factor_backend, "opencl")

    set.seed(20260409L)
    fac <- svd_factor(X, k = 32L, method = "subspace", n_oversamples = 8L, n_iter = 1L)
    ref <- diag_vals[seq_len(32L)]
    rel_sv <- abs(fac@d - ref) / pmax(abs(ref), 1e-12)

    expect_identical(fac@method, "subspace")
    expect_identical(fac@backend, "opencl")
    expect_identical(fac@precision, "fast")
    expect_true(fac@engine %in% c("gram", "qr", "svd_core"))
    expect_lt(max(rel_sv[seq_len(8L)]), 0.06)
  })
})

test_that("svd_factor exact supports OpenCL parity on dense inputs", {
  spec <- .opencl_svd_spec()
  skip_if_backend_package_missing(spec)

  register_backend <- .register_optional_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    dat <- make_dense(96L, 48L, seed = 20260413L)
    X <- adgeMatrix(dat$host, preferred_backend = "opencl", precision = "fast")

    plan <- amatrix:::.amatrix_svd_factor_plan(
      X,
      k = 12L,
      method = "exact",
      n_oversamples = 6L,
      n_iter = 1L
    )
    expect_identical(plan$method, "exact")

    fac <- svd_factor(X, k = 12L, method = "exact", n_oversamples = 6L, n_iter = 1L)
    ref <- base::svd(dat$host, nu = 12L, nv = 12L)

    expect_identical(fac@method, "exact")
    expect_identical(fac@backend, "opencl")
    expect_identical(fac@precision, "fast")
    expect_identical(fac@engine, "exact_svd")
    expect_equal(fac@d, ref$d[seq_len(12L)], tolerance = 1e-8)
    expect_equal(base::crossprod(fac@u), diag(12L), tolerance = 1e-8)
    expect_equal(base::crossprod(fac@v), diag(12L), tolerance = 1e-8)
  })
})

test_that("svd_factor cache key distinguishes exact, rsvd, and subspace factors", {
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
    rsvd = function(x, k, n_oversamples = 10L, n_iter = 2L, ...) {
      list(
        u = matrix(0, nrow = nrow(x), ncol = k),
        d = rep(7, k),
        v = matrix(0, nrow = ncol(x), ncol = k)
      )
    },
    .amatrix_subspace_svd = function(X, k, n_oversamples = 10L, n_iter = 2L, target_backend = NULL, eps_rank = 1e-8) {
      list(
        u = matrix(0, nrow = nrow(X), ncol = k),
        d = rep(5, k),
        v = matrix(0, nrow = ncol(X), ncol = k),
        rank_discovered = k,
        core_solver = "gram",
        diag_history = rep(1, k)
      )
    },
    .package = "amatrix"
  )

  fac_exact <- svd_factor(X, k = 4L, method = "exact")
  fac_rsvd <- svd_factor(X, k = 4L, method = "rsvd")
  fac_subspace <- svd_factor(X, k = 4L, method = "subspace")

  expect_equal(fac_exact@d, rep(11, 4L))
  expect_equal(fac_rsvd@d, rep(7, 4L))
  expect_equal(fac_subspace@d, rep(5, 4L))
  expect_identical(fac_exact@engine, "exact_svd")
  expect_identical(fac_rsvd@engine, "backend_rsvd")
  expect_identical(fac_subspace@engine, "gram")
  expect_false(identical(fac_exact, fac_rsvd))
  expect_false(identical(fac_exact, fac_subspace))
  expect_false(identical(fac_rsvd, fac_subspace))
})
