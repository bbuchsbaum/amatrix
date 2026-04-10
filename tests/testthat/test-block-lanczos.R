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

.opencl_block_spec <- function() {
  specs <- optional_backend_specs()
  specs[[match("opencl", vapply(specs, `[[`, character(1), "backend"))]]
}

.register_optional_backend <- function(spec) {
  ns <- optional_backend_namespace(spec$package)
  get(spec$register_fun, envir = ns, inherits = FALSE)
}

test_that("block_lanczos recovers a low-rank matrix accurately", {
  x <- low_rank_matrix()
  ref <- La.svd(x, nu = 5L, nv = 5L)

  set.seed(1001L)
  res <- block_lanczos(x, nv = 5L, nu = 5L, block_size = 4L, n_steps = 6L)

  expect_identical(dim(res$u), c(nrow(x), 5L))
  expect_identical(dim(res$v), c(ncol(x), 5L))
  expect_equal(res$d, ref$d[seq_len(5L)], tolerance = 1e-6)

  rel_err <- norm(svd_reconstruct(res) - x, type = "F") / max(norm(x, type = "F"), 1e-12)
  expect_lt(rel_err, 1e-6)
})

test_that("block_svd remains an alias for block_lanczos", {
  x <- low_rank_matrix(n = 60L, p = 40L, rank = 4L)

  set.seed(1002L)
  via_alias <- block_svd(x, k = 4L, block_size = 4L, n_steps = 6L)
  set.seed(1002L)
  via_named <- block_lanczos(x, nv = 4L, nu = 4L, block_size = 4L, n_steps = 6L)

  expect_equal(via_alias$d, via_named$d, tolerance = 1e-10)
  expect_equal(svd_reconstruct(via_alias), svd_reconstruct(via_named), tolerance = 1e-10)
  expect_identical(via_alias$iter, via_named$iter)
  expect_identical(via_alias$mprod, via_named$mprod)
})

test_that("block_lanczos supports oversampled block sizes above k", {
  x <- low_rank_matrix(n = 84L, p = 56L, rank = 6L)
  ref <- La.svd(x, nu = 4L, nv = 4L)

  set.seed(1005L)
  fit <- block_lanczos(x, nv = 4L, nu = 4L, block_size = 6L, n_steps = 5L)

  expect_identical(dim(fit$u), c(nrow(x), 4L))
  expect_identical(dim(fit$v), c(ncol(x), 4L))
  expect_equal(fit$d, ref$d[seq_len(4L)], tolerance = 1e-6)
})

test_that("irlba block implementation delegates to block_lanczos", {
  x <- low_rank_matrix(n = 72L, p = 48L, rank = 4L)

  set.seed(1003L)
  expect_warning(
    via_irlba <- irlba(
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
  via_named <- block_lanczos(x, nv = 4L, nu = 3L, block_size = 4L, n_steps = 6L)

  expect_equal(via_irlba$d, via_named$d, tolerance = 1e-10)
  expect_equal(via_irlba$u, via_named$u, tolerance = 1e-10)
  expect_equal(via_irlba$v, via_named$v, tolerance = 1e-10)
  expect_identical(dim(via_irlba$u), c(nrow(x), 3L))
  expect_identical(dim(via_irlba$v), c(ncol(x), 4L))
  expect_identical(via_irlba$iter, via_named$iter)
  expect_identical(via_irlba$mprod, via_named$mprod)
})

test_that("block_lanczos uses compiled resident plans for dense non-MLX backends", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "block_lanczos_dense_plan",
    make_recording_backend(
      counter,
      supported_ops = character(),
      cold_supported_ops = character(),
      resident_supported_ops = c("matmul", "crossprod")
    ),
    {
      x_host <- low_rank_matrix(n = 72L, p = 48L, rank = 4L)
      x_plan <- adgeMatrix(
        x_host,
        preferred_backend = "block_lanczos_dense_plan",
        precision = "strict"
      )

      set.seed(1007L)
      fit <- block_lanczos(x_plan, nv = 4L, nu = 4L, block_size = 4L, n_steps = 6L)

      ref <- La.svd(x_host, nu = 4L, nv = 4L)$d[seq_len(4L)]
      expect_equal(fit$d, ref, tolerance = 1e-6)
      expect_true(counter$resident_store >= 1L)
      expect_true(counter$matmul_resident >= 1L)
      expect_true(counter$crossprod_resident >= 1L)
      expect_true(counter$resident_drop >= 1L)
    }
  )
})

test_that("block_lanczos reuses a sparse lhs through compiled resident plans", {
  skip_if_not_installed("Matrix")

  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "block_lanczos_sparse_plan",
    make_recording_backend(
      counter,
      supported_ops = character(),
      cold_supported_ops = character(),
      resident_supported_ops = character(),
      supports_sparse_ops = c("matmul", "crossprod"),
      supports_sparse_resident = TRUE
    ),
    {
      diag_vals <- seq(from = 32, to = 1, length.out = 32L)
      x_host <- Matrix::sparseMatrix(
        i = seq_along(diag_vals),
        j = seq_along(diag_vals),
        x = diag_vals,
        dims = c(48L, 32L)
      )
      x_plan <- adgCMatrix(
        x_host,
        preferred_backend = "block_lanczos_sparse_plan",
        precision = "strict"
      )

      set.seed(1008L)
      fit <- block_lanczos(x_plan, nv = 4L, nu = 4L, block_size = 4L, n_steps = 5L)

      rel_sv_err <- max(abs(fit$d[seq_len(4L)] - diag_vals[seq_len(4L)]) / diag_vals[seq_len(4L)])
      expect_lt(rel_sv_err, 0.01)
      expect_true(counter$sparse_resident_store >= 1L)
      expect_true(counter$spmm_resident >= 1L)
      expect_true(counter$sparse_resident_drop >= 1L)
    }
  )
})

test_that("block_lanczos runs on MLX fast matrices", {
  skip_if_not_installed("amatrix.mlx")
  skip_if_not(
    isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE)),
    "mlx backend not available"
  )

  old <- options(amatrix.mlx.available = TRUE)
  on.exit(options(old), add = TRUE)
  amatrix.mlx::amatrix_mlx_register(overwrite = TRUE)

  x_host <- low_rank_matrix(n = 96L, p = 64L, rank = 6L)
  x_mlx <- adgeMatrix(x_host, preferred_backend = "mlx", precision = "fast")

  set.seed(1004L)
  fit <- block_lanczos(x_mlx, nv = 6L, nu = 6L, block_size = 8L, n_steps = 5L)

  expect_identical(dim(fit$u), c(nrow(x_host), 6L))
  expect_identical(dim(fit$v), c(ncol(x_host), 6L))
  expect_true(all(is.finite(fit$d)))

  ref <- La.svd(x_host, nu = 6L, nv = 6L)$d[seq_len(6L)]
  rel_sv_err <- max(abs(fit$d - ref) / pmax(abs(ref), 1e-12))
  expect_lt(rel_sv_err, 0.1)
})

test_that("block_lanczos runs on OpenCL fast matrices", {
  spec <- .opencl_block_spec()
  skip_if_backend_package_missing(spec)

  register_backend <- .register_optional_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    x_host <- low_rank_matrix(n = 96L, p = 64L, rank = 6L)
    x_opencl <- adgeMatrix(x_host, preferred_backend = "opencl", precision = "fast")

    set.seed(1006L)
    fit <- block_lanczos(x_opencl, nv = 6L, nu = 6L, block_size = 8L, n_steps = 5L)

    expect_identical(dim(fit$u), c(nrow(x_host), 6L))
    expect_identical(dim(fit$v), c(ncol(x_host), 6L))
    expect_true(all(is.finite(fit$d)))

    ref <- La.svd(x_host, nu = 6L, nv = 6L)$d[seq_len(6L)]
    rel_sv_err <- max(abs(fit$d - ref) / pmax(abs(ref), 1e-12))
    expect_lt(rel_sv_err, 1e-5)
  })
})
