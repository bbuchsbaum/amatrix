test_that("metal backend registers with amatrix", {
  amatrix_metal_register(overwrite = TRUE)
  expect_true("metal" %in% amatrix::amatrix_backend_names())
})

test_that("metal backend advertises sparse-only capabilities", {
  backend <- amatrix_metal_backend()
  old <- options(amatrix.metal.available = TRUE, amatrix.metal.spmm_min_nnz = 1L)
  on.exit(options(old), add = TRUE)

  sparse <- amatrix::adgCMatrix(
    Matrix::rsparsematrix(64, 64, density = 0.05),
    preferred_backend = "metal",
    precision = "fast"
  )
  dense <- amatrix::adgeMatrix(matrix(1, nrow = 64, ncol = 8), precision = "fast")

  expect_identical(backend$available(), amatrix_metal_is_available())
  expect_false(backend$supports("matmul", dense, y = matrix(1, nrow = 8, ncol = 8)))
  expect_true(backend$supports("matmul", sparse, y = matrix(1, nrow = 64, ncol = 8)))
  expect_true(backend$supports("crossprod", sparse, y = matrix(1, nrow = 64, ncol = 4)))
  expect_true(backend$supports("tcrossprod", sparse, y = matrix(1, nrow = 4, ncol = 64)))
})

test_that("metal resident support gates narrow SpMV but keeps SpMM hot paths", {
  backend <- amatrix_metal_backend()
  old <- options(amatrix.metal.available = TRUE)
  on.exit(options(old), add = TRUE)

  sparse <- amatrix::adgCMatrix(
    Matrix::rsparsematrix(512, 128, density = 0.05),
    preferred_backend = "metal",
    precision = "fast"
  )

  expect_false(backend$supports_resident("matmul", sparse, y = matrix(1, nrow = 128, ncol = 1)))
  expect_true(backend$supports_resident("matmul", sparse, y = matrix(1, nrow = 128, ncol = 4)))
  expect_false(backend$supports_resident("tcrossprod", sparse, y = matrix(1, nrow = 1, ncol = 128)))
  expect_true(backend$supports_resident("tcrossprod", sparse, y = matrix(1, nrow = 4, ncol = 128)))
})

test_that("metal bridge info reports coherent status", {
  info <- amatrix_metal_bridge_info()

  expect_type(info$compiled, "logical")
  expect_identical(info$native, amatrix_metal_native_available())
  expect_identical(info$available, amatrix_metal_is_available())
  expect_true(info$engine %in% c("mock-metal-bridge", "metal-unavailable", "metal-runtime"))
  expect_identical(info$capabilities, amatrix_metal_capabilities())
})

test_that("metal bridge exposes opt-in profile counters", {
  expect_true(amatrix_metal_profile_enable(TRUE, reset = TRUE))
  on.exit(amatrix_metal_profile_enable(FALSE), add = TRUE)

  profile <- amatrix_metal_profile()
  expect_named(profile, c(
    "enabled",
    "sparse_upload_ms",
    "sparse_upload_count",
    "sparse_upload_reuse_count",
    "dense_upload_ms",
    "dense_upload_count",
    "spmm_submit_ms",
    "spmm_submit_count",
    "spmm_wait_ms",
    "spmm_wait_count",
    "dense_sparse_submit_ms",
    "dense_sparse_submit_count",
    "dense_sparse_wait_ms",
    "dense_sparse_wait_count",
    "transpose_submit_ms",
    "transpose_submit_count",
    "transpose_wait_ms",
    "transpose_wait_count",
    "pending_wait_ms",
    "pending_wait_count",
    "materialize_ms",
    "materialize_count",
    "sparse_resident_count",
    "dense_resident_count",
    "profile_schema_version"
  ))
  expect_equal(unname(profile[["enabled"]]), 1)

  if (isTRUE(amatrix_metal_native_available(force = TRUE))) {
    set.seed(7)
    s_host <- Matrix::rsparsematrix(64, 32, density = 0.08)
    b_host <- matrix(rnorm(32 * 4), nrow = 32, ncol = 4)
    sp_key <- paste0("metal-sp-", sample.int(1e6, 1))
    y_key <- paste0("metal-y-", sample.int(1e6, 1))
    out_key <- paste0("metal-out-", sample.int(1e6, 1))

    on.exit(amatrix_metal_sparse_resident_drop(sp_key), add = TRUE)
    on.exit(amatrix_metal_resident_drop(y_key), add = TRUE)
    on.exit(amatrix_metal_resident_drop(out_key), add = TRUE)

    amatrix_metal_sparse_resident_store(sp_key, s_host)
    amatrix_metal_resident_store(y_key, b_host)
    amatrix_metal_backend()$spmm_resident_key(sp_key, y_key, out_key, trans_lhs = FALSE)

    profile <- amatrix_metal_profile()
    expect_gt(unname(profile[["sparse_upload_count"]]), 0)
    expect_gt(unname(profile[["dense_upload_count"]]), 0)
    expect_gt(unname(profile[["spmm_submit_count"]]), 0)
    expect_gt(unname(profile[["materialize_count"]]), 0)
  }

  amatrix_metal_profile_reset()
  profile <- amatrix_metal_profile()
  expect_equal(unname(profile[["spmm_submit_count"]]), 0)
})

test_that("metal sparse bridge matches Matrix on small products when native is available", {
  skip_if_not(isTRUE(amatrix_metal_native_available(force = TRUE)))

  set.seed(1)
  x_host <- Matrix::rsparsematrix(256, 128, density = 0.05)
  y_host <- matrix(rnorm(128 * 16), nrow = 128, ncol = 16)
  y_cross_host <- matrix(rnorm(256 * 9), nrow = 256, ncol = 9)
  ref <- as.matrix(x_host %*% y_host)

  direct <- .Call(
    "amatrix_metal_spmm_bridge",
    as.double(x_host@x),
    as.integer(x_host@p),
    as.integer(x_host@i),
    as.integer(x_host@Dim),
    y_host,
    FALSE,
    PACKAGE = "amatrix.metal"
  )

  key <- paste0("metal-sp-", sample.int(1e6, 1))
  invisible(.Call(
    "amatrix_metal_sparse_store_bridge",
    key,
    as.double(x_host@x),
    as.integer(x_host@p),
    as.integer(x_host@i),
    as.integer(x_host@Dim),
    PACKAGE = "amatrix.metal"
  ))
  resident <- .Call("amatrix_metal_spmm_resident_bridge", key, y_host, FALSE, PACKAGE = "amatrix.metal")
  direct_cross <- .Call(
    "amatrix_metal_spmm_bridge",
    as.double(x_host@x),
    as.integer(x_host@p),
    as.integer(x_host@i),
    as.integer(x_host@Dim),
    y_cross_host,
    TRUE,
    PACKAGE = "amatrix.metal"
  )
  invisible(.Call("amatrix_metal_sparse_drop_bridge", key, PACKAGE = "amatrix.metal"))

  expect_equal(direct, ref, tolerance = 1e-5)
  expect_equal(resident, ref, tolerance = 1e-5)
  expect_equal(direct_cross, as.matrix(Matrix::crossprod(x_host, y_cross_host)), tolerance = 1e-5)
})

test_that("metal backend can reuse sparse and dense resident keys", {
  set.seed(2)
  x_host <- Matrix::rsparsematrix(128, 64, density = 0.05)
  y_host <- matrix(rnorm(64 * 10), nrow = 64, ncol = 10)
  ref <- as.matrix(x_host %*% y_host)

  sp_key <- paste0("metal-sp-", sample.int(1e6, 1))
  y_key <- paste0("metal-y-", sample.int(1e6, 1))
  out_key <- paste0("metal-out-", sample.int(1e6, 1))

  on.exit(amatrix_metal_sparse_resident_drop(sp_key), add = TRUE)
  on.exit(amatrix_metal_resident_drop(y_key), add = TRUE)
  on.exit(amatrix_metal_resident_drop(out_key), add = TRUE)

  amatrix_metal_sparse_resident_store(sp_key, x_host)
  amatrix_metal_resident_store(y_key, y_host)

  value <- amatrix_metal_backend()$spmm_resident_key(sp_key, y_key, out_key, trans_lhs = FALSE)

  expect_true(amatrix_metal_resident_has(out_key))
  expect_equal(value, ref, tolerance = 1e-5)
  expect_equal(amatrix_metal_resident_materialize(out_key), ref, tolerance = 1e-5)
})

test_that("metal backend can reuse resident keys for sparse crossprod", {
  set.seed(3)
  x_host <- Matrix::rsparsematrix(96, 48, density = 0.06)
  y_host <- matrix(rnorm(96 * 7), nrow = 96, ncol = 7)
  ref <- as.matrix(Matrix::crossprod(x_host, y_host))

  sp_key <- paste0("metal-sp-", sample.int(1e6, 1))
  y_key <- paste0("metal-y-", sample.int(1e6, 1))
  out_key <- paste0("metal-out-", sample.int(1e6, 1))

  on.exit(amatrix_metal_sparse_resident_drop(sp_key), add = TRUE)
  on.exit(amatrix_metal_resident_drop(y_key), add = TRUE)
  on.exit(amatrix_metal_resident_drop(out_key), add = TRUE)

  amatrix_metal_sparse_resident_store(sp_key, x_host)
  amatrix_metal_resident_store(y_key, y_host)

  value <- amatrix_metal_backend()$spmm_resident_key(sp_key, y_key, out_key, trans_lhs = TRUE)

  expect_true(amatrix_metal_resident_has(out_key))
  expect_equal(value, ref, tolerance = 1e-5)
  expect_equal(amatrix_metal_resident_materialize(out_key), ref, tolerance = 1e-5)
})

test_that("metal backend can multiply dense resident lhs by sparse resident rhs", {
  set.seed(5)
  a_host <- matrix(rnorm(6 * 48), nrow = 6, ncol = 48)
  s_host <- Matrix::rsparsematrix(48, 24, density = 0.06)
  ref <- as.matrix(a_host %*% s_host)

  x_key <- paste0("metal-x-", sample.int(1e6, 1))
  sp_key <- paste0("metal-sp-", sample.int(1e6, 1))
  out_key <- paste0("metal-out-", sample.int(1e6, 1))

  on.exit(amatrix_metal_resident_drop(x_key), add = TRUE)
  on.exit(amatrix_metal_sparse_resident_drop(sp_key), add = TRUE)
  on.exit(amatrix_metal_resident_drop(out_key), add = TRUE)

  amatrix_metal_resident_store(x_key, a_host)
  amatrix_metal_sparse_resident_store(sp_key, s_host)

  value <- amatrix_metal_backend()$dense_sparse_matmul_resident_key(x_key, sp_key, out_key, defer = FALSE)

  expect_true(amatrix_metal_resident_has(out_key))
  expect_equal(value, ref, tolerance = 1e-5)
  expect_equal(amatrix_metal_resident_materialize(out_key), ref, tolerance = 1e-5)
})

test_that("metal backend can transpose a resident dense matrix", {
  x_host <- matrix(seq_len(24), nrow = 6, ncol = 4)
  x_key <- paste0("metal-x-", sample.int(1e6, 1))
  out_key <- paste0("metal-t-", sample.int(1e6, 1))

  on.exit(amatrix_metal_resident_drop(x_key), add = TRUE)
  on.exit(amatrix_metal_resident_drop(out_key), add = TRUE)

  amatrix_metal_resident_store(x_key, x_host)
  amatrix_metal_backend()$transpose_resident(x_key, out_key)

  expect_true(amatrix_metal_resident_has(out_key))
  expect_equal(amatrix_metal_resident_materialize(out_key), t(x_host), tolerance = 0)
})

test_that("metal sparse resident key can defer host materialization", {
  skip_if_not(isTRUE(amatrix_metal_native_available(force = TRUE)))

  old <- options(
    amatrix.defer_host = TRUE
  )
  on.exit(options(old), add = TRUE)

  set.seed(4)
  s_host <- Matrix::rsparsematrix(128, 64, density = 0.05)
  b_host <- matrix(rnorm(64 * 8), nrow = 64, ncol = 8)
  ref <- as.matrix(s_host %*% b_host)

  sp_key <- paste0("metal-sp-", sample.int(1e6, 1))
  y_key <- paste0("metal-y-", sample.int(1e6, 1))
  out_key <- paste0("metal-out-", sample.int(1e6, 1))

  on.exit(amatrix_metal_sparse_resident_drop(sp_key), add = TRUE)
  on.exit(amatrix_metal_resident_drop(y_key), add = TRUE)
  on.exit(amatrix_metal_resident_drop(out_key), add = TRUE)

  amatrix_metal_sparse_resident_store(sp_key, s_host)
  amatrix_metal_resident_store(y_key, b_host)

  value <- amatrix_metal_backend()$spmm_resident_key(sp_key, y_key, out_key, trans_lhs = FALSE, defer = TRUE)

  expect_null(value)
  expect_true(amatrix_metal_resident_has(out_key))
  expect_equal(amatrix_metal_resident_materialize(out_key), ref, tolerance = 1e-5)
})

test_that("metal deferred resident sparse product can chain into resident transpose", {
  skip_if_not(isTRUE(amatrix_metal_native_available(force = TRUE)))

  set.seed(6)
  s_host <- Matrix::rsparsematrix(96, 48, density = 0.07)
  b_host <- matrix(rnorm(48 * 6), nrow = 48, ncol = 6)
  ref <- t(as.matrix(s_host %*% b_host))

  sp_key <- paste0("metal-sp-", sample.int(1e6, 1))
  y_key <- paste0("metal-y-", sample.int(1e6, 1))
  prod_key <- paste0("metal-prod-", sample.int(1e6, 1))
  out_key <- paste0("metal-out-", sample.int(1e6, 1))

  on.exit(amatrix_metal_sparse_resident_drop(sp_key), add = TRUE)
  on.exit(amatrix_metal_resident_drop(y_key), add = TRUE)
  on.exit(amatrix_metal_resident_drop(prod_key), add = TRUE)
  on.exit(amatrix_metal_resident_drop(out_key), add = TRUE)

  amatrix_metal_sparse_resident_store(sp_key, s_host)
  amatrix_metal_resident_store(y_key, b_host)

  value <- amatrix_metal_backend()$spmm_resident_key(
    sp_key,
    y_key,
    prod_key,
    trans_lhs = FALSE,
    defer = TRUE
  )
  expect_null(value)

  amatrix_metal_backend()$transpose_resident(prod_key, out_key)

  expect_true(amatrix_metal_resident_has(prod_key))
  expect_true(amatrix_metal_resident_has(out_key))
  expect_equal(amatrix_metal_resident_materialize(out_key), ref, tolerance = 1e-5)
})
