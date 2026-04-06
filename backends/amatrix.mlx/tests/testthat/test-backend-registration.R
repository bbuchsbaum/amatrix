test_that("mlx backend registers with amatrix", {
  expect_true("mlx" %in% amatrix::amatrix_backend_names())
})

test_that("mlx backend advertises dense-first capabilities", {
  backend <- amatrix_mlx_backend()

  expect_identical(backend$available(), amatrix_mlx_is_available())
  expect_false(backend$supports("matmul", amatrix::adgeMatrix(matrix(1:4, nrow = 2))))
  expect_true(backend$supports("matmul", amatrix::adgeMatrix(matrix(1, nrow = 128, ncol = 128), precision = "fast")))
  expect_true(backend$supports("ewise", amatrix::adgeMatrix(matrix(1:4, nrow = 2), precision = "fast")))
  expect_true(backend$supports("qr", amatrix::adgeMatrix(matrix(1, nrow = 256, ncol = 256), precision = "fast")))
  expect_false(backend$supports("solve", amatrix::adgeMatrix(matrix(1:4, nrow = 2))))
  expect_false(backend$supports("matmul", amatrix::adgeMatrix(matrix(1, nrow = 128, ncol = 128), precision = "strict")))
  expect_false(backend$supports("matmul", amatrix::adgCMatrix(matrix(c(1, 0, 0, 1), nrow = 2))))
})

test_that("mlx capability list is stable and explicit", {
  expect_identical(
    amatrix_mlx_capabilities(),
    c("matmul", "crossprod", "tcrossprod", "ewise", "rowSums", "colSums",
      "qr", "rsvd", "chol", "chol_gpu", "batched_trsm", "eigen")
  )
})

test_that("mlx registration exposes capabilities through core helpers", {
  expect_identical(
    amatrix::amatrix_backend_capabilities("mlx"),
    amatrix_mlx_capabilities()
  )

  status <- amatrix::amatrix_backend_status("mlx")
  expect_identical(status$name, "mlx")
  expect_identical(status$available, amatrix_mlx_is_available())
  expect_identical(status$precision_modes, "fast")
  expect_identical(status$capabilities, paste(amatrix_mlx_capabilities(), collapse = ","))
})

test_that("mlx bridge boundary reports coherent native status", {
  info <- amatrix_mlx_bridge_info()

  expect_true(info$compiled)
  expect_identical(info$native, amatrix_mlx_native_available())
  expect_identical(info$available, amatrix_mlx_is_available())
  expect_true(info$engine %in% c("mock-c-bridge", "mlx-c"))
  expect_identical(info$capabilities, amatrix_mlx_capabilities())
})

test_that("mlx bridge boundary is callable", {
  backend <- amatrix_mlx_backend()
  x <- matrix(c(1, 2, 3, 4), nrow = 2)

  expect_equal(backend$matmul(x, diag(2)), x)
  expect_equal(backend$crossprod(x), crossprod(x))
  expect_equal(backend$crossprod(x, x), crossprod(x, x))
  expect_equal(backend$tcrossprod(x), tcrossprod(x))
  expect_equal(backend$tcrossprod(x, x), tcrossprod(x, x))
  expect_equal(backend$ewise(x, lhs = x, rhs = 2, op = "*"), x * 2, tolerance = 1e-5)
  expect_equal(backend$rowSums(x), rowSums(x), tolerance = 1e-5)
  expect_equal(backend$colSums(x), colSums(x), tolerance = 1e-5)
  fac_base <- qr(x)
  qr_fit <- backend$qr(x)
  expect_true(is.list(qr_fit))
  expect_equal(unname(qr_fit$r), unname(qr.R(fac_base)), tolerance = 1e-4)
  expect_true(is.character(qr_fit$q_key))
  expect_length(qr_fit$q_key, 1L)
  expect_null(qr_fit$q)
  expect_equal(unname(amatrix_mlx_resident_materialize(qr_fit$q_key)), unname(qr.Q(fac_base)), tolerance = 1e-4)
  expect_null(qr_fit$factor)
  expect_identical(qr_fit$factor_source, "reconstructable")
  expect_identical(qr_fit$representation, "explicit_qr")
})

test_that("mlx qr bridge can stamp compact representation mode", {
  old <- options(amatrix.mlx.qr_helper_mode = "compact")
  on.exit(options(old), add = TRUE)

  qr_fit <- amatrix_mlx_qr(matrix(c(1, 2, 3, 4), nrow = 2))

  expect_identical(qr_fit$representation, "mlx_compact_qr")
  expect_null(qr_fit$factor)
  expect_true(is.function(qr_fit$factor_builder))
  expect_identical(qr_fit$factor_source, "host_compact")
})

test_that("mlx compact qr can use tsqr-blocked factorization", {
  old <- options(amatrix.mlx.qr_helper_mode = "compact", amatrix.mlx.qr_compact_method = "tsqr", amatrix.mlx.qr_tsqr_block_rows = 4L)
  on.exit(options(old), add = TRUE)

  x <- matrix(rnorm(24 * 4), nrow = 24, ncol = 4)
  qr_fit <- amatrix_mlx_qr(x)

  expect_identical(qr_fit$representation, "mlx_compact_qr")
  expect_identical(qr_fit$factor_source, "tsqr_blocked")
  expect_true(inherits(qr_fit$factor, "amatrix_mlx_tsqr_factor"))
  expect_type(qr_fit$factor$top_q_key, "character")
  expect_length(qr_fit$factor$top_q_key, 1L)
  expect_type(qr_fit$factor$top_r_key, "character")
  expect_length(qr_fit$factor$top_r_key, 1L)
  expect_length(qr_fit$factor$block_q_keys, 6L)
  expect_type(qr_fit$factor$r_stack_key, "character")
  expect_length(qr_fit$factor$r_stack_key, 1L)
  expect_null(qr_fit$q)
  expect_null(qr_fit$q_key)
  expect_equal(
    unname(abs(diag(amatrix_mlx_resident_materialize(qr_fit$factor$top_r_key)))),
    unname(abs(diag(qr.R(base::qr(x))))),
    tolerance = 1e-5
  )
})

test_that("mlx qr cache signature distinguishes helper strategies", {
  old <- options(amatrix.mlx.qr_tsqr_block_rows = 64L)
  on.exit(options(old), add = TRUE)

  native <- {
    options(amatrix.mlx.qr_helper_mode = "native")
    amatrix_mlx_qr_cache_signature(c(1024L, 128L))
  }
  compact_bridge <- {
    options(amatrix.mlx.qr_helper_mode = "compact", amatrix.mlx.qr_compact_method = "bridge")
    amatrix_mlx_qr_cache_signature(c(1024L, 128L))
  }
  compact_tsqr <- {
    options(amatrix.mlx.qr_helper_mode = "compact", amatrix.mlx.qr_compact_method = "tsqr")
    amatrix_mlx_qr_cache_signature(c(1024L, 128L))
  }

  expect_identical(native, "mlx:native")
  expect_identical(compact_bridge, "mlx:compact:bridge")
  expect_identical(compact_tsqr, "mlx:compact:tsqr:64")
})

test_that("mlx availability can be enabled for routing tests", {
  old <- getOption("amatrix.mlx.available")
  options(amatrix.mlx.available = TRUE)
  on.exit(options(amatrix.mlx.available = old), add = TRUE)

  x <- amatrix::adgeMatrix(matrix(1:4, nrow = 2), preferred_backend = "mlx", precision = "fast")
  dense_plan <- amatrix::amatrix_backend_plan(x, "matmul", y = diag(2))
  unsupported_plan <- amatrix::amatrix_backend_plan(x, "solve")
  sparse_plan <- amatrix::amatrix_backend_plan(
    amatrix::adgCMatrix(matrix(c(1, 0, 0, 1), nrow = 2), preferred_backend = "mlx", precision = "fast"),
    "matmul",
    y = diag(2)
  )

  expect_true(amatrix_mlx_is_available())
  expect_identical(dense_plan$chosen, "mlx")
  expect_identical(unsupported_plan$chosen, "cpu")
  expect_identical(sparse_plan$chosen, "cpu")
})

test_that("forced availability bypasses size heuristics for backend tests", {
  old <- getOption("amatrix.mlx.available")
  options(amatrix.mlx.available = TRUE)
  on.exit(options(amatrix.mlx.available = old), add = TRUE)

  backend <- amatrix_mlx_backend()
  x <- amatrix::adgeMatrix(matrix(1:4, nrow = 2), precision = "fast")

  expect_true(backend$supports("matmul", x))
  expect_true(backend$supports("crossprod", x))
  expect_true(backend$supports("tcrossprod", x))
  expect_true(backend$supports("qr", x))
})
