test_that("arrayfire backend registers with amatrix", {
  expect_true("arrayfire" %in% amatrix::amatrix_backend_names())
})

test_that("arrayfire backend advertises dense-first capabilities", {
  backend <- amatrix_arrayfire_backend()

  expect_identical(backend$available(), amatrix_arrayfire_is_available())
  expect_false(backend$supports("matmul", amatrix::adgeMatrix(matrix(1:4, nrow = 2))))
  expect_true(backend$supports("matmul", amatrix::adgeMatrix(matrix(1, nrow = 512, ncol = 512), precision = "fast")))
  expect_false(backend$supports("ewise", amatrix::adgeMatrix(matrix(1:4, nrow = 2), precision = "fast")))
  expect_false(backend$supports("ewise", amatrix::adgeMatrix(matrix(1, nrow = 1024, ncol = 1024), precision = "fast")))
  expect_false(backend$supports("rowSums", amatrix::adgeMatrix(matrix(1, nrow = 1024, ncol = 1024), precision = "fast")))
  expect_false(backend$supports("colSums", amatrix::adgeMatrix(matrix(1, nrow = 1024, ncol = 1024), precision = "fast")))
  old_backend <- amatrix_arrayfire_active_backend()
  amatrix_arrayfire_set_backend("cpu")
  on.exit(amatrix_arrayfire_set_backend(if (identical(old_backend, 4L)) "opencl" else "cpu"), add = TRUE)
  expect_true(backend$supports("qr", amatrix::adgeMatrix(matrix(1, nrow = 512, ncol = 512), precision = "fast")))
  expect_false(backend$supports("solve", amatrix::adgeMatrix(matrix(1:4, nrow = 2))))
  expect_false(backend$supports("matmul", amatrix::adgeMatrix(matrix(1, nrow = 512, ncol = 512), precision = "strict")))
  expect_false(backend$supports("matmul", amatrix::adgCMatrix(matrix(c(1, 0, 0, 1), nrow = 2))))
})

test_that("arrayfire capability list is stable and explicit", {
  expect_identical(
    amatrix_arrayfire_capabilities(),
    c("matmul", "crossprod", "tcrossprod", "ewise", "rowSums", "colSums", "qr")
  )
})

test_that("arrayfire registration exposes capabilities through core helpers", {
  expect_identical(
    amatrix::amatrix_backend_capabilities("arrayfire"),
    amatrix_arrayfire_capabilities()
  )

  status <- amatrix::amatrix_backend_status("arrayfire")
  expect_identical(status$name, "arrayfire")
  expect_identical(status$available, amatrix_arrayfire_is_available())
  expect_identical(status$precision_modes, "fast")
  expect_identical(status$capabilities, paste(amatrix_arrayfire_capabilities(), collapse = ","))
})

test_that("arrayfire bridge boundary reports coherent native status", {
  info <- amatrix_arrayfire_bridge_info()

  expect_true(info$compiled)
  expect_identical(info$native, amatrix_arrayfire_native_available())
  expect_identical(info$available, amatrix_arrayfire_is_available())
  expect_true(info$engine %in% c("mock-c-bridge", "arrayfire-c"))
  expect_identical(info$capabilities, amatrix_arrayfire_capabilities())
})

test_that("arrayfire bridge boundary is callable", {
  backend <- amatrix_arrayfire_backend()
  x <- matrix(c(1, 2, 3, 4), nrow = 2)
  old_backend <- amatrix_arrayfire_active_backend()
  amatrix_arrayfire_set_backend("cpu")
  on.exit(amatrix_arrayfire_set_backend(if (identical(old_backend, 4L)) "opencl" else "cpu"), add = TRUE)

  expect_equal(backend$matmul(x, diag(2)), x, tolerance = 1e-5)
  expect_equal(backend$crossprod(x), crossprod(x), tolerance = 1e-5)
  expect_equal(backend$crossprod(x, x), crossprod(x, x), tolerance = 1e-5)
  expect_equal(backend$tcrossprod(x), tcrossprod(x), tolerance = 1e-5)
  expect_equal(backend$tcrossprod(x, x), tcrossprod(x, x), tolerance = 1e-5)
  expect_equal(backend$ewise(x, lhs = x, rhs = 2, op = "*"), x * 2, tolerance = 1e-5)
  expect_equal(backend$ewise(x, lhs = x, rhs = x, op = "+"), x + x, tolerance = 1e-5)
  expect_equal(backend$rowSums(x), rowSums(x), tolerance = 1e-5)
  expect_equal(backend$colSums(x), colSums(x), tolerance = 1e-5)
  fac_base <- qr(x)
  qr_fit <- backend$qr(x)
  expect_true(is.list(qr_fit))
  expect_equal(unname(qr_fit$q), unname(qr.Q(fac_base)), tolerance = 1e-4)
  expect_equal(unname(qr_fit$r), unname(qr.R(fac_base)), tolerance = 1e-4)
  expect_true(amatrix_arrayfire_native_available())
})

test_that("arrayfire availability can be enabled for routing tests", {
  old <- getOption("amatrix.arrayfire.available")
  options(amatrix.arrayfire.available = TRUE)
  on.exit(options(amatrix.arrayfire.available = old), add = TRUE)

  x <- amatrix::adgeMatrix(matrix(1:4, nrow = 2), preferred_backend = "arrayfire", precision = "fast")
  dense_plan <- amatrix::amatrix_backend_plan(x, "matmul", y = diag(2))
  unsupported_plan <- amatrix::amatrix_backend_plan(x, "solve")
  sparse_plan <- amatrix::amatrix_backend_plan(
    amatrix::adgCMatrix(matrix(c(1, 0, 0, 1), nrow = 2), preferred_backend = "arrayfire", precision = "fast"),
    "matmul",
    y = diag(2)
  )

  expect_true(amatrix_arrayfire_is_available())
  expect_identical(dense_plan$chosen, "arrayfire")
  expect_identical(unsupported_plan$chosen, "cpu")
  expect_identical(sparse_plan$chosen, "cpu")
})

test_that("forced availability bypasses size heuristics for backend tests", {
  old <- getOption("amatrix.arrayfire.available")
  options(amatrix.arrayfire.available = TRUE)
  on.exit(options(amatrix.arrayfire.available = old), add = TRUE)

  backend <- amatrix_arrayfire_backend()
  x <- amatrix::adgeMatrix(matrix(1:4, nrow = 2), precision = "fast")
  old_backend <- amatrix_arrayfire_active_backend()
  amatrix_arrayfire_set_backend("cpu")
  on.exit(amatrix_arrayfire_set_backend(if (identical(old_backend, 4L)) "opencl" else "cpu"), add = TRUE)

  expect_true(backend$supports("matmul", x))
  expect_true(backend$supports("ewise", x))
  expect_true(backend$supports("crossprod", x))
  expect_true(backend$supports("tcrossprod", x))
  expect_true(backend$supports("rowSums", x))
  expect_true(backend$supports("colSums", x))
  expect_true(backend$supports("qr", x))
})
