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
    c("matmul", "crossprod", "tcrossprod", "ewise", "rowSums", "colSums",
      "qr", "rsvd", "chol", "solve", "covariance", "svd")
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
  dense_plan     <- amatrix::amatrix_backend_plan(x, "matmul", y = diag(2))
  # Use an op not in AF capabilities at all — must always fall back to cpu.
  not_in_caps_plan <- amatrix::amatrix_backend_plan(x, "fft")
  sparse_plan <- amatrix::amatrix_backend_plan(
    amatrix::adgCMatrix(matrix(c(1, 0, 0, 1), nrow = 2), preferred_backend = "arrayfire", precision = "fast"),
    "matmul",
    y = diag(2)
  )

  expect_true(amatrix_arrayfire_is_available())
  expect_identical(dense_plan$chosen, "arrayfire")
  expect_identical(not_in_caps_plan$chosen, "cpu")
  expect_identical(sparse_plan$chosen, "cpu")
})

test_that("bdc_svd gives accurate reconstruction for square and rectangular matrices", {
  old_backend <- amatrix_arrayfire_active_backend()
  amatrix_arrayfire_set_backend("cpu")
  on.exit(amatrix_arrayfire_set_backend(
    if (identical(old_backend, 4L)) "opencl" else "cpu"), add = TRUE)

  # Square matrix
  set.seed(42)
  A <- matrix(rnorm(64 * 64), 64, 64); storage.mode(A) <- "double"
  res <- amatrix_arrayfire_bdc_svd(A, nu = 64L, nv = 64L)
  recon_err <- max(abs(res$u %*% diag(res$d) %*% t(res$v) - A)) /
               max(abs(A), 1e-12)
  expect_lt(recon_err, 1e-4, label = "bdc_svd reconstruction (square)")
  ref_d <- base::svd(A, nu = 0L, nv = 0L)$d
  expect_equal(sort(res$d, decreasing = TRUE), ref_d, tolerance = 1e-4,
               label = "bdc_svd singular values (square)")

  # Near-square tall matrix (aspect 2:1 — still routed to BDC by dispatcher)
  set.seed(7)
  B <- matrix(rnorm(80 * 48), 80, 48); storage.mode(B) <- "double"
  res2 <- amatrix_arrayfire_bdc_svd(B, nu = 48L, nv = 48L)
  recon2 <- res2$u %*% diag(res2$d) %*% t(res2$v)
  recon_err2 <- max(abs(recon2 - B)) / max(abs(B), 1e-12)
  expect_lt(recon_err2, 1e-4, label = "bdc_svd reconstruction (near-square)")

  # Thin SVD: nu=nv=5 only
  set.seed(3)
  C <- matrix(rnorm(64 * 64), 64, 64); storage.mode(C) <- "double"
  res3 <- amatrix_arrayfire_bdc_svd(C, nu = 5L, nv = 5L)
  expect_equal(length(res3$d), 5L)
  expect_equal(dim(res3$u), c(64L, 5L))
  expect_equal(dim(res3$v), c(64L, 5L))
  ref_d5 <- base::svd(C, nu = 0L, nv = 0L)$d[seq_len(5L)]
  expect_equal(sort(res3$d, decreasing = TRUE), ref_d5, tolerance = 1e-4,
               label = "bdc_svd top-5 singular values")
})

test_that("svd dispatcher routes square-ish matrices to bdc_svd when native unsafe", {
  old_probe <- getOption("amatrix.arrayfire.native_svd_available")
  options(amatrix.arrayfire.native_svd_available = FALSE)
  old_bdc_n <- getOption("amatrix.arrayfire.bdc_min_n")
  options(amatrix.arrayfire.bdc_min_n = 32L)   # lower threshold for fast test
  old_backend <- amatrix_arrayfire_active_backend()
  amatrix_arrayfire_set_backend("cpu")
  on.exit({
    options(amatrix.arrayfire.native_svd_available = old_probe)
    options(amatrix.arrayfire.bdc_min_n = old_bdc_n)
    amatrix_arrayfire_set_backend(if (identical(old_backend, 4L)) "opencl" else "cpu")
  }, add = TRUE)

  set.seed(99)
  x <- matrix(rnorm(64 * 64), 64, 64); storage.mode(x) <- "double"
  res <- amatrix_arrayfire_svd(x, nu = 64L, nv = 64L)
  ref_d <- base::svd(x, nu = 0L, nv = 0L)$d
  expect_equal(sort(res$d, decreasing = TRUE), ref_d, tolerance = 1e-4,
               label = "dispatcher BDC path singular values")

  # Tall-thin should still go to ts_svd (aspect >= 4 bypasses BDC)
  set.seed(11)
  y <- matrix(rnorm(256 * 16), 256, 16); storage.mode(y) <- "double"
  res_ts <- amatrix_arrayfire_svd(y, nu = 16L, nv = 16L)
  ref_ts <- base::svd(y, nu = 0L, nv = 0L)$d
  expect_equal(sort(res_ts$d, decreasing = TRUE), ref_ts, tolerance = 1e-4,
               label = "dispatcher ts_svd path (tall-thin)")
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
