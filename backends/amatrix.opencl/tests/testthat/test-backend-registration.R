test_that("opencl backend registers with amatrix", {
  amatrix_opencl_register(overwrite = TRUE)
  expect_true("opencl" %in% amatrix::amatrix_backend_names())
})

test_that("opencl backend advertises dense-first scaffold capabilities", {
  backend <- amatrix_opencl_backend()
  small_fast <- amatrix::adgeMatrix(matrix(1:4, nrow = 2), precision = "fast")

  expect_false(backend$supports("matmul", amatrix::adgeMatrix(matrix(1:4, nrow = 2), precision = "fast")))

  old <- options(amatrix.opencl.available = TRUE)
  on.exit(options(old), add = TRUE)

  expect_identical(backend$available(), amatrix_opencl_is_available())
  expect_false(backend$supports("matmul", amatrix::adgeMatrix(matrix(1:4, nrow = 2), precision = "fast")))
  expect_false(backend$supports("matmul", amatrix::adgeMatrix(matrix(1, nrow = 256, ncol = 32), precision = "fast")))
  expect_true(backend$supports("matmul", amatrix::adgeMatrix(matrix(1, nrow = 128, ncol = 128), precision = "fast")))
  expect_true(backend$supports("ewise", amatrix::adgeMatrix(matrix(1:4, nrow = 2), precision = "fast")))
  expect_true(backend$supports("solve", amatrix::adgeMatrix(diag(2), precision = "fast")))
  expect_true(backend$supports("chol", amatrix::adgeMatrix(diag(2), precision = "fast")))
  expect_true(backend$supports("qr", amatrix::adgeMatrix(matrix(1:4, nrow = 2), precision = "fast")))
  expect_true(backend$supports("svd", small_fast))
  expect_true(backend$supports("eigen", amatrix::adgeMatrix(diag(2), precision = "fast")))
  expect_false(backend$supports("eigen", amatrix::adgeMatrix(matrix(1:6, nrow = 2), precision = "fast")))
  expect_true(backend$supports("covariance", amatrix::adgeMatrix(matrix(1:4, nrow = 2), precision = "fast")))
  expect_false(backend$supports("matmul", amatrix::adgeMatrix(matrix(1, nrow = 128, ncol = 128), precision = "strict")))
  expect_false(backend$supports("svd", amatrix::adgeMatrix(matrix(1:4, nrow = 2), precision = "strict")))
  expect_false(backend$supports("matmul", amatrix::adgCMatrix(matrix(c(1, 0, 0, 1), nrow = 2), precision = "fast")))
})

test_that("opencl factor GPU support is opt-in and size-gated", {
  backend <- amatrix_opencl_backend()
  x <- amatrix::adgeMatrix(diag(2), precision = "fast")
  nonsym <- amatrix::adgeMatrix(matrix(c(4, 1, 0, 3), nrow = 2), precision = "fast")
  old <- options(
    amatrix.opencl.available = TRUE,
    amatrix.opencl.factor_gpu = TRUE,
    amatrix.opencl.factor_min_dim = 1L
  )
  on.exit(options(old), add = TRUE)

  expect_true(is.function(backend$supports_resident))
  expect_true(backend$supports_resident("solve", x))
  expect_false(backend$supports_resident("solve", nonsym))
  expect_true(backend$supports_resident("chol", x))
})

test_that("opencl capability list is stable and explicit", {
  backend <- amatrix_opencl_backend()
  expect_identical(
    amatrix_opencl_capabilities(),
    c("matmul", "crossprod", "tcrossprod", "ewise",
      "broadcast_ewise", "rowSums", "colSums", "solve", "chol", "qr", "svd", "eigen", "covariance")
  )
  expect_true("sparse_spmm" %in% backend$features())
})

test_that("opencl registration exposes capabilities through core helpers", {
  old <- options(amatrix.opencl.available = TRUE)
  on.exit(options(old), add = TRUE)

  amatrix_opencl_register(overwrite = TRUE)

  expect_identical(
    amatrix::amatrix_backend_capabilities("opencl"),
    amatrix_opencl_capabilities()
  )

  status <- amatrix::amatrix_backend_status("opencl")
  expect_identical(status$name, "opencl")
  expect_identical(status$available, amatrix_opencl_is_available())
  expect_identical(status$precision_modes, "fast")
  expect_identical(status$capabilities, paste(amatrix_opencl_capabilities(), collapse = ","))
})

test_that("opencl bridge boundary reports coherent scaffold status", {
  old_env <- Sys.getenv("AMATRIX_OPENCL_PROBE_GPU", unset = NA_character_)
  Sys.unsetenv("AMATRIX_OPENCL_PROBE_GPU")
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("AMATRIX_OPENCL_PROBE_GPU") else Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = old_env)
  }, add = TRUE)

  info <- amatrix_opencl_bridge_info()

  expect_type(info$compiled, "logical")
  expect_type(info$clblast, "logical")
  expect_identical(info$native, amatrix_opencl_native_available())
  expect_identical(info$available, amatrix_opencl_is_available())
  expect_true(info$engine %in% c("mock-c-bridge", "opencl-runtime", "opencl-clblast-scaffold"))
  expect_identical(info$capabilities, amatrix_opencl_capabilities())
})

test_that("opencl bridge boundary is callable in scaffold mode", {
  old <- options(amatrix.opencl.available = TRUE)
  on.exit(options(old), add = TRUE)

  backend <- amatrix_opencl_backend()
  x <- matrix(c(1, 2, 3, 4), nrow = 2)
  spd <- crossprod(matrix(c(1, 2, 3, 5), nrow = 2)) + diag(2)
  rhs <- matrix(c(1, -1, 2, 0), nrow = 2)

  expect_equal(backend$matmul(x, diag(2)), x)
  expect_equal(backend$crossprod(x), crossprod(x))
  expect_equal(backend$crossprod(x, x), crossprod(x, x))
  expect_equal(backend$tcrossprod(x), tcrossprod(x))
  expect_equal(backend$tcrossprod(x, x), tcrossprod(x, x))
  expect_equal(backend$ewise(x, lhs = x, rhs = 2, op = "*"), x * 2, tolerance = 1e-12)
  expect_equal(backend$ewise(x, lhs = x, rhs = matrix(1, nrow = 2, ncol = 2), op = "+"), x + 1, tolerance = 1e-12)
  expect_equal(backend$broadcast_ewise(x, lhs = x, v = c(10, 20), margin = 1L, op = "+"), sweep(x, 1L, c(10, 20), "+"), tolerance = 1e-12)
  expect_equal(backend$rowSums(x), rowSums(x), tolerance = 1e-12)
  expect_equal(backend$colSums(x), colSums(x), tolerance = 1e-12)
  expect_equal(backend$chol(spd), chol(spd), tolerance = 1e-12)
  qr_template <- amatrix::adgeMatrix(x, preferred_backend = "opencl", precision = "fast")
  qr_fit <- amatrix:::.amatrix_wrap_qr(backend$qr(x), qr_template, method = "fast")
  qr_ref <- qr(x)
  svd_fit <- backend$svd(x, nu = 2L, nv = 2L)
  svd_ref <- base::svd(x, nu = 2L, nv = 2L)
  expect_equal(backend$solve(spd), solve(spd), tolerance = 1e-12)
  expect_equal(backend$solve(spd, rhs), solve(spd, rhs), tolerance = 1e-12)
  expect_equal(backend$solve(matrix(c(4, 1, 0, 3), nrow = 2), rhs), solve(matrix(c(4, 1, 0, 3), nrow = 2), rhs), tolerance = 1e-12)
  expect_equal(backend$eigen(spd, symmetric = TRUE)$values, eigen(spd, symmetric = TRUE)$values, tolerance = 1e-12)
  expect_equal(
    sort(Re(backend$eigen(matrix(c(2, 1, 0, 3), nrow = 2), symmetric = FALSE)$values)),
    sort(Re(eigen(matrix(c(2, 1, 0, 3), nrow = 2), symmetric = FALSE)$values)),
    tolerance = 1e-12
  )
  expect_equal(backend$chol_solve_factor(chol(spd), rhs), solve(spd, rhs), tolerance = 1e-12)
  expect_equal(
    backend$solve_triangular_factor(chol(spd), rhs, lower = FALSE, transpose = FALSE),
    backsolve(chol(spd), rhs),
    tolerance = 1e-12
  )
  expect_equal(backend$covariance(x), stats::cov(x), tolerance = 1e-12)
  expect_true(nzchar(amatrix:::.amatrix_qr_q_key(qr_fit)))
  expect_equal(unname(amatrix:::.amatrix_qr_q(qr_fit)), unname(qr.Q(qr_ref)), tolerance = 1e-12)
  expect_equal(unname(amatrix:::.amatrix_qr_r(qr_fit)), unname(qr.R(qr_ref)), tolerance = 1e-12)
  expect_equal(svd_fit$d, svd_ref$d, tolerance = 1e-12)
  expect_equal(
    svd_fit$u %*% diag(svd_fit$d, nrow = length(svd_fit$d)) %*% t(svd_fit$v),
    x,
    tolerance = 1e-12
  )

  backend$resident_store("x", x)
  backend$resident_store("y", diag(2))
  backend$resident_store("spd", spd)
  backend$resident_store("rhs", rhs)
  on.exit({
    backend$resident_drop("out")
    backend$resident_drop("swept")
    backend$resident_drop("chol")
    backend$resident_drop("sol")
    backend$resident_drop("tri")
    backend$resident_drop("q")
    backend$resident_drop("x")
    backend$resident_drop("y")
    backend$resident_drop("spd")
    backend$resident_drop("rhs")
  }, add = TRUE)

  expect_true(backend$resident_has("x"))
  expect_equal(backend$resident_materialize("x"), x, tolerance = 1e-12)
  expect_equal(backend$rowSums_resident("x"), rowSums(x), tolerance = 1e-12)
  expect_equal(backend$colSums_resident("x"), colSums(x), tolerance = 1e-12)

  backend$matmul_resident("x", "y", "out")
  expect_equal(backend$resident_materialize("out"), x %*% diag(2), tolerance = 1e-12)

  backend$broadcast_ewise_resident("x", c(10, 20), 1L, "+", "swept")
  expect_equal(backend$resident_materialize("swept"), sweep(x, 1L, c(10, 20), "+"), tolerance = 1e-12)

  backend$chol_resident("spd", "chol")
  backend$qr_Q_resident("x", "q")
  backend$solve_resident("spd", "rhs", "sol")
  backend$solve_triangular_resident("chol", "rhs", "tri", lower = FALSE, transpose = FALSE)
  expect_equal(backend$resident_materialize("chol"), chol(spd), tolerance = 1e-12)
  expect_equal(backend$resident_materialize("q"), qr.Q(qr_ref), tolerance = 1e-12)
  expect_equal(backend$resident_materialize("sol"), solve(spd, rhs), tolerance = 1e-12)
  expect_equal(backend$resident_materialize("tri"), backsolve(chol(spd), rhs), tolerance = 1e-12)
})

test_that("opencl sparse routing is threshold-gated and host-backed", {
  backend <- amatrix_opencl_backend()
  sparse <- amatrix::adgCMatrix(
    Matrix::rsparsematrix(64, 64, density = 0.05),
    preferred_backend = "opencl",
    precision = "fast"
  )
  rhs <- matrix(rnorm(64 * 8), nrow = 64, ncol = 8)

  old <- options(
    amatrix.opencl.available = TRUE,
    amatrix.opencl.spmv_min_nnz = Inf,
    amatrix.opencl.spmm_min_nnz = Inf
  )
  on.exit(options(old), add = TRUE)

  expect_false(backend$supports("matmul", sparse, y = rhs))

  options(amatrix.opencl.spmv_min_nnz = 1L, amatrix.opencl.spmm_min_nnz = 1L)

  expect_true(backend$supports("matmul", sparse, y = rhs))
  expect_true(backend$supports_resident("matmul", sparse, y = rhs))
  expect_true(backend$supports("crossprod", sparse, y = rhs))
  expect_true(backend$supports("tcrossprod", sparse, y = t(rhs)))

  host_sparse <- as(amatrix::amatrix_materialize_host(sparse), "dgCMatrix")
  expect_equal(backend$matmul(host_sparse, rhs), as.matrix(host_sparse %*% rhs), tolerance = 1e-12)

  backend$sparse_resident_store("sp", host_sparse)
  backend$resident_store("rhs", rhs)
  on.exit({
    backend$resident_drop("rhs")
    backend$resident_drop("out")
    backend$sparse_resident_drop("sp")
  }, add = TRUE)

  expect_true(backend$sparse_resident_has("sp"))
  backend$spmm_resident_key("sp", "rhs", "out", trans_lhs = FALSE, defer = TRUE)
  expect_equal(backend$resident_materialize("out"), as.matrix(host_sparse %*% rhs), tolerance = 1e-12)
})

test_that("opencl diagnostics remain safe without probe enablement", {
  old_env <- Sys.getenv("AMATRIX_OPENCL_PROBE_GPU", unset = NA_character_)
  Sys.unsetenv("AMATRIX_OPENCL_PROBE_GPU")
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("AMATRIX_OPENCL_PROBE_GPU") else Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = old_env)
  }, add = TRUE)

  diag <- amatrix_opencl_diagnostics()

  expect_type(diag$compiled, "logical")
  expect_type(diag$clblast, "logical")
  expect_false(diag$probe_enabled)
  expect_true(diag$engine %in% c("mock-c-bridge", "opencl-runtime", "opencl-clblast-scaffold"))
  expect_true(diag$resident_entries >= 0L)
  expect_true(diag$resident_device_entries >= 0L)
  expect_true(diag$resident_host_entries >= 0L)
  expect_type(diag$device_name, "character")
})

test_that("opencl package load plus backend listing stays safe with probe disabled", {
  old_env <- Sys.getenv("AMATRIX_OPENCL_PROBE_GPU", unset = NA_character_)
  Sys.unsetenv("AMATRIX_OPENCL_PROBE_GPU")
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("AMATRIX_OPENCL_PROBE_GPU") else Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = old_env)
  }, add = TRUE)

  expect_no_error(loadNamespace("amatrix.opencl"))
  expect_no_error(amatrix::amatrix_backend_names())
})

test_that("opencl resident registry round-trips on device when native runtime is available", {
  old_env <- Sys.getenv("AMATRIX_OPENCL_PROBE_GPU", unset = NA_character_)
  Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = "1")
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("AMATRIX_OPENCL_PROBE_GPU") else Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = old_env)
  }, add = TRUE)

  skip_if_not(isTRUE(amatrix_opencl_native_available(force = TRUE)))

  backend <- amatrix_opencl_backend()
  x <- matrix(as.double(1:6), nrow = 2)

  backend$resident_store("gpu_x", x)
  on.exit(backend$resident_drop("gpu_x"), add = TRUE)

  diag <- amatrix_opencl_diagnostics()
  expect_true(isTRUE(backend$resident_has("gpu_x")))
  expect_true(diag$resident_entries >= 1L)
  expect_true(diag$resident_device_entries >= 1L)
  expect_equal(backend$resident_materialize("gpu_x"), x, tolerance = 1e-12)
})

test_that("opencl resident matmul stays on device when CLBlast runtime is available", {
  old_env <- Sys.getenv("AMATRIX_OPENCL_PROBE_GPU", unset = NA_character_)
  Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = "1")
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("AMATRIX_OPENCL_PROBE_GPU") else Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = old_env)
  }, add = TRUE)

  skip_if_not(isTRUE(amatrix_opencl_native_available(force = TRUE)))
  skip_if_not(isTRUE(amatrix_opencl_bridge_info()$clblast))

  backend <- amatrix_opencl_backend()
  x <- matrix(as.double(1:6), nrow = 2)
  y <- diag(3)

  backend$resident_store("gpu_x", x)
  backend$resident_store("gpu_y", y)
  on.exit({
    backend$resident_drop("gpu_out")
    backend$resident_drop("gpu_x")
    backend$resident_drop("gpu_y")
  }, add = TRUE)

  backend$matmul_resident("gpu_x", "gpu_y", "gpu_out", defer = TRUE)
  diag <- amatrix_opencl_diagnostics()

  expect_true(isTRUE(backend$resident_has("gpu_out")))
  expect_true(diag$resident_device_entries >= 3L)
  expect_equal(backend$resident_materialize("gpu_out"), x %*% y, tolerance = 1e-5)
})

test_that("opencl resident ewise, broadcast, and reductions stay coherent on device", {
  old_env <- Sys.getenv("AMATRIX_OPENCL_PROBE_GPU", unset = NA_character_)
  Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = "1")
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("AMATRIX_OPENCL_PROBE_GPU") else Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = old_env)
  }, add = TRUE)

  skip_if_not(isTRUE(amatrix_opencl_native_available(force = TRUE)))
  skip_if_not(isTRUE(amatrix_opencl_bridge_info()$clblast))

  backend <- amatrix_opencl_backend()
  x <- matrix(as.double(1:6), nrow = 2)
  y <- matrix(1, nrow = 2, ncol = 3)

  backend$resident_store("gpu_x", x)
  backend$resident_store("gpu_y", y)
  on.exit({
    backend$resident_drop("gpu_add")
    backend$resident_drop("gpu_row")
    backend$resident_drop("gpu_swept")
    backend$resident_drop("gpu_scaled")
    backend$resident_drop("gpu_x")
    backend$resident_drop("gpu_y")
  }, add = TRUE)

  backend$ewise_resident("gpu_x", "gpu_y", "+", "gpu_add", defer = TRUE)
  backend$broadcast_ewise_resident("gpu_x", c(10, 20), 1L, "+", "gpu_swept", defer = TRUE)
  backend$rowSums_resident_key("gpu_x", "gpu_row")
  backend$broadcast_ewise_resident_key("gpu_x", "gpu_row", 1L, "/", "gpu_scaled", defer = TRUE)
  backend$broadcast_ewise_resident_inplace("gpu_x", c(2, 3), 1L, "*")
  backend$broadcast_ewise_resident_inplace_key("gpu_x", "gpu_row", 1L, "/")

  diag <- amatrix_opencl_diagnostics()
  expect_true(diag$resident_device_entries >= 6L)
  expect_equal(backend$resident_materialize("gpu_add"), x + y, tolerance = 1e-5)
  expect_equal(backend$resident_materialize("gpu_swept"), sweep(x, 1L, c(10, 20), "+"), tolerance = 1e-5)
  expect_equal(backend$resident_materialize("gpu_row"), matrix(rowSums(x), ncol = 1L), tolerance = 1e-5)
  expect_equal(backend$resident_materialize("gpu_scaled"), sweep(x, 1L, rowSums(x), "/"), tolerance = 1e-5)
  expect_equal(
    backend$resident_materialize("gpu_x"),
    sweep(sweep(x, 1L, c(2, 3), "*"), 1L, rowSums(x), "/"),
    tolerance = 1e-5
  )
  expect_equal(
    backend$rowSums_resident("gpu_x"),
    rowSums(sweep(sweep(x, 1L, c(2, 3), "*"), 1L, rowSums(x), "/")),
    tolerance = 1e-5
  )
  expect_equal(
    backend$colSums_resident("gpu_x"),
    colSums(sweep(sweep(x, 1L, c(2, 3), "*"), 1L, rowSums(x), "/")),
    tolerance = 1e-5
  )
})

test_that("opencl rectangular self-crossprod and self-tcrossprod stay symmetric", {
  old_env <- Sys.getenv("AMATRIX_OPENCL_PROBE_GPU", unset = NA_character_)
  Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = "1")
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("AMATRIX_OPENCL_PROBE_GPU") else Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = old_env)
  }, add = TRUE)

  skip_if_not(isTRUE(amatrix_opencl_native_available(force = TRUE)))
  skip_if_not(isTRUE(amatrix_opencl_bridge_info()$clblast))

  backend <- amatrix_opencl_backend()
  x <- matrix(as.double(1:15), nrow = 3, ncol = 5)

  cp <- backend$crossprod(x)
  tcp <- backend$tcrossprod(x)

  expect_equal(cp, crossprod(x), tolerance = 1e-5)
  expect_equal(tcp, tcrossprod(x), tolerance = 1e-5)
  expect_equal(cp, t(cp), tolerance = 1e-7)
  expect_equal(tcp, t(tcp), tolerance = 1e-7)
})

test_that("opencl resident_handle sinkhorn chain converges for 50 iterations", {
  old_env <- Sys.getenv("AMATRIX_OPENCL_PROBE_GPU", unset = NA_character_)
  Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = "1")
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("AMATRIX_OPENCL_PROBE_GPU") else Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = old_env)
  }, add = TRUE)

  skip_if_not(isTRUE(amatrix_opencl_native_available(force = TRUE)))

  set.seed(42)
  a <- matrix(rexp(36), nrow = 6) + 1e-3
  h <- amatrix:::resident_handle(amatrix::adgeMatrix(a, preferred_backend = "opencl", precision = "fast"))
  on.exit({
    if (isTRUE(h$active) && !is.null(h$resident_key)) {
      try(amatrix_opencl_backend()$resident_drop(h$resident_key), silent = TRUE)
    }
  }, add = TRUE)

  for (iter in seq_len(50L)) {
    rs <- amatrix:::rh_rowSums(h)
    amatrix::am_sweep_inplace(h, 1L, pmax(rs, 1e-15), "/")
    cs <- amatrix:::rh_colSums(h)
    amatrix::am_sweep_inplace(h, 2L, pmax(cs, 1e-15), "/")
  }

  mat <- as.matrix(h)
  expect_lt(max(abs(rowSums(mat) - 1)), 1e-4)
  expect_lt(max(abs(colSums(mat) - 1)), 1e-4)
})
