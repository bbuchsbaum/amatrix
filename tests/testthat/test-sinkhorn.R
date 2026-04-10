.sinkhorn_positive_matrix <- function(n = 24L, seed = 42L) {
  set.seed(seed)
  mat <- exp(matrix(rnorm(n * n), nrow = n, ncol = n))
  storage.mode(mat) <- "double"
  mat
}

.sinkhorn_backend_spec <- function(name) {
  specs <- optional_backend_specs()
  specs[[match(name, vapply(specs, `[[`, character(1), "backend"))]]
}

.register_optional_backend <- function(spec) {
  ns <- optional_backend_namespace(spec$package)
  get(spec$register_fun, envir = ns, inherits = FALSE)
}

.sinkhorn_recording_backend <- function(counter) {
  cpu <- amatrix:::.amatrix_cpu_backend()
  resident <- new.env(parent = emptyenv())

  bump <- function(name) {
    if (is.null(counter[[name]])) {
      counter[[name]] <- 0L
    }
    counter[[name]] <- counter[[name]] + 1L
  }

  list(
    capabilities = function() c("broadcast_ewise", "colSums", "rowSums"),
    features = function() "dense_f64",
    precision_modes = function() c("strict", "fast"),
    available = function() TRUE,
    supports = function(op, x, y = NULL) {
      inherits(x, "adgeMatrix") &&
        (x@precision %in% c("strict", "fast")) &&
        op %in% c("broadcast_ewise", "colSums", "rowSums")
    },
    matmul = function(x, y) cpu$matmul(x, y),
    crossprod = function(x, y = NULL) cpu$crossprod(x, y),
    tcrossprod = function(x, y = NULL) cpu$tcrossprod(x, y),
    ewise = function(lhs, rhs, op) cpu$ewise(lhs, rhs, op),
    rowSums = function(x, na.rm = FALSE, dims = 1L) {
      bump("rowSums")
      base::rowSums(as.matrix(x), na.rm = na.rm, dims = dims)
    },
    colSums = function(x, na.rm = FALSE, dims = 1L) {
      bump("colSums")
      base::colSums(as.matrix(x), na.rm = na.rm, dims = dims)
    },
    broadcast_ewise = function(lhs, v, margin, op) {
      bump("broadcast_ewise")
      base::sweep(as.matrix(lhs), MARGIN = margin, STATS = v, FUN = op)
    },
    resident_has = function(key) exists(key, envir = resident, inherits = FALSE),
    resident_store = function(key, x) {
      bump("resident_store")
      assign(key, as.matrix(x), envir = resident)
      invisible(key)
    },
    resident_drop = function(key) {
      bump("resident_drop")
      if (exists(key, envir = resident, inherits = FALSE)) {
        rm(list = key, envir = resident)
      }
      invisible(key)
    },
    resident_materialize = function(key) {
      bump("resident_materialize")
      get(key, envir = resident, inherits = FALSE)
    },
    rowSums_resident = function(x_key, na.rm = FALSE, dims = 1L) {
      bump("rowSums_resident")
      base::rowSums(get(x_key, envir = resident, inherits = FALSE), na.rm = na.rm, dims = dims)
    },
    colSums_resident = function(x_key, na.rm = FALSE, dims = 1L) {
      bump("colSums_resident")
      base::colSums(get(x_key, envir = resident, inherits = FALSE), na.rm = na.rm, dims = dims)
    },
    broadcast_ewise_resident = function(lhs_key, v, margin, op, out_key, defer = FALSE) {
      bump("broadcast_ewise_resident")
      value <- base::sweep(
        get(lhs_key, envir = resident, inherits = FALSE),
        MARGIN = margin,
        STATS = v,
        FUN = op
      )
      assign(out_key, value, envir = resident)
      value
    }
  )
}

test_that("sinkhorn returns a doubly stochastic adgeMatrix on the host path", {
  A <- .sinkhorn_positive_matrix(n = 20L, seed = 101L)
  fit <- sinkhorn(A, tol = 1e-10, return_info = TRUE)

  expect_true(inherits(fit$result, "adgeMatrix"))
  expect_identical(fit$backend, "cpu")
  expect_identical(fit$method, "host")
  expect_true(isTRUE(fit$converged))

  mat <- as.matrix(fit$result)
  expect_equal(base::rowSums(mat), rep(1, nrow(mat)), tolerance = 1e-8)
  expect_equal(base::colSums(mat), rep(1, ncol(mat)), tolerance = 1e-8)
  expect_true(all(mat >= 0))
})

test_that("sinkhorn is invariant to positive scalar rescaling", {
  A <- .sinkhorn_positive_matrix(n = 18L, seed = 202L)

  fit_a <- sinkhorn(A, tol = 1e-10)
  fit_b <- sinkhorn(7 * A, tol = 1e-10)

  expect_equal(as.matrix(fit_a), as.matrix(fit_b), tolerance = 1e-10)
})

test_that("sinkhorn validates shape and positivity constraints", {
  expect_error(sinkhorn(matrix(1:6, nrow = 2L)), "square matrix")
  expect_error(sinkhorn(matrix(c(1, -1, 2, 3), nrow = 2L)), "non-negative")
  expect_error(sinkhorn(matrix(c(1, 0, 0, 0), nrow = 2L)), "strictly positive row sums and column sums")
})

test_that("sinkhorn uses resident sweep kernels when the backend supports them", {
  counter <- new.env(parent = emptyenv())
  backend <- .sinkhorn_recording_backend(counter)

  with_registered_backend("sinkhorn_resident_backend", backend, {
    A_host <- .sinkhorn_positive_matrix(n = 12L, seed = 303L)
    A <- adgeMatrix(A_host, preferred_backend = "sinkhorn_resident_backend", precision = "fast")

    fit <- sinkhorn(A, tol = 1e-10, return_info = TRUE)

    expect_identical(fit$backend, "sinkhorn_resident_backend")
    expect_identical(fit$method, "resident")
    expect_true(isTRUE(fit$result@finalizer_env$host_deferred))
    expect_gt(counter$rowSums_resident, 0L)
    expect_gt(counter$colSums_resident, 0L)
    expect_gt(counter$broadcast_ewise_resident, 0L)

    mat <- as.matrix(fit$result)
    expect_equal(base::rowSums(mat), rep(1, nrow(mat)), tolerance = 1e-10)
    expect_equal(base::colSums(mat), rep(1, ncol(mat)), tolerance = 1e-10)
  })
})

test_that("sinkhorn runs on OpenCL fast matrices when the backend is available", {
  spec <- .sinkhorn_backend_spec("opencl")
  skip_if_backend_package_missing(spec)

  register_backend <- .register_optional_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    A_host <- .sinkhorn_positive_matrix(n = 24L, seed = 404L)
    A <- adgeMatrix(A_host, preferred_backend = "opencl", precision = "fast")
    fit <- sinkhorn(A, tol = 1e-5, return_info = TRUE)

    expect_identical(fit$backend, "opencl")
    expect_identical(fit$method, "resident")
    expect_true(isTRUE(fit$result@finalizer_env$host_deferred))

    mat <- as.matrix(fit$result)
    expect_equal(base::rowSums(mat), rep(1, nrow(mat)), tolerance = 5e-5)
    expect_equal(base::colSums(mat), rep(1, ncol(mat)), tolerance = 5e-5)
  })
})

test_that("sinkhorn runs on MLX fast matrices when the backend is available", {
  spec <- .sinkhorn_backend_spec("mlx")
  skip_if_backend_package_missing(spec)
  skip_if_not(isTRUE(backend_package_available(spec)), "mlx backend not available")

  register_backend <- .register_optional_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    A_host <- .sinkhorn_positive_matrix(n = 24L, seed = 505L)
    A <- adgeMatrix(A_host, preferred_backend = "mlx", precision = "fast")
    fit <- sinkhorn(A, tol = 1e-5, return_info = TRUE)

    expect_identical(fit$backend, "mlx")
    expect_identical(fit$method, "resident")
    expect_true(isTRUE(fit$result@finalizer_env$host_deferred))

    mat <- as.matrix(fit$result)
    expect_equal(base::rowSums(mat), rep(1, nrow(mat)), tolerance = 5e-5)
    expect_equal(base::colSums(mat), rep(1, ncol(mat)), tolerance = 5e-5)
  })
})
