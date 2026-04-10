.amatrix_opencl_probe_var <- "AMATRIX_OPENCL_PROBE_GPU"
.amatrix_opencl_state <- new.env(parent = emptyenv())
.amatrix_opencl_sparse_host_resident <- new.env(parent = emptyenv())

.amatrix_opencl_probe_enabled <- function() {
  identical(Sys.getenv(.amatrix_opencl_probe_var, unset = ""), "1")
}

.amatrix_opencl_probe_cache_get <- function() {
  get0("native_available", envir = .amatrix_opencl_state, inherits = FALSE)
}

.amatrix_opencl_probe_cache_set <- function(value) {
  assign("native_available", isTRUE(value), envir = .amatrix_opencl_state)
  invisible(isTRUE(value))
}

.amatrix_opencl_probe_cache_clear <- function() {
  if (exists("native_available", envir = .amatrix_opencl_state, inherits = FALSE)) {
    rm("native_available", envir = .amatrix_opencl_state)
  }
  invisible(NULL)
}

.amatrix_opencl_dense_host <- function(x) {
  x_mat <- as.matrix(x)
  if (!is.double(x_mat)) {
    storage.mode(x_mat) <- "double"
  }
  x_mat
}

.amatrix_opencl_sparse_host <- function(x) {
  methods::as(x, "dgCMatrix")
}

.amatrix_opencl_rhs_width <- function(y) {
  if (is.null(y)) {
    return(NA_integer_)
  }

  dims <- dim(y)
  if (is.null(dims) || length(dims) <= 1L) {
    return(1L)
  }

  as.integer(dims[[2L]])
}

.amatrix_opencl_rhs_arg <- function(rhs) {
  if (is.null(rhs)) {
    return(NULL)
  }
  if (is.character(rhs) && length(rhs) == 1L) {
    return(as.character(rhs))
  }
  if (is.matrix(rhs)) {
    return(.amatrix_opencl_dense_host(rhs))
  }
  if (inherits(rhs, "Matrix") || inherits(rhs, "denseMatrix")) {
    return(.amatrix_opencl_dense_host(rhs))
  }
  if (is.numeric(rhs) && length(rhs) == 1L) {
    return(as.double(rhs))
  }
  rhs_mat <- as.matrix(rhs)
  if (!is.double(rhs_mat)) {
    storage.mode(rhs_mat) <- "double"
  }
  rhs_mat
}

.amatrix_opencl_broadcast_arg <- function(v) {
  if (is.character(v) && length(v) == 1L) {
    return(as.character(v))
  }
  if (!is.numeric(v)) {
    stop("broadcast operand must be a numeric vector or resident key", call. = FALSE)
  }
  as.double(v)
}

.amatrix_opencl_temp_key <- function(prefix = "tmp") {
  counter <- get0("temp_key_counter", envir = .amatrix_opencl_state, inherits = FALSE, ifnotfound = 0L) + 1L
  assign("temp_key_counter", counter, envir = .amatrix_opencl_state)
  sprintf("%s:%d", prefix, counter)
}

.amatrix_opencl_factor_gpu_enabled <- function() {
  isTRUE(getOption("amatrix.opencl.factor_gpu", FALSE))
}

.amatrix_opencl_experimental_qr_solve_enabled <- function() {
  isTRUE(getOption("amatrix.opencl.experimental_qr_solve", FALSE))
}

.amatrix_opencl_is_square_matrix <- function(x) {
  dims <- dim(x)
  !is.null(dims) && length(dims) == 2L && identical(dims[[1L]], dims[[2L]])
}

.amatrix_opencl_is_symmetric_matrix <- function(x, tol = sqrt(.Machine$double.eps)) {
  .amatrix_opencl_is_square_matrix(x) &&
    isTRUE(isSymmetric(.amatrix_opencl_dense_host(x), tol = tol))
}

.amatrix_opencl_can_try_device_spd_solve <- function(a) {
  .amatrix_opencl_device_linalg_available() &&
    .amatrix_opencl_is_symmetric_matrix(a)
}

.amatrix_opencl_can_try_device_qr_solve <- function(a, b = NULL) {
  rhs_width <- .amatrix_opencl_rhs_width(b)

  .amatrix_opencl_experimental_qr_solve_enabled() &&
    !is.null(b) &&
    !is.na(rhs_width) &&
    .amatrix_opencl_device_linalg_available() &&
    .amatrix_opencl_is_square_matrix(a) &&
    max(dim(a)) >= getOption("amatrix.opencl.solve_qr_min_dim", 1536L) &&
    rhs_width >= getOption("amatrix.opencl.solve_qr_min_rhs", 64L)
}

.amatrix_opencl_dense_product_supported <- function(op, x, y = NULL) {
  threshold <- as.integer(getOption("amatrix.opencl.matmul_min_dim", 128L))
  dims <- dim(x)

  if (identical(op, "matmul") && !is.null(y) && !is.null(dim(y))) {
    dims <- c(dims, dim(y))
  }

  min(as.integer(dims)) >= threshold
}

.amatrix_opencl_sparse_product_supported <- function(op, x, y = NULL) {
  if (!inherits(x, "adgCMatrix")) {
    return(FALSE)
  }

  if (!(op %in% c("matmul", "crossprod", "tcrossprod"))) {
    return(FALSE)
  }

  if (op %in% c("crossprod", "tcrossprod") && is.null(y)) {
    return(FALSE)
  }

  rhs_width <- .amatrix_opencl_rhs_width(y)
  min_nnz <- if (!is.na(rhs_width) && rhs_width <= 1L) {
    getOption("amatrix.opencl.spmv_min_nnz", Inf)
  } else {
    getOption("amatrix.opencl.spmm_min_nnz", Inf)
  }

  if (!is.finite(min_nnz)) {
    return(FALSE)
  }

  length(x@x) >= as.integer(min_nnz)
}

.amatrix_opencl_qr_use_cholqr <- function(x_mat) {
  if (!isTRUE(.amatrix_opencl_device_linalg_available())) {
    return(FALSE)
  }

  n <- nrow(x_mat)
  p <- ncol(x_mat)
  n >= getOption("amatrix.opencl.qr_min_n", 4096L) &&
    p >= getOption("amatrix.opencl.qr_min_p", 16L) &&
    p <= getOption("amatrix.opencl.qr_max_p", 256L) &&
    n >= (4L * p)
}

.amatrix_opencl_device_linalg_available <- function(force = FALSE) {
  .amatrix_opencl_factor_gpu_enabled() &&
    isTRUE(amatrix_opencl_native_available(force = force)) &&
    isTRUE(amatrix_opencl_bridge_info()$clblast)
}

.amatrix_opencl_drop_if_present <- function(key) {
  if (isTRUE(amatrix_opencl_resident_has(key))) {
    amatrix_opencl_resident_drop(key)
  }
  invisible(NULL)
}

.amatrix_opencl_temp_dense_result <- function(x, prefix, worker) {
  in_key <- .amatrix_opencl_temp_key(paste0(prefix, "-in"))
  out_key <- .amatrix_opencl_temp_key(paste0(prefix, "-out"))
  on.exit(.amatrix_opencl_drop_if_present(in_key), add = TRUE)
  on.exit(.amatrix_opencl_drop_if_present(out_key), add = TRUE)

  amatrix_opencl_resident_store(in_key, .amatrix_opencl_dense_host(x))
  worker(in_key, out_key)
  amatrix_opencl_resident_materialize(out_key)
}

.amatrix_opencl_temp_solve_result <- function(a, b = NULL, prefix = "solve", worker) {
  a_key <- .amatrix_opencl_temp_key(paste0(prefix, "-a"))
  out_key <- .amatrix_opencl_temp_key(paste0(prefix, "-out"))
  b_key <- if (is.null(b)) NULL else .amatrix_opencl_temp_key(paste0(prefix, "-b"))

  on.exit(.amatrix_opencl_drop_if_present(a_key), add = TRUE)
  on.exit(.amatrix_opencl_drop_if_present(out_key), add = TRUE)
  if (!is.null(b_key)) {
    on.exit(.amatrix_opencl_drop_if_present(b_key), add = TRUE)
  }

  amatrix_opencl_resident_store(a_key, .amatrix_opencl_dense_host(a))
  if (!is.null(b_key)) {
    amatrix_opencl_resident_store(b_key, .amatrix_opencl_dense_host(b))
  }

  worker(a_key, b_key, out_key)
  amatrix_opencl_resident_materialize(out_key)
}

amatrix_opencl_capabilities <- function() {
  c("matmul", "crossprod", "tcrossprod", "ewise",
    "broadcast_ewise", "rowSums", "colSums", "solve", "chol", "qr", "svd", "eigen", "covariance")
}

amatrix_opencl_features <- function() {
  c("dense_f32", "resident_dense", "custom_ops", "solve", "chol", "qr", "svd", "eigen_sym", "covariance", "sparse_spmm")
}

amatrix_opencl_precision_modes <- function() {
  "fast"
}

amatrix_opencl_native_available <- function(force = FALSE) {
  if (!force) {
    cached <- .amatrix_opencl_probe_cache_get()
    if (!is.null(cached)) {
      return(cached)
    }
  }

  if (!.amatrix_opencl_probe_enabled()) {
    return(FALSE)
  }

  available <- isTRUE(.Call("amatrix_opencl_native_available_bridge", PACKAGE = "amatrix.opencl"))
  .amatrix_opencl_probe_cache_set(available)
}

amatrix_opencl_enable_probe <- function(register = TRUE) {
  Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = "1")
  options(amatrix.enable_opencl = TRUE)
  .amatrix_opencl_probe_cache_clear()

  available <- amatrix_opencl_native_available(force = TRUE)
  if (isTRUE(register)) {
    try(amatrix_opencl_register(overwrite = TRUE), silent = TRUE)
  }

  invisible(available)
}

amatrix_opencl_is_available <- function() {
  isTRUE(getOption("amatrix.opencl.available", FALSE)) || isTRUE(amatrix_opencl_native_available())
}

amatrix_opencl_bridge_info <- function() {
  info <- .Call("amatrix_opencl_bridge_info_bridge", PACKAGE = "amatrix.opencl")
  info$available <- amatrix_opencl_is_available()
  info$capabilities <- amatrix_opencl_capabilities()
  info
}

amatrix_opencl_diagnostics <- function() {
  diag <- .Call("amatrix_opencl_diagnostics_bridge", PACKAGE = "amatrix.opencl")
  diag$available <- amatrix_opencl_is_available()
  diag
}

amatrix_opencl_matmul <- function(x, y) {
  if (inherits(x, "dgCMatrix")) {
    return(as.matrix(.amatrix_opencl_sparse_host(x) %*% .amatrix_opencl_dense_host(y)))
  }
  .Call(
    "amatrix_opencl_matmul_bridge",
    .amatrix_opencl_dense_host(x),
    .amatrix_opencl_dense_host(y),
    PACKAGE = "amatrix.opencl"
  )
}

amatrix_opencl_crossprod <- function(x, y = NULL) {
  if (inherits(x, "dgCMatrix")) {
    if (is.null(y)) {
      return(as.matrix(Matrix::crossprod(.amatrix_opencl_sparse_host(x))))
    }
    return(as.matrix(Matrix::crossprod(.amatrix_opencl_sparse_host(x), .amatrix_opencl_dense_host(y))))
  }
  rhs <- if (is.null(y)) NULL else .amatrix_opencl_dense_host(y)
  .Call(
    "amatrix_opencl_crossprod_bridge",
    .amatrix_opencl_dense_host(x),
    rhs,
    PACKAGE = "amatrix.opencl"
  )
}

amatrix_opencl_tcrossprod <- function(x, y = NULL) {
  if (inherits(x, "dgCMatrix")) {
    if (is.null(y)) {
      return(as.matrix(Matrix::tcrossprod(.amatrix_opencl_sparse_host(x))))
    }
    return(as.matrix(Matrix::tcrossprod(.amatrix_opencl_sparse_host(x), .amatrix_opencl_dense_host(y))))
  }
  rhs <- if (is.null(y)) NULL else .amatrix_opencl_dense_host(y)
  .Call(
    "amatrix_opencl_tcrossprod_bridge",
    .amatrix_opencl_dense_host(x),
    rhs,
    PACKAGE = "amatrix.opencl"
  )
}

amatrix_opencl_ewise <- function(lhs, rhs = NULL, op) {
  .Call(
    "amatrix_opencl_ewise_bridge",
    .amatrix_opencl_dense_host(lhs),
    .amatrix_opencl_rhs_arg(rhs),
    as.character(op),
    PACKAGE = "amatrix.opencl"
  )
}

amatrix_opencl_broadcast_ewise <- function(lhs, v, margin, op) {
  .Call(
    "amatrix_opencl_broadcast_ewise_bridge",
    .amatrix_opencl_dense_host(lhs),
    as.double(v),
    as.integer(margin),
    as.character(op),
    PACKAGE = "amatrix.opencl"
  )
}

amatrix_opencl_axis_sums <- function(x, axis) {
  .Call(
    "amatrix_opencl_sum_axis_bridge",
    .amatrix_opencl_dense_host(x),
    as.integer(axis),
    PACKAGE = "amatrix.opencl"
  )
}

.amatrix_opencl_chol_host <- function(x) {
  x_mat <- .amatrix_opencl_dense_host(x)
  if (nrow(x_mat) != ncol(x_mat)) {
    stop("x must be a square matrix", call. = FALSE)
  }
  result <- base::chol(x_mat)
  result[lower.tri(result)] <- 0
  result
}

.amatrix_opencl_solve_host <- function(a, b = NULL) {
  a_mat <- .amatrix_opencl_dense_host(a)
  if (nrow(a_mat) != ncol(a_mat)) {
    stop("a must be a square matrix", call. = FALSE)
  }

  if (is.null(b)) {
    return(base::solve(a_mat))
  }

  b_mat <- .amatrix_opencl_dense_host(b)
  base::solve(a_mat, b_mat)
}

.amatrix_opencl_solve_triangular_host <- function(R, B, lower = FALSE, transpose = FALSE) {
  r_mat <- .amatrix_opencl_dense_host(R)
  b_mat <- .amatrix_opencl_dense_host(B)

  if (nrow(r_mat) != ncol(r_mat)) {
    stop("triangular factor must be square", call. = FALSE)
  }
  if (nrow(b_mat) != nrow(r_mat)) {
    stop("triangular solve rhs has incompatible dimensions", call. = FALSE)
  }

  if (isTRUE(transpose) && isTRUE(lower)) {
    return(backsolve(t(r_mat), b_mat))
  }
  if (isTRUE(transpose) && !isTRUE(lower)) {
    return(forwardsolve(t(r_mat), b_mat))
  }
  if (isTRUE(lower)) {
    return(forwardsolve(r_mat, b_mat))
  }
  backsolve(r_mat, b_mat)
}

amatrix_opencl_chol <- function(x) {
  if (.amatrix_opencl_device_linalg_available()) {
    result <- tryCatch(
      .amatrix_opencl_temp_dense_result(
        x,
        prefix = "chol",
        worker = function(in_key, out_key) {
          amatrix_opencl_chol_resident(in_key, out_key, defer = TRUE)
        }
      ),
      error = function(e) NULL
    )
    if (!is.null(result)) {
      return(result)
    }
  }
  .amatrix_opencl_chol_host(x)
}

amatrix_opencl_solve <- function(a, b = NULL) {
  if (.amatrix_opencl_can_try_device_spd_solve(a)) {
    result <- tryCatch(
      .amatrix_opencl_temp_solve_result(
        a,
        b = b,
        prefix = "solve",
        worker = function(a_key, b_key, out_key) {
          amatrix_opencl_solve_resident(a_key, b_key = b_key, out_key = out_key, defer = TRUE)
        }
      ),
      error = function(e) NULL
    )
    if (!is.null(result)) {
      return(result)
    }
  }
  if (.amatrix_opencl_can_try_device_qr_solve(a, b = b)) {
    result <- tryCatch(
      .amatrix_opencl_qr_solve_rhs(a, b),
      error = function(e) NULL
    )
    if (!is.null(result)) {
      return(result)
    }
  }
  .amatrix_opencl_solve_host(a, b = b)
}

amatrix_opencl_qr <- function(x, ...) {
  x_host <- .amatrix_opencl_dense_host(x)
  extra_args <- list(...)

  if (length(extra_args) == 0L && .amatrix_opencl_qr_use_cholqr(x_host)) {
    qr_fast <- tryCatch(
      amatrix_opencl_qr_cholqr(x_host),
      error = function(e) NULL
    )
    if (!is.null(qr_fast)) {
      return(qr_fast)
    }
  }

  qr_host <- base::qr(x_host, ...)
  list(
    qr = qr_host,
    factor = qr_host,
    factor_source = "native",
    backend_ops = "opencl"
  )
}

amatrix_opencl_qr_cholqr <- function(x) {
  x_mat <- .amatrix_opencl_dense_host(x)
  n <- nrow(x_mat)
  p <- ncol(x_mat)
  diag_tol <- sqrt(.Machine$double.eps)
  success <- FALSE

  x_key <- .amatrix_opencl_temp_key("cholqr-x")
  gram1_key <- .amatrix_opencl_temp_key("cholqr-gram1")
  r1_key <- .amatrix_opencl_temp_key("cholqr-r1")
  inv1_key <- .amatrix_opencl_temp_key("cholqr-inv1")
  q1_key <- .amatrix_opencl_temp_key("cholqr-q1")
  gram2_key <- .amatrix_opencl_temp_key("cholqr-gram2")
  r2_key <- .amatrix_opencl_temp_key("cholqr-r2")
  inv2_key <- .amatrix_opencl_temp_key("cholqr-inv2")
  q_key <- .amatrix_opencl_temp_key("cholqr-q")
  cleanup_keys <- c(x_key, gram1_key, r1_key, inv1_key, q1_key, gram2_key, r2_key, inv2_key)

  on.exit({
    for (key in cleanup_keys) {
      .amatrix_opencl_drop_if_present(key)
    }
    if (!success) {
      .amatrix_opencl_drop_if_present(q_key)
    }
  }, add = TRUE)

  amatrix_opencl_resident_store(x_key, x_mat)
  amatrix_opencl_crossprod_resident(x_key, out_key = gram1_key, defer = TRUE)
  gram1 <- amatrix_opencl_resident_materialize(gram1_key)
  gram1_sv <- tryCatch(base::svd(gram1, nu = 0L, nv = 0L)$d, error = function(e) numeric())
  gram_tol <- getOption("amatrix.opencl.qr_gram_rank_tol", 1e-8)
  if (length(gram1_sv) != p || any(!is.finite(gram1_sv)) || min(gram1_sv) <= gram_tol * max(1, max(gram1_sv))) {
    stop("OpenCL CholeskyQR2 rejected rank-deficient input", call. = FALSE)
  }
  amatrix_opencl_chol_resident(gram1_key, out_key = r1_key, defer = TRUE)
  r1 <- amatrix_opencl_resident_materialize(r1_key)
  inv_r1 <- backsolve(r1, diag(p))
  amatrix_opencl_resident_store(inv1_key, inv_r1)
  amatrix_opencl_matmul_resident(x_key, inv1_key, q1_key, defer = TRUE)

  amatrix_opencl_crossprod_resident(q1_key, out_key = gram2_key, defer = TRUE)
  amatrix_opencl_chol_resident(gram2_key, out_key = r2_key, defer = TRUE)
  r2 <- amatrix_opencl_resident_materialize(r2_key)
  inv_r2 <- backsolve(r2, diag(p))
  amatrix_opencl_resident_store(inv2_key, inv_r2)
  amatrix_opencl_matmul_resident(q1_key, inv2_key, q_key, defer = TRUE)

  r <- r2 %*% r1
  diag_r <- abs(diag(r))
  sv_tol <- getOption("amatrix.opencl.qr_rank_tol", 1e-5)
  if (length(diag_r) != p || any(!is.finite(diag_r)) || any(diag_r <= diag_tol * max(1, max(diag_r)))) {
    stop("OpenCL CholeskyQR2 failed numerical rank check", call. = FALSE)
  }
  sv <- tryCatch(base::svd(r, nu = 0L, nv = 0L)$d, error = function(e) numeric())
  if (length(sv) != p || any(!is.finite(sv)) || min(sv) <= sv_tol * max(1, max(sv))) {
    stop("OpenCL CholeskyQR2 failed full-rank verification", call. = FALSE)
  }

  success <- TRUE
  list(
    representation = "explicit_qr",
    q_key = q_key,
    r = r,
    rank = p,
    pivot = seq_len(p),
    factor = NULL,
    factor_source = "cholqr2",
    backend_ops = "opencl"
  )
}

amatrix_opencl_qr_qty_key <- function(q_key, y) {
  .Call(
    "amatrix_opencl_crossprod_resident_host_bridge",
    as.character(q_key),
    .amatrix_opencl_dense_host(y),
    PACKAGE = "amatrix.opencl"
  )
}

amatrix_opencl_qr_qy_key <- function(q_key, y) {
  .Call(
    "amatrix_opencl_matmul_resident_host_bridge",
    as.character(q_key),
    .amatrix_opencl_dense_host(y),
    PACKAGE = "amatrix.opencl"
  )
}

amatrix_opencl_matmul_resident_host_into <- function(x_key, y, out_key) {
  invisible(.Call(
    "amatrix_opencl_matmul_resident_host_into_bridge",
    as.character(x_key),
    .amatrix_opencl_dense_host(y),
    as.character(out_key),
    PACKAGE = "amatrix.opencl"
  ))
}

amatrix_opencl_crossprod_resident_host_into <- function(x_key, y, out_key) {
  invisible(.Call(
    "amatrix_opencl_crossprod_resident_host_into_bridge",
    as.character(x_key),
    .amatrix_opencl_dense_host(y),
    as.character(out_key),
    PACKAGE = "amatrix.opencl"
  ))
}

amatrix_opencl_qr_fitted_key <- function(q_key, y) {
  qty_key <- .amatrix_opencl_temp_key("qr-qty")
  fitted_key <- .amatrix_opencl_temp_key("qr-fitted")

  on.exit(.amatrix_opencl_drop_if_present(qty_key), add = TRUE)
  on.exit(.amatrix_opencl_drop_if_present(fitted_key), add = TRUE)

  amatrix_opencl_crossprod_resident_host_into(q_key, y, qty_key)
  amatrix_opencl_matmul_resident(q_key, qty_key, fitted_key, defer = TRUE)
  amatrix_opencl_resident_materialize(fitted_key)
}

.amatrix_opencl_wrap_qr <- function(qr_value, x_host) {
  amatrix_ns <- asNamespace("amatrix")
  wrap_qr <- get(".amatrix_wrap_qr", envir = amatrix_ns, inherits = FALSE)
  wrap_qr(
    qr_value,
    amatrix::adgeMatrix(x_host, preferred_backend = "opencl", precision = "fast"),
    method = "fast"
  )
}

.amatrix_opencl_qr_solve_rhs <- function(a, b, tol = 1e-07) {
  amatrix_ns <- asNamespace("amatrix")
  qr_solve_value <- get(".amatrix_qr_solve_value", envir = amatrix_ns, inherits = FALSE)
  qr_fit <- .amatrix_opencl_wrap_qr(amatrix_opencl_qr(a), a)
  qr_solve_value(qr_fit, b = .amatrix_opencl_dense_host(b), tol = tol)
}

amatrix_opencl_svd <- function(x, nu, nv, LINPACK = FALSE, ...) {
  if (isTRUE(LINPACK)) {
    stop("LINPACK is not supported", call. = FALSE)
  }

  base::svd(.amatrix_opencl_dense_host(x), nu = nu, nv = nv, ...)
}

amatrix_opencl_eigen <- function(x, symmetric, only.values = FALSE, EISPACK = FALSE) {
  x_host <- .amatrix_opencl_dense_host(x)
  if (nrow(x_host) != ncol(x_host)) {
    stop("x must be a square matrix", call. = FALSE)
  }

  if (!isTRUE(symmetric)) {
    return(base::eigen(x_host, symmetric = FALSE, only.values = only.values, EISPACK = EISPACK))
  }

  base::eigen(x_host, symmetric = TRUE, only.values = only.values, EISPACK = EISPACK)
}

amatrix_opencl_solve_triangular_resident <- function(factor_key, rhs_key, out_key, lower = FALSE, transpose = FALSE, defer = FALSE) {
  result <- tryCatch({
    .Call(
      "amatrix_opencl_solve_triangular_resident_bridge",
      as.character(factor_key),
      as.character(rhs_key),
      as.character(out_key),
      as.logical(lower),
      as.logical(transpose),
      PACKAGE = "amatrix.opencl"
    )
    if (defer) NULL else amatrix_opencl_resident_materialize(out_key)
  }, error = function(e) {
    factor_host <- amatrix_opencl_resident_materialize(factor_key)
    rhs_host <- amatrix_opencl_resident_materialize(rhs_key)
    host_result <- .amatrix_opencl_solve_triangular_host(
      factor_host,
      rhs_host,
      lower = lower,
      transpose = transpose
    )
    amatrix_opencl_resident_store(out_key, host_result)
    if (defer) NULL else host_result
  })
  result
}

amatrix_opencl_solve_triangular_factor <- function(R, B, lower = FALSE, transpose = FALSE) {
  if (.amatrix_opencl_device_linalg_available()) {
    factor_key <- .amatrix_opencl_temp_key("tri-factor")
    rhs_key <- .amatrix_opencl_temp_key("tri-rhs")
    out_key <- .amatrix_opencl_temp_key("tri-out")
    on.exit(.amatrix_opencl_drop_if_present(factor_key), add = TRUE)
    on.exit(.amatrix_opencl_drop_if_present(rhs_key), add = TRUE)
    on.exit(.amatrix_opencl_drop_if_present(out_key), add = TRUE)

    result <- tryCatch({
      amatrix_opencl_resident_store(factor_key, .amatrix_opencl_dense_host(R))
      amatrix_opencl_resident_store(rhs_key, .amatrix_opencl_dense_host(B))
      amatrix_opencl_solve_triangular_resident(
        factor_key,
        rhs_key,
        out_key,
        lower = lower,
        transpose = transpose,
        defer = TRUE
      )
      amatrix_opencl_resident_materialize(out_key)
    }, error = function(e) NULL)

    if (!is.null(result)) {
      return(result)
    }
  }

  .amatrix_opencl_solve_triangular_host(R, B, lower = lower, transpose = transpose)
}

amatrix_opencl_chol_solve_factor <- function(R, B) {
  z <- amatrix_opencl_solve_triangular_factor(R, B, lower = FALSE, transpose = TRUE)
  amatrix_opencl_solve_triangular_factor(R, z, lower = FALSE, transpose = FALSE)
}

amatrix_opencl_covariance <- function(x, center = TRUE, denom = NULL) {
  x_mat <- .amatrix_opencl_dense_host(x)
  n <- nrow(x_mat)
  d <- if (is.null(denom)) n - 1L else as.double(denom)
  if (!is.finite(d) || length(d) != 1L || d <= 0) {
    stop("denom must be a single positive value", call. = FALSE)
  }

  if (isTRUE(center)) {
    x_mat <- sweep(x_mat, 2L, colMeans(x_mat), FUN = "-")
  }

  amatrix_opencl_crossprod(x_mat) / d
}

amatrix_opencl_resident_store <- function(key, x) {
  invisible(.Call(
    "amatrix_opencl_resident_store_bridge",
    as.character(key),
    .amatrix_opencl_dense_host(x),
    PACKAGE = "amatrix.opencl"
  ))
}

amatrix_opencl_resident_has <- function(key) {
  isTRUE(.Call(
    "amatrix_opencl_resident_has_bridge",
    as.character(key),
    PACKAGE = "amatrix.opencl"
  ))
}

amatrix_opencl_resident_drop <- function(key) {
  invisible(.Call(
    "amatrix_opencl_resident_drop_bridge",
    as.character(key),
    PACKAGE = "amatrix.opencl"
  ))
}

amatrix_opencl_resident_materialize <- function(key) {
  .Call(
    "amatrix_opencl_resident_materialize_bridge",
    as.character(key),
    PACKAGE = "amatrix.opencl"
  )
}

amatrix_opencl_sparse_resident_store <- function(key, x_sp) {
  assign(as.character(key), .amatrix_opencl_sparse_host(x_sp), envir = .amatrix_opencl_sparse_host_resident)
  invisible(TRUE)
}

amatrix_opencl_sparse_resident_has <- function(key) {
  exists(as.character(key), envir = .amatrix_opencl_sparse_host_resident, inherits = FALSE)
}

amatrix_opencl_sparse_resident_drop <- function(key) {
  key <- as.character(key)
  if (exists(key, envir = .amatrix_opencl_sparse_host_resident, inherits = FALSE)) {
    rm(list = key, envir = .amatrix_opencl_sparse_host_resident)
  }
  invisible(TRUE)
}

amatrix_opencl_spmm_resident <- function(sp_key, B, trans_lhs = FALSE) {
  sp_key <- as.character(sp_key)
  if (!exists(sp_key, envir = .amatrix_opencl_sparse_host_resident, inherits = FALSE)) {
    stop(sprintf("unknown opencl sparse resident key '%s'", sp_key), call. = FALSE)
  }

  sparse_x <- get(sp_key, envir = .amatrix_opencl_sparse_host_resident, inherits = FALSE)
  rhs <- .amatrix_opencl_dense_host(B)
  if (isTRUE(trans_lhs)) {
    return(as.matrix(Matrix::crossprod(sparse_x, rhs)))
  }
  as.matrix(sparse_x %*% rhs)
}

amatrix_opencl_spmm_resident_key <- function(sp_key, y_key, out_key, trans_lhs = FALSE, defer = FALSE) {
  host_result <- amatrix_opencl_spmm_resident(sp_key, amatrix_opencl_resident_materialize(y_key), trans_lhs = trans_lhs)
  amatrix_opencl_resident_store(out_key, host_result)
  if (isTRUE(defer)) {
    return(NULL)
  }
  host_result
}

amatrix_opencl_solve_resident <- function(a_key, b_key = NULL, out_key, defer = FALSE) {
  result <- tryCatch({
    .Call(
      "amatrix_opencl_solve_resident_bridge",
      as.character(a_key),
      if (is.null(b_key)) NULL else as.character(b_key),
      as.character(out_key),
      PACKAGE = "amatrix.opencl"
    )
    if (defer) NULL else amatrix_opencl_resident_materialize(out_key)
  }, error = function(e) {
    a_host <- amatrix_opencl_resident_materialize(a_key)
    host_result <- if (is.null(b_key)) {
      .amatrix_opencl_solve_host(a_host)
    } else {
      b_host <- amatrix_opencl_resident_materialize(b_key)
      .amatrix_opencl_solve_host(a_host, b = b_host)
    }
    amatrix_opencl_resident_store(out_key, host_result)
    if (defer) NULL else host_result
  })
  result
}

amatrix_opencl_chol_resident <- function(x_key, out_key, defer = FALSE) {
  result <- tryCatch({
    .Call(
      "amatrix_opencl_chol_resident_bridge",
      as.character(x_key),
      as.character(out_key),
      PACKAGE = "amatrix.opencl"
    )
    if (defer) NULL else amatrix_opencl_resident_materialize(out_key)
  }, error = function(e) {
    host_result <- .amatrix_opencl_chol_host(amatrix_opencl_resident_materialize(x_key))
    amatrix_opencl_resident_store(out_key, host_result)
    if (defer) NULL else host_result
  })
  result
}

amatrix_opencl_qr_Q_resident <- function(x_key, q_out_key) {
  x_host <- amatrix_opencl_resident_materialize(x_key)
  q_host <- qr.Q(base::qr(x_host), complete = FALSE)
  amatrix_opencl_resident_store(q_out_key, q_host)
  invisible(q_host)
}

amatrix_opencl_matmul_resident <- function(x_key, y_key, out_key, defer = FALSE) {
  .Call(
    "amatrix_opencl_matmul_resident_bridge",
    as.character(x_key),
    as.character(y_key),
    as.character(out_key),
    PACKAGE = "amatrix.opencl"
  )
  if (defer) NULL else amatrix_opencl_resident_materialize(out_key)
}

amatrix_opencl_crossprod_resident <- function(x_key, y_key = NULL, out_key, defer = FALSE) {
  rhs_key <- if (is.null(y_key)) NULL else as.character(y_key)
  .Call(
    "amatrix_opencl_crossprod_resident_bridge",
    as.character(x_key),
    rhs_key,
    as.character(out_key),
    PACKAGE = "amatrix.opencl"
  )
  if (defer) NULL else amatrix_opencl_resident_materialize(out_key)
}

amatrix_opencl_tcrossprod_resident <- function(x_key, y_key = NULL, out_key, defer = FALSE) {
  rhs_key <- if (is.null(y_key)) NULL else as.character(y_key)
  .Call(
    "amatrix_opencl_tcrossprod_resident_bridge",
    as.character(x_key),
    rhs_key,
    as.character(out_key),
    PACKAGE = "amatrix.opencl"
  )
  if (defer) NULL else amatrix_opencl_resident_materialize(out_key)
}

amatrix_opencl_matmul_resident_host <- function(x_key, y) {
  .Call(
    "amatrix_opencl_matmul_resident_host_bridge",
    as.character(x_key),
    .amatrix_opencl_dense_host(y),
    PACKAGE = "amatrix.opencl"
  )
}

amatrix_opencl_crossprod_resident_host <- function(x_key, y) {
  .Call(
    "amatrix_opencl_crossprod_resident_host_bridge",
    as.character(x_key),
    .amatrix_opencl_dense_host(y),
    PACKAGE = "amatrix.opencl"
  )
}

amatrix_opencl_ewise_resident <- function(lhs_key, rhs, op, out_key, defer = FALSE) {
  .Call(
    "amatrix_opencl_ewise_resident_bridge",
    as.character(lhs_key),
    .amatrix_opencl_rhs_arg(rhs),
    as.character(op),
    as.character(out_key),
    PACKAGE = "amatrix.opencl"
  )
  if (defer) NULL else amatrix_opencl_resident_materialize(out_key)
}

amatrix_opencl_broadcast_ewise_resident <- function(lhs_key, v, margin, op, out_key, defer = FALSE) {
  .Call(
    "amatrix_opencl_broadcast_ewise_resident_bridge",
    as.character(lhs_key),
    .amatrix_opencl_broadcast_arg(v),
    as.integer(margin),
    as.character(op),
    as.character(out_key),
    PACKAGE = "amatrix.opencl"
  )
  if (defer) NULL else amatrix_opencl_resident_materialize(out_key)
}

amatrix_opencl_broadcast_ewise_resident_inplace <- function(lhs_key, v, margin, op) {
  .Call(
    "amatrix_opencl_broadcast_ewise_resident_inplace_bridge",
    as.character(lhs_key),
    .amatrix_opencl_broadcast_arg(v),
    as.integer(margin),
    as.character(op),
    PACKAGE = "amatrix.opencl"
  )
  invisible(as.character(lhs_key))
}

amatrix_opencl_sum_axis_resident_key <- function(x_key, axis, out_key) {
  .Call(
    "amatrix_opencl_sum_axis_resident_key_bridge",
    as.character(x_key),
    as.integer(axis),
    as.character(out_key),
    PACKAGE = "amatrix.opencl"
  )
  invisible(as.character(out_key))
}

amatrix_opencl_rowSums_resident_key <- function(x_key, out_key, na.rm = FALSE, dims = 1L) {
  if (isTRUE(na.rm) || !identical(dims, 1L)) {
    values <- base::rowSums(amatrix_opencl_resident_materialize(x_key), na.rm = na.rm, dims = dims)
    amatrix_opencl_resident_store(out_key, matrix(as.double(values), ncol = 1L))
    return(invisible(as.character(out_key)))
  }
  amatrix_opencl_sum_axis_resident_key(x_key, axis = 0L, out_key = out_key)
}

amatrix_opencl_colSums_resident_key <- function(x_key, out_key, na.rm = FALSE, dims = 1L) {
  if (isTRUE(na.rm) || !identical(dims, 1L)) {
    values <- base::colSums(amatrix_opencl_resident_materialize(x_key), na.rm = na.rm, dims = dims)
    amatrix_opencl_resident_store(out_key, matrix(as.double(values), ncol = 1L))
    return(invisible(as.character(out_key)))
  }
  amatrix_opencl_sum_axis_resident_key(x_key, axis = 1L, out_key = out_key)
}

amatrix_opencl_rowSums_resident <- function(x_key, na.rm = FALSE, dims = 1L) {
  if (isTRUE(na.rm) || !identical(dims, 1L)) {
    return(base::rowSums(amatrix_opencl_resident_materialize(x_key), na.rm = na.rm, dims = dims))
  }
  .Call(
    "amatrix_opencl_sum_axis_resident_bridge",
    as.character(x_key),
    as.integer(0L),
    PACKAGE = "amatrix.opencl"
  )
}

amatrix_opencl_colSums_resident <- function(x_key, na.rm = FALSE, dims = 1L) {
  if (isTRUE(na.rm) || !identical(dims, 1L)) {
    return(base::colSums(amatrix_opencl_resident_materialize(x_key), na.rm = na.rm, dims = dims))
  }
  .Call(
    "amatrix_opencl_sum_axis_resident_bridge",
    as.character(x_key),
    as.integer(1L),
    PACKAGE = "amatrix.opencl"
  )
}

amatrix_opencl_backend <- function() {
  caps <- amatrix_opencl_capabilities()
  list(
    capabilities = function() caps,
    features = function() amatrix_opencl_features(),
    precision_modes = function() amatrix_opencl_precision_modes(),
    available = function() amatrix_opencl_is_available(),
    supports_resident = function(op, x, y = NULL) {
      if (!isTRUE(amatrix_opencl_is_available())) {
        return(FALSE)
      }
      if (!inherits(x, "aMatrix")) {
        return(FALSE)
      }
      if (!(x@precision %in% amatrix_opencl_precision_modes())) {
        return(FALSE)
      }
      if (!op %in% caps) {
        return(FALSE)
      }
      if (inherits(x, "adgCMatrix")) {
        return(.amatrix_opencl_sparse_product_supported(op, x, y = y))
      }
      if (op %in% c("solve", "chol")) {
        return(
          .amatrix_opencl_factor_gpu_enabled() &&
            .amatrix_opencl_is_square_matrix(x) &&
            max(dim(x)) >= getOption("amatrix.opencl.factor_min_dim", 1024L) &&
            (identical(op, "chol") || .amatrix_opencl_is_symmetric_matrix(x))
        )
      }
      TRUE
    },
    supports = function(op, x, y = NULL) {
      if (!isTRUE(amatrix_opencl_is_available())) {
        return(FALSE)
      }
      if (!op %in% caps) {
        return(FALSE)
      }
      if (inherits(x, "adgCMatrix")) {
        if (!(x@precision %in% amatrix_opencl_precision_modes())) {
          return(FALSE)
        }
        return(.amatrix_opencl_sparse_product_supported(op, x, y = y))
      }
      if (!inherits(x, "adgeMatrix")) {
        return(FALSE)
      }
      if (!(x@precision %in% amatrix_opencl_precision_modes())) {
        return(FALSE)
      }

      if (op %in% c("matmul", "crossprod", "tcrossprod")) {
        return(.amatrix_opencl_dense_product_supported(op, x, y = y))
      }

      if (op %in% c("solve", "chol")) {
        if (!.amatrix_opencl_is_square_matrix(x)) {
          return(FALSE)
        }
        if (identical(op, "chol")) {
          return(TRUE)
        }
        return(
          .amatrix_opencl_is_symmetric_matrix(x) ||
            .amatrix_opencl_can_try_device_qr_solve(x, b = y)
        )
      }

      if (identical(op, "eigen")) {
        return(.amatrix_opencl_is_square_matrix(x))
      }

      if (identical(op, "svd")) {
        return(TRUE)
      }

      TRUE
    },
    matmul = function(x, y) {
      amatrix_opencl_matmul(x, y)
    },
    crossprod = function(x, y = NULL, ...) {
      amatrix_opencl_crossprod(x, y = y)
    },
    tcrossprod = function(x, y = NULL, ...) {
      amatrix_opencl_tcrossprod(x, y = y)
    },
    ewise = function(x, lhs, rhs = NULL, op, ...) {
      amatrix_opencl_ewise(lhs, rhs = rhs, op = op)
    },
    broadcast_ewise = function(x, lhs, v, margin, op, ...) {
      amatrix_opencl_broadcast_ewise(lhs, v = v, margin = margin, op = op)
    },
    rowSums = function(x, na.rm = FALSE, dims = 1L) {
      if (isTRUE(na.rm) || !identical(dims, 1L)) {
        return(base::rowSums(as.matrix(x), na.rm = na.rm, dims = dims))
      }
      amatrix_opencl_axis_sums(x, axis = 0L)
    },
    colSums = function(x, na.rm = FALSE, dims = 1L) {
      if (isTRUE(na.rm) || !identical(dims, 1L)) {
        return(base::colSums(as.matrix(x), na.rm = na.rm, dims = dims))
      }
      amatrix_opencl_axis_sums(x, axis = 1L)
    },
    solve = function(x, b = NULL, ...) {
      amatrix_opencl_solve(x, b = b)
    },
    chol = function(x, ...) {
      amatrix_opencl_chol(x)
    },
    qr = function(x, ...) {
      amatrix_opencl_qr(x, ...)
    },
    svd = function(x, nu, nv, LINPACK = FALSE, ...) {
      amatrix_opencl_svd(x, nu = nu, nv = nv, LINPACK = LINPACK, ...)
    },
    eigen = function(x, symmetric, only.values = FALSE, EISPACK = FALSE) {
      amatrix_opencl_eigen(x, symmetric = symmetric, only.values = only.values, EISPACK = EISPACK)
    },
    covariance = function(x, center = TRUE, denom = NULL, ...) {
      amatrix_opencl_covariance(x, center = center, denom = denom)
    },
    resident_has = function(key) {
      amatrix_opencl_resident_has(key)
    },
    resident_store = function(key, x) {
      amatrix_opencl_resident_store(key, x)
    },
    resident_drop = function(key) {
      amatrix_opencl_resident_drop(key)
    },
    resident_materialize = function(key) {
      amatrix_opencl_resident_materialize(key)
    },
    sparse_resident_store = function(key, x_sp) {
      amatrix_opencl_sparse_resident_store(key, x_sp)
    },
    sparse_resident_has = function(key) {
      isTRUE(amatrix_opencl_sparse_resident_has(key))
    },
    sparse_resident_drop = function(key) {
      amatrix_opencl_sparse_resident_drop(key)
    },
    spmm_resident = function(sp_key, B, trans_lhs = FALSE) {
      amatrix_opencl_spmm_resident(sp_key, B, trans_lhs = trans_lhs)
    },
    spmm_resident_key = function(sp_key, y_key, out_key, trans_lhs = FALSE, defer = FALSE) {
      amatrix_opencl_spmm_resident_key(sp_key, y_key, out_key, trans_lhs = trans_lhs, defer = defer)
    },
    matmul_resident = function(x_key, y_key, out_key, defer = FALSE) {
      amatrix_opencl_matmul_resident(x_key, y_key, out_key, defer = defer)
    },
    crossprod_resident = function(x_key, y_key = NULL, out_key, defer = FALSE) {
      amatrix_opencl_crossprod_resident(x_key, y_key = y_key, out_key = out_key, defer = defer)
    },
    tcrossprod_resident = function(x_key, y_key = NULL, out_key, defer = FALSE) {
      amatrix_opencl_tcrossprod_resident(x_key, y_key = y_key, out_key = out_key, defer = defer)
    },
    ewise_resident = function(lhs_key, rhs, op, out_key, defer = FALSE) {
      amatrix_opencl_ewise_resident(lhs_key, rhs, op, out_key, defer = defer)
    },
    broadcast_ewise_resident = function(lhs_key, v, margin, op, out_key, defer = FALSE) {
      amatrix_opencl_broadcast_ewise_resident(lhs_key, v, margin, op, out_key, defer = defer)
    },
    broadcast_ewise_resident_inplace = function(lhs_key, v, margin, op) {
      amatrix_opencl_broadcast_ewise_resident_inplace(lhs_key, v, margin, op)
    },
    broadcast_ewise_resident_key = function(lhs_key, v_key, margin, op, out_key, defer = FALSE) {
      amatrix_opencl_broadcast_ewise_resident(lhs_key, v_key, margin, op, out_key, defer = defer)
    },
    broadcast_ewise_resident_inplace_key = function(lhs_key, v_key, margin, op) {
      amatrix_opencl_broadcast_ewise_resident_inplace(lhs_key, v_key, margin, op)
    },
    rowSums_resident = function(x_key, na.rm = FALSE, dims = 1L) {
      amatrix_opencl_rowSums_resident(x_key, na.rm = na.rm, dims = dims)
    },
    rowSums_resident_key = function(x_key, out_key, na.rm = FALSE, dims = 1L) {
      amatrix_opencl_rowSums_resident_key(x_key, out_key, na.rm = na.rm, dims = dims)
    },
    colSums_resident = function(x_key, na.rm = FALSE, dims = 1L) {
      amatrix_opencl_colSums_resident(x_key, na.rm = na.rm, dims = dims)
    },
    colSums_resident_key = function(x_key, out_key, na.rm = FALSE, dims = 1L) {
      amatrix_opencl_colSums_resident_key(x_key, out_key, na.rm = na.rm, dims = dims)
    },
    solve_resident = function(a_key, b_key = NULL, out_key) {
      amatrix_opencl_solve_resident(a_key, b_key = b_key, out_key = out_key)
    },
    solve_triangular_resident = function(factor_key, rhs_key, out_key, lower = FALSE, transpose = FALSE, defer = FALSE) {
      amatrix_opencl_solve_triangular_resident(
        factor_key,
        rhs_key,
        out_key,
        lower = lower,
        transpose = transpose,
        defer = defer
      )
    },
    chol_resident = function(x_key, out_key) {
      amatrix_opencl_chol_resident(x_key, out_key = out_key)
    },
    qr_Q_resident = function(x_key, q_out_key) {
      amatrix_opencl_qr_Q_resident(x_key, q_out_key = q_out_key)
    },
    chol_solve_factor = function(R, B) {
      amatrix_opencl_chol_solve_factor(R, B)
    },
    solve_triangular_factor = function(R, B, lower = FALSE, transpose = FALSE) {
      amatrix_opencl_solve_triangular_factor(R, B, lower = lower, transpose = transpose)
    }
  )
}

amatrix_opencl_register <- function(overwrite = TRUE) {
  amatrix_register_backend("opencl", amatrix_opencl_backend(), overwrite = overwrite)
  invisible("opencl")
}
