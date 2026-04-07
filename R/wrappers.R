.amatrix_rewrap_like <- function(template, value) {
  if (inherits(template, "adgeMatrix")) {
    return(new_adgeMatrix(
      value,
      preferred_backend = template@preferred_backend,
      policy = template@policy,
      precision = template@precision
    ))
  }

  if (inherits(template, "adgCMatrix")) {
    if (inherits(value, "sparseMatrix")) {
      return(new_adgCMatrix(
        as(value, "dgCMatrix"),
        preferred_backend = template@preferred_backend,
        policy = template@policy,
        precision = template@precision
      ))
    }
    return(new_adgeMatrix(
      as.matrix(value),
      preferred_backend = template@preferred_backend,
      policy = template@policy,
      precision = template@precision
    ))
  }

  value
}

.amatrix_host_arg <- function(value) {
  if (inherits(value, "aMatrix")) {
    return(amatrix_materialize_host(value))
  }
  value
}

.amatrix_is_numeric_matrix_value <- function(value) {
  if (inherits(value, "Matrix")) {
    if ("x" %in% slotNames(value)) {
      return(is.numeric(value@x))
    }
    return(FALSE)
  }

  if (is.matrix(value)) {
    return(is.numeric(value))
  }

  FALSE
}

.amatrix_rewrap_value <- function(template, value) {
  if ((inherits(value, "Matrix") || is.matrix(value)) && .amatrix_is_numeric_matrix_value(value)) {
    return(.amatrix_rewrap_like(template, value))
  }
  value
}

.amatrix_template <- function(e1, e2 = NULL) {
  if (inherits(e1, "aMatrix")) {
    return(e1)
  }

  if (inherits(e2, "aMatrix")) {
    return(e2)
  }

  NULL
}

.amatrix_is_dense_matrix_like <- function(value) {
  inherits(value, "adgeMatrix") || inherits(value, "dgeMatrix") || is.matrix(value)
}

.amatrix_prepare_resident_arg <- function(value, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_residency_capable(backend)) {
    return(NULL)
  }

  if (inherits(value, "adgeMatrix")) {
    resident_key <- .amatrix_resident_key(value, backend = backend_name)
    if (!is.null(resident_key) && isTRUE(backend$resident_has(resident_key))) {
      return(list(key = resident_key, temporary = FALSE, tracked = TRUE))
    }

    resident_key <- .amatrix_next_resident_key(backend_name)
    backend$resident_store(resident_key, amatrix_materialize_host(value))
    .amatrix_bind_resident(value, backend_name, resident_key)
    return(list(key = resident_key, temporary = FALSE, tracked = TRUE))
  }

  if (inherits(value, "aMatrix")) {
    return(NULL)
  }

  if (.amatrix_is_dense_matrix_like(value)) {
    resident_key <- .amatrix_next_resident_key(backend_name)
    backend$resident_store(resident_key, .amatrix_host_arg(value))
    return(list(key = resident_key, temporary = TRUE, tracked = FALSE))
  }

  NULL
}

.amatrix_cleanup_temp_resident <- function(args, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_residency_capable(backend)) {
    return(invisible(NULL))
  }

  for (arg in args) {
    if (!is.null(arg) && isTRUE(arg$temporary) && isTRUE(backend$resident_has(arg$key))) {
      backend$resident_drop(arg$key)
    }
  }

  invisible(NULL)
}

.amatrix_try_resident_matmul <- function(x, y, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_supports_resident_op(backend, "matmul")) {
    return(NULL)
  }

  lhs <- .amatrix_prepare_resident_arg(x, backend_name)
  rhs <- .amatrix_prepare_resident_arg(y, backend_name)
  if (is.null(lhs) || is.null(rhs)) {
    .amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name)
    return(NULL)
  }

  out_key <- .amatrix_next_resident_key(backend_name)
  value <- backend$matmul_resident(lhs$key, rhs$key, out_key)
  .amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name)

  list(value = value, backend = backend_name, resident_key = out_key)
}

.amatrix_try_resident_crossprod <- function(x, y, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_supports_resident_op(backend, "crossprod")) {
    return(NULL)
  }

  lhs <- .amatrix_prepare_resident_arg(x, backend_name)
  rhs <- if (is.null(y)) NULL else .amatrix_prepare_resident_arg(y, backend_name)
  if (is.null(lhs) || (!is.null(y) && is.null(rhs))) {
    .amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name)
    return(NULL)
  }

  out_key <- .amatrix_next_resident_key(backend_name)
  rhs_key <- if (is.null(rhs)) NULL else rhs$key
  value <- backend$crossprod_resident(lhs$key, rhs_key, out_key)
  .amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name)

  list(value = value, backend = backend_name, resident_key = out_key)
}

.amatrix_try_resident_tcrossprod <- function(x, y, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_supports_resident_op(backend, "tcrossprod")) {
    return(NULL)
  }

  lhs <- .amatrix_prepare_resident_arg(x, backend_name)
  rhs <- if (is.null(y)) NULL else .amatrix_prepare_resident_arg(y, backend_name)
  if (is.null(lhs) || (!is.null(y) && is.null(rhs))) {
    .amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name)
    return(NULL)
  }

  out_key <- .amatrix_next_resident_key(backend_name)
  rhs_key <- if (is.null(rhs)) NULL else rhs$key
  value <- backend$tcrossprod_resident(lhs$key, rhs_key, out_key)
  .amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name)

  list(value = value, backend = backend_name, resident_key = out_key)
}

.amatrix_try_resident_ewise <- function(op, e1, e2, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_supports_resident_op(backend, "ewise")) {
    return(NULL)
  }

  template <- .amatrix_template(e1, e2)
  if (is.null(template) || !inherits(template, "adgeMatrix")) {
    return(NULL)
  }

  lhs <- if (inherits(e1, "adgeMatrix")) .amatrix_prepare_resident_arg(e1, backend_name) else .amatrix_prepare_resident_arg(e2, backend_name)
  if (is.null(lhs)) {
    return(NULL)
  }

  rhs_arg <- if (inherits(e1, "adgeMatrix")) e2 else e1
  rhs <- NULL
  rhs_payload <- rhs_arg

  if (inherits(rhs_arg, "adgeMatrix") || .amatrix_is_dense_matrix_like(rhs_arg)) {
    rhs <- .amatrix_prepare_resident_arg(rhs_arg, backend_name)
    if (is.null(rhs)) {
      .amatrix_cleanup_temp_resident(list(lhs), backend_name)
      return(NULL)
    }
    rhs_payload <- rhs$key
  } else if (!is.null(rhs_arg) && is.numeric(rhs_arg) && length(rhs_arg) == 1L) {
    rhs_payload <- as.double(rhs_arg)
  } else if (!is.null(rhs_arg)) {
    .amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name)
    return(NULL)
  }

  out_key <- .amatrix_next_resident_key(backend_name)
  value <- backend$ewise_resident(lhs$key, rhs_payload, op, out_key)
  .amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name)

  list(value = value, backend = backend_name, resident_key = out_key)
}

# rowSums/colSums: output is a vector — no resident binding on result.
.amatrix_try_resident_rowSums <- function(x, na.rm, dims, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_supports_resident_op(backend, "rowSums")) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(x, backend_name)
  if (is.null(lhs)) return(NULL)
  result <- backend$rowSums_resident(lhs$key, na.rm, dims)
  .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  result
}

.amatrix_try_resident_colSums <- function(x, na.rm, dims, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_supports_resident_op(backend, "colSums")) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(x, backend_name)
  if (is.null(lhs)) return(NULL)
  result <- backend$colSums_resident(lhs$key, na.rm, dims)
  .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  result
}

# am_solve: output is a matrix; store at out_key and bind resident.
.amatrix_try_resident_solve <- function(a, b, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_supports_resident_op(backend, "solve")) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(a, backend_name)
  if (is.null(lhs)) return(NULL)
  out_key <- .amatrix_next_resident_key(backend_name)
  if (is.null(b)) {
    value <- backend$solve_resident(lhs$key, NULL, out_key)
    .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  } else {
    b_arg <- if (is.vector(b)) matrix(b, ncol = 1L) else b
    rhs <- .amatrix_prepare_resident_arg(b_arg, backend_name)
    if (is.null(rhs)) {
      .amatrix_cleanup_temp_resident(list(lhs), backend_name)
      return(NULL)
    }
    value <- backend$solve_resident(lhs$key, rhs$key, out_key)
    .amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name)
  }
  list(value = value, backend = backend_name, resident_key = out_key)
}

# chol: output is a matrix; store at out_key and bind resident.
.amatrix_try_resident_chol <- function(x, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_supports_resident_op(backend, "chol")) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(x, backend_name)
  if (is.null(lhs)) return(NULL)
  out_key <- .amatrix_next_resident_key(backend_name)
  value <- backend$chol_resident(lhs$key, out_key)
  .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  list(value = value, backend = backend_name, resident_key = out_key)
}

matmul <- function(x, y) {
  # irlba's hot path passes a plain numeric vector for v (A %*% v).
  # _prepare_resident_arg only accepts matrices, so without promotion the
  # resident path silently fails: _try_resident_matmul returns NULL,
  # amatrix_dispatch_op drops the resident binding, and A is re-uploaded on
  # every Lanczos step.  Promoting to a column matrix fixes that.
  y_vec <- is.numeric(y) && is.null(dim(y))
  y_eff <- if (y_vec) matrix(y, ncol = 1L) else y

  choice <- .amatrix_backend_for(x, "matmul", y = y_eff)
  resident <- .amatrix_try_resident_matmul(x, y_eff, choice$name)
  if (!is.null(resident)) {
    if (y_vec) {
      # Result is an m×1 matrix stored at out_key on device. Materialize to host,
      # squeeze to vector, then free the out_key (not useful as a resident matrix).
      # A's key is marked non-temporary so it stays resident for the next call.
      bk <- .amatrix_get_backend(resident$backend)
      mat_result <- bk$resident_materialize(resident$resident_key)
      if (isTRUE(bk$resident_has(resident$resident_key))) {
        bk$resident_drop(resident$resident_key)
      }
      return(drop(mat_result))
    }
    value <- .amatrix_rewrap_like(x, resident$value)
    return(.amatrix_bind_resident(value, resident$backend, resident$resident_key))
  }

  .amatrix_rewrap_like(
    x,
    amatrix_dispatch_op(
      x = x,
      op = "matmul",
      method = "matmul",
      y = y,
      args = list(y = .amatrix_host_arg(y)),
      fallback = function() amatrix_materialize_host(x) %*% .amatrix_host_arg(y)
    )
  )
}

am_crossprod <- function(x, y = NULL, ...) {
  if (!inherits(x, "aMatrix")) {
    if (is.null(y)) return(base::crossprod(x, ...))
    return(base::crossprod(x, y = y, ...))
  }
  choice <- .amatrix_backend_for(x, "crossprod", y = y)
  resident <- .amatrix_try_resident_crossprod(x, y, choice$name)
  if (!is.null(resident)) {
    value <- .amatrix_rewrap_like(x, resident$value)
    return(.amatrix_bind_resident(value, resident$backend, resident$resident_key))
  }

  .amatrix_rewrap_like(
    x,
    amatrix_dispatch_op(
      x = x,
      op = "crossprod",
      method = "crossprod",
      y = y,
      args = list(y = .amatrix_host_arg(y), ...),
      fallback = function() {
        if (is.null(y)) {
          return(base::crossprod(as.matrix(amatrix_materialize_host(x)), ...))
        }
        base::crossprod(as.matrix(amatrix_materialize_host(x)), y = .amatrix_host_arg(y), ...)
      }
    )
  )
}

am_tcrossprod <- function(x, y = NULL, ...) {
  if (!inherits(x, "aMatrix")) {
    if (is.null(y)) return(base::tcrossprod(x, ...))
    return(base::tcrossprod(x, y = y, ...))
  }
  choice <- .amatrix_backend_for(x, "tcrossprod", y = y)
  resident <- .amatrix_try_resident_tcrossprod(x, y, choice$name)
  if (!is.null(resident)) {
    value <- .amatrix_rewrap_like(x, resident$value)
    return(.amatrix_bind_resident(value, resident$backend, resident$resident_key))
  }

  .amatrix_rewrap_like(
    x,
    amatrix_dispatch_op(
      x = x,
      op = "tcrossprod",
      method = "tcrossprod",
      y = y,
      args = list(y = .amatrix_host_arg(y), ...),
      fallback = function() {
        if (is.null(y)) {
          return(base::tcrossprod(as.matrix(amatrix_materialize_host(x)), ...))
        }
        base::tcrossprod(as.matrix(amatrix_materialize_host(x)), y = .amatrix_host_arg(y), ...)
      }
    )
  )
}

# gemm: full BLAS DGEMM-style control surface.
# Computes alpha * op(A) %*% op(B) + beta * C.
# op(X) = t(X) when transX = TRUE, otherwise X.
# Routes each case to the most efficient resident operation:
#   transA only  → crossprod_resident  (t(A) %*% B)
#   transB only  → tcrossprod_resident (A %*% t(B))
#   transA+B     → t(B %*% A)  (tcrossprod identity; no host copy of t(B))
#   neither      → matmul_resident     (A %*% B)
# Use matmul / am_crossprod / am_tcrossprod for plain operator idioms.
gemm <- function(A, B, C = NULL, alpha = 1.0, beta = 1.0,
                    transA = FALSE, transB = FALSE) {
  AB <- if (transA && transB) {
    # t(A) %*% t(B) = t(B %*% A): use identity to avoid materialising t(B) to host.
    am_transpose(matmul(B, A))
  } else if (transA) {
    am_crossprod(A, B)
  } else if (transB) {
    am_tcrossprod(A, B)
  } else {
    matmul(A, B)
  }

  if (alpha != 1.0) AB <- ewise("*", AB, alpha)

  if (!is.null(C)) {
    C_scaled <- if (beta != 1.0) ewise("*", C, beta) else C
    ewise("+", AB, C_scaled)
  } else {
    AB
  }
}

rowsums <- function(x, na.rm = FALSE, dims = 1L) {
  choice <- .amatrix_backend_for(x, "rowSums")
  resident <- .amatrix_try_resident_rowSums(x, na.rm, dims, choice$name)
  if (!is.null(resident)) return(resident)
  amatrix_dispatch_op(
    x = x,
    op = "rowSums",
    method = "rowSums",
    args = list(na.rm = na.rm, dims = dims),
    fallback = function() Matrix::rowSums(amatrix_materialize_host(x), na.rm = na.rm, dims = dims)
  )
}

colsums <- function(x, na.rm = FALSE, dims = 1L) {
  choice <- .amatrix_backend_for(x, "colSums")
  resident <- .amatrix_try_resident_colSums(x, na.rm, dims, choice$name)
  if (!is.null(resident)) return(resident)
  amatrix_dispatch_op(
    x = x,
    op = "colSums",
    method = "colSums",
    args = list(na.rm = na.rm, dims = dims),
    fallback = function() Matrix::colSums(amatrix_materialize_host(x), na.rm = na.rm, dims = dims)
  )
}

am_transpose <- function(x) {
  if (inherits(x, "adgeMatrix")) {
    return(.new_aTransposeView(x))
  }
  # Sparse and other types: materialize and transpose on host
  .amatrix_rewrap_like(x, t(as.matrix(amatrix_materialize_host(x))))
}

am_subset <- function(x, i, j, ..., drop = TRUE) {
  value <- amatrix_materialize_host(x)[i, j, ..., drop = drop]
  .amatrix_rewrap_value(x, value)
}

am_subassign <- function(x, i, j, ..., value) {
  host_x <- amatrix_materialize_host(x)
  host_value <- .amatrix_host_arg(value)
  host_x[i, j, ...] <- host_value
  .amatrix_rewrap_value(x, host_x)
}

am_solve <- function(a, b = NULL, ...) {
  b_arg <- if (missing(b)) NULL else b
  if (!inherits(a, "aMatrix")) {
    if (is.null(b_arg)) return(base::solve(a, ...))
    return(base::solve(a, b = b_arg, ...))
  }
  b_was_vector <- !is.null(b_arg) && is.numeric(b_arg) && is.null(dim(b_arg))
  choice <- .amatrix_backend_for(a, "solve", y = b_arg)
  resident <- .amatrix_try_resident_solve(a, b_arg, choice$name)
  if (!is.null(resident)) {
    if (b_was_vector) return(as.vector(resident$value))
    value <- .amatrix_rewrap_value(a, resident$value)
    return(.amatrix_bind_resident(value, resident$backend, resident$resident_key))
  }

  if (is.null(b_arg)) {
    return(.amatrix_rewrap_value(
      a,
      amatrix_dispatch_op(
        x = a,
        op = "solve",
        method = "solve",
        args = list(...),
        fallback = function() base::solve(as.matrix(amatrix_materialize_host(a)), ...)
      )
    ))
  }

  result <- amatrix_dispatch_op(
    x = a,
    op = "solve",
    method = "solve",
    y = b_arg,
    args = list(b = .amatrix_host_arg(b_arg), ...),
    fallback = function() base::solve(as.matrix(amatrix_materialize_host(a)), as.matrix(.amatrix_host_arg(b_arg)), ...)
  )
  # When b was a plain vector, base::solve() returns a named numeric vector.
  # Preserve that contract: don't wrap back to adgeMatrix.
  if (b_was_vector) as.vector(result) else .amatrix_rewrap_value(a, result)
}

am_chol <- function(x, ...) {
  choice <- .amatrix_backend_for(x, "chol")
  resident <- .amatrix_try_resident_chol(x, choice$name)
  if (!is.null(resident)) {
    value <- .amatrix_rewrap_value(x, resident$value)
    return(.amatrix_bind_resident(value, resident$backend, resident$resident_key))
  }
  .amatrix_rewrap_value(
    x,
    amatrix_dispatch_op(
      x = x,
      op = "chol",
      method = "chol",
      args = list(...),
      fallback = function() base::chol(as.matrix(amatrix_materialize_host(x)), ...)
    )
  )
}

am_qr <- function(x, ...) {
  qr_value <- amatrix_dispatch_op(
    x = x,
    op = "qr",
    method = "qr",
    args = list(...),
    fallback = function() base::qr(as.matrix(amatrix_materialize_host(x)), ...)
  )

  .amatrix_wrap_qr(qr_value, x)
}

# ── Resident QR helper ────────────────────────────────────────────────────────
#
# Given an adgeMatrix z_am that is already (or can be) resident on the backend,
# run QR and return Q as a new resident adgeMatrix — no R-memory round-trip for
# the QR step itself.
#
# Falls back to qr.Q(qr(as.matrix(z_am))) + re-upload when the backend does
# not support qr_Q_resident (e.g. the cpu backend).
.amatrix_try_resident_qr_Q <- function(z_am) {
  if (!inherits(z_am, "adgeMatrix")) {
    z_am <- adgeMatrix(as.matrix(z_am))
  }
  choice       <- .amatrix_backend_for(z_am, "qr")
  backend_name <- choice$name
  backend      <- choice$backend

  if (is.function(backend$qr_Q_resident) &&
      .amatrix_backend_residency_capable(backend)) {
    z_info <- .amatrix_prepare_resident_arg(z_am, backend_name)
    if (!is.null(z_info)) {
      q_key <- .amatrix_next_resident_key(backend_name)
      backend$qr_Q_resident(z_info$key, q_key)
      .amatrix_cleanup_temp_resident(list(z_info), backend_name)
      # Materialize Q for the host slot (dims/fallback paths).
      # Q is p×k_over — small relative to the data matrix.
      q_mat <- backend$resident_materialize(q_key)
      q_am  <- new_adgeMatrix(q_mat,
                  preferred_backend = backend_name,
                  policy            = z_am@policy,
                  precision         = z_am@precision)
      .amatrix_bind_resident(q_am, backend_name, q_key)
      return(q_am)
    }
  }
  {
    q_mat <- qr.Q(qr(as.matrix(amatrix_materialize_host(z_am))))
    adgeMatrix(q_mat,
      preferred_backend = backend_name,
      policy            = z_am@policy,
      precision         = z_am@precision)
  }
}

am_svd <- function(x, nu = min(dim(x)), nv = min(dim(x)), LINPACK = FALSE, ...) {
  amatrix_dispatch_op(
    x = x,
    op = "svd",
    method = "svd",
    args = list(nu = nu, nv = nv, LINPACK = LINPACK, ...),
    fallback = function() {
      if (isTRUE(LINPACK)) {
        stop("LINPACK is not supported", call. = FALSE)
      }
      base::svd(as.matrix(amatrix_materialize_host(x)), nu = nu, nv = nv, ...)
    }
  )
}

am_eigen <- function(x, symmetric = NULL, only.values = FALSE, EISPACK = FALSE) {
  # Mirror base::eigen behaviour: if symmetric is not supplied, auto-detect
  # from the host matrix so callers don't have to know the structure.
  if (is.null(symmetric)) {
    x_host <- as.matrix(amatrix_materialize_host(x))
    symmetric <- isSymmetric(x_host)
  }
  amatrix_dispatch_op(
    x = x,
    op = "eigen",
    method = "eigen",
    args = list(symmetric = symmetric, only.values = only.values, EISPACK = EISPACK),
    fallback = function() base::eigen(as.matrix(amatrix_materialize_host(x)), symmetric = symmetric, only.values = only.values, EISPACK = EISPACK)
  )
}

eigh <- function(x) {
  am_eigen(x, symmetric = TRUE)
}

# ── Weighted am_crossprod helpers ─────────────────────────────────────────────

# X' diag(w) X  (p x p)
crossprod_weighted <- function(X, w) {
  X_arg <- .amatrix_model_dense_arg(X)
  w <- as.double(w)
  if (length(w) != nrow(X_arg)) {
    stop("length(w) must equal nrow(X)", call. = FALSE)
  }
  sqrt_w <- sqrt(w)
  x_host <- as.matrix(amatrix_materialize_host(X_arg))
  x_scaled <- x_host * sqrt_w
  am_crossprod(.amatrix_rewrap_like(X_arg, x_scaled))
}

# X diag(w) X'  (n x n)
tcrossprod_weighted <- function(X, w) {
  X_arg <- .amatrix_model_dense_arg(X)
  w <- as.double(w)
  if (length(w) != nrow(X_arg)) {
    stop("length(w) must equal nrow(X)", call. = FALSE)
  }
  sqrt_w <- sqrt(w)
  x_host <- as.matrix(amatrix_materialize_host(X_arg))
  x_scaled <- x_host * sqrt_w
  am_tcrossprod(.amatrix_rewrap_like(X_arg, x_scaled))
}

# X' diag(w) y  (p x k)
xty_weighted <- function(X, w, y) {
  X_arg <- .amatrix_model_dense_arg(X)
  w <- as.double(w)
  if (length(w) != nrow(X_arg)) {
    stop("length(w) must equal nrow(X)", call. = FALSE)
  }
  sqrt_w <- sqrt(w)
  x_host <- as.matrix(amatrix_materialize_host(X_arg))
  y_mat <- if (is.vector(y)) matrix(y, ncol = 1L) else as.matrix(y)
  if (nrow(y_mat) != nrow(X_arg)) {
    stop("nrow(y) must equal nrow(X)", call. = FALSE)
  }
  x_scaled <- x_host * sqrt_w
  y_scaled <- y_mat * sqrt_w
  am_crossprod(
    .amatrix_rewrap_like(X_arg, x_scaled),
    .amatrix_rewrap_like(X_arg, y_scaled)
  )
}

# ── Diagonal scaling ──────────────────────────────────────────────────────────

# diag(d) %*% X  — scale row i by d[i]
rowscale <- function(X, d) {
  X_arg <- .amatrix_model_dense_arg(X)
  d <- as.double(d)
  if (length(d) != nrow(X_arg)) {
    stop("length(d) must equal nrow(X)", call. = FALSE)
  }
  x_host <- as.matrix(amatrix_materialize_host(X_arg))
  .amatrix_rewrap_like(X_arg, x_host * d)
}

# X %*% diag(d)  — scale col j by d[j]
colscale <- function(X, d) {
  X_arg <- .amatrix_model_dense_arg(X)
  d <- as.double(d)
  if (length(d) != ncol(X_arg)) {
    stop("length(d) must equal ncol(X)", call. = FALSE)
  }
  x_host <- as.matrix(amatrix_materialize_host(X_arg))
  .amatrix_rewrap_like(X_arg, t(t(x_host) * d))
}

am_diag <- function(x, nrow, ncol, names = TRUE) {
  # Extract mode: diag(matrix_x) → numeric vector of diagonal elements
  extract_mode <- (missing(nrow) && missing(ncol)) &&
                  (is.matrix(x) || inherits(x, "aMatrix"))
  if (extract_mode) {
    x_host <- as.matrix(amatrix_materialize_host(x))
    return(base::diag(x_host, names = names))
  }
  # Create mode: diag(d) → diagonal adgeMatrix; nrow/ncol set the size
  nrow <- if (missing(nrow)) NULL else nrow
  ncol <- if (missing(ncol)) NULL else ncol
  amatrix_dispatch_op(
    x = x,
    op = "diag",
    method = "diag",
    args = list(nrow = nrow, ncol = ncol, names = names),
    fallback = function() {
      args <- list(as.matrix(amatrix_materialize_host(x)), names = names)
      if (!is.null(nrow)) args$nrow <- nrow
      if (!is.null(ncol)) args$ncol <- ncol
      do.call(base::diag, args)
    }
  )
}

# ── Fused crossprod + diagonal add ───────────────────────────────────────────

# X'X + lambda*I  or  X'X + diag(d)
crossprod_add_diag <- function(X, lambda) {
  X_arg  <- .amatrix_model_dense_arg(X)
  xtx    <- am_crossprod(X_arg)
  p      <- ncol(X_arg)
  xtx_m  <- as.matrix(amatrix_materialize_host(xtx))
  if (length(lambda) == 1L) {
    diag(xtx_m) <- diag(xtx_m) + as.double(lambda)
  } else {
    if (length(lambda) != p)
      stop("lambda must be a scalar or length ncol(X)", call. = FALSE)
    diag(xtx_m) <- diag(xtx_m) + as.double(lambda)
  }
  .amatrix_rewrap_like(X_arg, xtx_m)
}

# ── Matrix functions (via symmetric eigendecomposition) ───────────────────────

.mat_fun <- function(X, f, check_positive = TRUE) {
  X_arg <- .amatrix_model_dense_arg(X)
  res   <- eigh(X_arg)
  lam   <- res$values
  if (isTRUE(check_positive) && any(lam <= 0))
    warning("matrix has non-positive eigenvalues; result may be complex or NaN",
            call. = FALSE)
  Q     <- res$vectors
  new_lam <- f(lam)
  .amatrix_rewrap_like(X_arg, Q %*% diag(new_lam) %*% t(Q))
}

mat_sqrt <- function(X) .mat_fun(X, sqrt)
mat_pow  <- function(X, p) .mat_fun(X, function(lam) lam^p)
mat_log  <- function(X) .mat_fun(X, log)

# ── Stochastic trace estimator (Hutchinson) ───────────────────────────────────

# Estimates tr(A) or tr(A^{-1}) via k Rademacher probes.
# For tr(A):         trace_estim(A, k)
# For tr(K^{-1}):   trace_estim(solve_fn = function(v) chol_solve(L, v), n, k)
trace_estim <- function(A = NULL, k = 30L, seed = NULL,
                        solve_fn = NULL, n = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (!is.null(solve_fn)) {
    if (is.null(n)) stop("n must be supplied when using solve_fn", call. = FALSE)
    probes <- matrix(sample(c(-1, 1), n * k, replace = TRUE), n, k)
    sols   <- solve_fn(probes)
    return(mean(colSums(probes * as.matrix(sols))))
  }

  A_host <- as.matrix(amatrix_materialize_host(A))
  n      <- nrow(A_host)
  probes <- matrix(sample(c(-1, 1), n * k, replace = TRUE), n, k)
  mean(colSums(probes * (A_host %*% probes)))
}

# ── Row / column means ────────────────────────────────────────────────────────

rowmeans <- function(x, na.rm = FALSE) {
  x_host <- as.matrix(amatrix_materialize_host(x))
  base::rowMeans(x_host, na.rm = na.rm)
}

colmeans <- function(x, na.rm = FALSE) {
  x_host <- as.matrix(amatrix_materialize_host(x))
  base::colMeans(x_host, na.rm = na.rm)
}

# ── Matrix trace ──────────────────────────────────────────────────────────────

trace <- function(x) {
  x_host <- as.matrix(amatrix_materialize_host(x))
  sum(base::diag(x_host))
}

# ── Symmetry enforcement ──────────────────────────────────────────────────────

sym <- function(x) {
  X_arg <- .amatrix_model_dense_arg(x)
  x_host <- as.matrix(amatrix_materialize_host(X_arg))
  .amatrix_rewrap_like(X_arg, (x_host + t(x_host)) / 2)
}

# ── Inner product ─────────────────────────────────────────────────────────────

dot <- function(x, y) {
  x_host <- .amatrix_host_arg(x)
  y_host <- .amatrix_host_arg(y)
  sum(x_host * y_host)
}

# ── segment_sum / segment_mean (amatrix-ylo) ──────────────────────────────
# First-class grouped-reduction primitives.  GPU path stores result as a
# resident adgeMatrix (no data downloaded). CPU path uses base::rowsum.

.am_segment_sum_cpu <- function(X_mat, labels, K) {
  sums_raw <- rowsum(X_mat, labels, reorder = FALSE)
  out <- matrix(0, K, ncol(X_mat))
  idx <- as.integer(rownames(sums_raw))
  valid <- idx >= 1L & idx <= K
  if (any(valid)) out[idx[valid], ] <- sums_raw[valid, , drop = FALSE]
  out
}

.am_segment_mean_cpu <- function(X_mat, labels, K) {
  sums_raw <- rowsum(X_mat, labels, reorder = FALSE)
  idx    <- as.integer(rownames(sums_raw))
  counts <- tabulate(labels, nbins = K)
  out    <- matrix(NA_real_, K, ncol(X_mat))
  valid  <- idx >= 1L & idx <= K
  if (any(valid)) {
    k  <- idx[valid]
    nz <- counts[k] > 0L
    if (any(nz))
      out[k[nz], ] <- sums_raw[valid, , drop = FALSE][nz, , drop = FALSE] /
                      counts[k[nz]]
  }
  out
}

.am_try_resident_segment_op <- function(x, labels, K, backend_name, op_name) {
  backend <- .amatrix_get_backend(backend_name)
  fn_name <- paste0(op_name, "_resident")
  if (!is.function(backend[[fn_name]])) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(x, backend_name)
  if (is.null(lhs)) return(NULL)
  out_key <- .amatrix_next_resident_key(backend_name)
  value   <- backend[[fn_name]](lhs$key, labels, K, out_key)
  .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  # Bridge returns a plain R matrix; fall back to resident_materialize if not
  if (!is.matrix(value)) value <- backend$resident_materialize(out_key)
  list(value = value, key = out_key)
}

.am_segment_resident_wrap <- function(x, resident, choice_name) {
  wrapped <- new_adgeMatrix(resident$value,
                            preferred_backend = x@preferred_backend,
                            precision        = x@precision,
                            policy           = x@policy)
  .amatrix_bind_resident(wrapped, choice_name, resident$key)
}

segment_sum <- function(x, labels, K) {
  labels <- as.integer(labels)
  K      <- as.integer(K)
  if (!inherits(x, "adgeMatrix"))
    return(.am_segment_sum_cpu(as.matrix(x), labels, K))
  choice   <- .amatrix_backend_for(x, "segment_sum")
  resident <- .am_try_resident_segment_op(x, labels, K, choice$name, "segment_sum")
  if (is.null(resident))
    return(.am_segment_sum_cpu(as.matrix(amatrix_materialize_host(x)), labels, K))
  .am_segment_resident_wrap(x, resident, choice$name)
}

segment_mean <- function(x, labels, K) {
  labels <- as.integer(labels)
  K      <- as.integer(K)
  if (!inherits(x, "adgeMatrix"))
    return(.am_segment_mean_cpu(as.matrix(x), labels, K))
  choice   <- .amatrix_backend_for(x, "segment_mean")
  resident <- .am_try_resident_segment_op(x, labels, K, choice$name, "segment_mean")
  if (is.null(resident))
    return(.am_segment_mean_cpu(as.matrix(amatrix_materialize_host(x)), labels, K))
  .am_segment_resident_wrap(x, resident, choice$name)
}

# ── addmm (amatrix-uaj) ─────────────────────────────────────────────────
# alpha*(A%*%B) + beta*C  — BLAS-3 fused scaled matmul with optional bias.
# A: n×p adgeMatrix (resident if GPU); B: p×k R matrix; C: n×k R matrix or NULL.
# GPU path uses mlx_addmm directly; CPU path uses plain R arithmetic.

.am_addmm_cpu <- function(A_mat, B_mat, C_mat, alpha, beta) {
  result <- alpha * (A_mat %*% B_mat)
  if (!is.null(C_mat) && beta != 0) result <- result + beta * C_mat
  result
}

.am_try_addmm_gpu <- function(A, B_mat, C_mat, alpha, beta, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!is.function(backend$addmm_resident)) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(A, backend_name)
  if (is.null(lhs)) return(NULL)
  out_key <- .amatrix_next_resident_key(backend_name)
  value   <- backend$addmm_resident(lhs$key, B_mat, C_mat, alpha, beta, out_key)
  .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  if (!is.matrix(value)) value <- backend$resident_materialize(out_key)
  list(value = value, key = out_key)
}

#' Scaled matrix multiply with optional bias: alpha*(A\%*\%B) + beta*C
#'
#' @param A  n×p \code{adgeMatrix} or plain matrix.
#' @param B  p×k numeric matrix.
#' @param C  n×k numeric matrix or \code{NULL} (treated as zeros).
#' @param alpha Scalar multiplier for \code{A\%*\%B} (default 1).
#' @param beta  Scalar multiplier for \code{C} (default 1).
#' @return \code{adgeMatrix} if A is resident, otherwise plain matrix.
#' @export
addmm <- function(A, B, C = NULL, alpha = 1.0, beta = 1.0) {
  B_mat <- as.matrix(B); storage.mode(B_mat) <- "double"
  C_mat <- if (!is.null(C)) { m <- as.matrix(C); storage.mode(m) <- "double"; m } else NULL

  if (!inherits(A, "adgeMatrix")) {
    A_mat <- as.matrix(A); storage.mode(A_mat) <- "double"
    return(.am_addmm_cpu(A_mat, B_mat, C_mat, alpha, beta))
  }

  choice   <- .amatrix_backend_for(A, "addmm")
  resident <- .am_try_addmm_gpu(A, B_mat, C_mat, alpha, beta, choice$name)
  if (!is.null(resident)) {
    wrapped <- new_adgeMatrix(resident$value,
                              preferred_backend = A@preferred_backend,
                              precision = A@precision,
                              policy = A@policy)
    return(.amatrix_bind_resident(wrapped, choice$name, resident$key))
  }

  A_mat <- as.matrix(amatrix_materialize_host(A)); storage.mode(A_mat) <- "double"
  .am_addmm_cpu(A_mat, B_mat, C_mat, alpha, beta)
}

# ── pairwise_sqdist_argmin (amatrix-zas) ───────────────────────────────────
# Fused nearest-centroid assignment via the squared-distance identity:
#   D[i,k] = ||xi||^2 - 2*(X@Ct)[i,k] + ||ck||^2
# GPU path chains resident operations (no intermediate host round-trips):
#   1. cross  = matmul_resident(X, Ct)         [n×K]
#   2. neg2   = ewise_resident(cross, -2, "*") [n×K]
#   3. d1     = broadcast_ewise(neg2, x_norms, margin=1, "+") [add row norms]
#   4. d      = broadcast_ewise(d1,   c_norms, margin=2, "+") [add col norms]
#   5. labels = rowargmin_resident(d) + 1L     [0→1-indexed]
# CPU fallback: base R distance matrix + max.col.
.pairwise_sqdist_argmin_cpu <- function(X_mat, Ct_mat, x_norms, c_norms) {
  cross <- X_mat %*% Ct_mat                          # n×K
  D     <- -2 * cross + x_norms + rep(c_norms, each = nrow(X_mat))
  max.col(-D, ties.method = "first")
}

.pairwise_sqdist_argmin_gpu <- function(X, Ct_mat, x_norms, c_norms,
                                            backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  # All required resident ops must be present
  needed <- c("matmul_resident", "ewise_resident",
               "broadcast_ewise_resident", "rowargmin_resident")
  if (!all(vapply(needed, function(f) is.function(backend[[f]]), logical(1L))))
    return(NULL)

  lhs_X  <- .amatrix_prepare_resident_arg(X, backend_name)
  if (is.null(lhs_X)) return(NULL)

  # Upload Ct (p×K plain matrix) as temporary resident
  Ct_key  <- .amatrix_next_resident_key(backend_name)
  backend$resident_store(Ct_key, Ct_mat)
  temps   <- list(lhs_X, list(key = Ct_key, is_temp = TRUE))

  cross_key  <- .amatrix_next_resident_key(backend_name)
  neg2_key   <- .amatrix_next_resident_key(backend_name)
  d1_key     <- .amatrix_next_resident_key(backend_name)
  d_key      <- .amatrix_next_resident_key(backend_name)

  tryCatch({
    backend$matmul_resident(lhs_X$key, Ct_key, cross_key)
    backend$ewise_resident(cross_key, -2.0, "*", neg2_key)
    backend$broadcast_ewise_resident(neg2_key, as.double(x_norms), 1L, "+", d1_key)
    backend$broadcast_ewise_resident(d1_key, as.double(c_norms), 2L, "+", d_key)
    labels0 <- backend$rowargmin_resident(d_key)
    labels0 + 1L   # 0-indexed → 1-indexed
  }, error = function(e) NULL,
  finally = {
    .amatrix_cleanup_temp_resident(temps, backend_name)
    for (k in c(cross_key, neg2_key, d1_key, d_key))
      tryCatch(backend$resident_drop(k), error = function(e) invisible(NULL))
  })
}

#' Nearest-centroid assignment via fused squared-distance computation
#'
#' Computes \eqn{D[i,k] = \|x_i\|^2 - 2 x_i^\top c_k + \|c_k\|^2} and
#' returns \eqn{\arg\min_k D[i,k]} for each row \eqn{i}, 1-indexed.
#' GPU path avoids host round-trips by chaining resident operations.
#'
#' @param X       n×p \code{adgeMatrix} or plain matrix (query points).
#' @param Ct      p×K numeric matrix (centroids, transposed — columns are centroids).
#' @param x_norms Optional n-vector of precomputed \eqn{\|x_i\|^2}. Computed if \code{NULL}.
#' @param c_norms Optional K-vector of precomputed \eqn{\|c_k\|^2}. Computed if \code{NULL}.
#' @return Integer vector of length n, 1-indexed nearest centroid per row.
#' @export
pairwise_sqdist_argmin <- function(X, Ct, x_norms = NULL, c_norms = NULL) {
  Ct_mat  <- as.matrix(Ct);  storage.mode(Ct_mat) <- "double"
  if (is.null(x_norms)) {
    X_mat   <- if (inherits(X, "adgeMatrix")) as.matrix(amatrix_materialize_host(X))
               else as.matrix(X)
    x_norms <- rowSums(X_mat^2)
  } else {
    x_norms <- as.double(x_norms)
    X_mat   <- NULL
  }
  if (is.null(c_norms)) c_norms <- colSums(Ct_mat^2)
  else c_norms <- as.double(c_norms)

  if (!inherits(X, "adgeMatrix")) {
    X_mat <- if (is.null(X_mat)) as.matrix(X) else X_mat
    storage.mode(X_mat) <- "double"
    return(.pairwise_sqdist_argmin_cpu(X_mat, Ct_mat, x_norms, c_norms))
  }

  choice <- .amatrix_backend_for(X, "matmul")
  result <- .pairwise_sqdist_argmin_gpu(X, Ct_mat, x_norms, c_norms, choice$name)
  if (!is.null(result)) return(result)

  # Fallback: materialize and use CPU
  X_mat <- if (is.null(X_mat)) as.matrix(amatrix_materialize_host(X)) else X_mat
  storage.mode(X_mat) <- "double"
  .pairwise_sqdist_argmin_cpu(X_mat, Ct_mat, x_norms, c_norms)
}

.am_scatter_mean_cpu <- function(X_mat, labels, K) {
  p <- ncol(X_mat)
  centroids <- matrix(NA_real_, K, p)
  for (k in seq_len(K)) {
    idx <- which(labels == k)
    if (length(idx) > 0L)
      centroids[k, ] <- colMeans(X_mat[idx, , drop = FALSE])
  }
  centroids
}

.amatrix_try_resident_scatter_mean <- function(x, labels, K, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!is.function(backend$scatter_mean_resident)) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(x, backend_name)
  if (is.null(lhs)) return(NULL)
  result <- backend$scatter_mean_resident(lhs$key, labels, K)
  .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  result
}

am_scatter_mean <- function(x, labels, K) {
  labels <- as.integer(labels)
  K      <- as.integer(K)
  counts <- tabulate(labels, nbins = K)   # O(n), always on CPU

  if (!inherits(x, "adgeMatrix")) {
    return(.am_scatter_mean_cpu(as.matrix(x), labels, K))
  }

  choice <- .amatrix_backend_for(x, "scatter_mean")
  sums   <- .amatrix_try_resident_scatter_mean(x, labels, K, choice$name)

  if (is.null(sums)) {
    return(.am_scatter_mean_cpu(as.matrix(amatrix_materialize_host(x)), labels, K))
  }

  # sums is K×p; divide each row by cluster count
  means <- sums
  nonzero <- counts > 0L
  if (any(nonzero))
    means[nonzero, ] <- sums[nonzero, , drop = FALSE] / counts[nonzero]
  means[!nonzero, ] <- NA_real_
  means
}

.amatrix_try_resident_broadcast_ewise <- function(x, v, margin, op, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_supports_resident_op(backend, "broadcast_ewise")) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(x, backend_name)
  if (is.null(lhs)) return(NULL)
  out_key <- .amatrix_next_resident_key(backend_name)
  value <- backend$broadcast_ewise_resident(lhs$key, v, margin, op, out_key)
  .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  list(value = value, backend = backend_name, resident_key = out_key)
}

am_sweep <- function(x, MARGIN, STATS, FUN = "+") {
  if (!inherits(x, "adgeMatrix")) {
    return(base::sweep(as.matrix(x), MARGIN = MARGIN, STATS = STATS, FUN = FUN))
  }
  op <- if (is.character(FUN) && length(FUN) == 1L) FUN else NULL
  if (is.null(op) || !op %in% c("+", "-", "*", "/")) {
    return(base::sweep(as.matrix(amatrix_materialize_host(x)), MARGIN = MARGIN,
                       STATS = STATS, FUN = FUN))
  }
  if (!is.numeric(STATS) || is.matrix(STATS)) {
    return(base::sweep(as.matrix(amatrix_materialize_host(x)), MARGIN = MARGIN,
                       STATS = STATS, FUN = FUN))
  }
  v <- as.double(STATS)
  choice <- .amatrix_backend_for(x, "broadcast_ewise")
  resident <- .amatrix_try_resident_broadcast_ewise(x, v, MARGIN, op, choice$name)
  if (!is.null(resident)) {
    value <- .amatrix_rewrap_value(x, resident$value)
    return(.amatrix_bind_resident(value, resident$backend, resident$resident_key))
  }
  result <- amatrix_dispatch_op(
    x = x,
    op = "broadcast_ewise",
    method = "broadcast_ewise",
    args = list(lhs = as.matrix(amatrix_materialize_host(x)), v = v,
                margin = MARGIN, op = op),
    fallback = function() base::sweep(as.matrix(amatrix_materialize_host(x)),
                                      MARGIN = MARGIN, STATS = STATS, FUN = FUN)
  )
  .amatrix_rewrap_value(x, result)
}

.amatrix_try_resident_argreduce <- function(x, kind, backend_name) {
  backend  <- .amatrix_get_backend(backend_name)
  fn_name  <- paste0(kind, "_resident")
  if (!is.function(backend[[fn_name]])) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(x, backend_name)
  if (is.null(lhs)) return(NULL)
  result <- backend[[fn_name]](lhs$key)
  .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  result
}

.am_argreduce_cpu <- function(x, kind) {
  mat <- as.matrix(x)
  switch(kind,
    rowargmax = max.col(mat,  ties.method = "first"),
    rowargmin = max.col(-mat, ties.method = "first"),
    colargmax = max.col(t(mat),  ties.method = "first"),
    colargmin = max.col(-t(mat), ties.method = "first")
  )
}

.am_argreduce <- function(x, kind) {
  if (!inherits(x, "adgeMatrix")) return(.am_argreduce_cpu(x, kind))
  choice <- .amatrix_backend_for(x, "argmax")
  result <- .amatrix_try_resident_argreduce(x, kind, choice$name)
  if (!is.null(result)) return(result)
  .am_argreduce_cpu(x, kind)
}

am_rowargmax <- function(x) .am_argreduce(x, "rowargmax")
am_rowargmin <- function(x) .am_argreduce(x, "rowargmin")
am_colargmax <- function(x) .am_argreduce(x, "colargmax")
am_colargmin <- function(x) .am_argreduce(x, "colargmin")

ewise <- function(op, e1, e2 = NULL) {
  template <- .amatrix_template(e1, e2)
  host_e1 <- .amatrix_host_arg(e1)
  host_e2 <- .amatrix_host_arg(e2)

  if (!is.null(template)) {
    choice <- .amatrix_backend_for(template, "ewise", y = e2)
    resident <- .amatrix_try_resident_ewise(op, e1, e2, choice$name)
    if (!is.null(resident)) {
      value <- .amatrix_rewrap_value(template, resident$value)
      return(.amatrix_bind_resident(value, resident$backend, resident$resident_key))
    }
  }

  value <- amatrix_dispatch_op(
    x = template,
    op = "ewise",
    method = "ewise",
    y = e2,
    args = list(lhs = host_e1, rhs = host_e2, op = op),
    fallback = function() {
      if (is.null(e2)) {
        return(do.call(op, list(host_e1)))
      }
      do.call(op, list(host_e1, host_e2))
    }
  )

  if (is.null(template)) {
    return(value)
  }

  .amatrix_rewrap_value(template, value)
}

am_set_dimnames <- function(x, value) {
  host_x <- amatrix_materialize_host(x)
  dimnames(host_x) <- value
  .amatrix_rewrap_like(x, host_x)
}

# ── Distance / Kernel helpers ──────────────────────────────────────────────

.am_as_double_matrix <- function(x) {
  if (inherits(x, c("adgeMatrix", "adgCMatrix")))
    x <- amatrix_materialize_host(x)
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.double(x)) storage.mode(x) <- "double"
  x
}

# GPU dispatch helpers for distance/kernel computation.
# .dist_matrix_sq_gpu: returns squared Euclidean distance matrix [m×n].
# .am_kernel_gpu:  returns kernel matrix [m×n].
# Both use dedicated column-major AF bridges (all dims) or MLX.

.am_af_ok <- function() {
  tryCatch(
    requireNamespace("amatrix.arrayfire", quietly = TRUE) &&
      amatrix.arrayfire::amatrix_arrayfire_is_available(),
    error = function(e) FALSE
  )
}
.am_mlx_ok <- function() {
  tryCatch(
    requireNamespace("amatrix.mlx", quietly = TRUE) &&
      amatrix.mlx::amatrix_mlx_is_available(),
    error = function(e) FALSE
  )
}

# MLX is preferred on Apple Silicon (39x vs 9x speedup in benchmarks).
# AF is the fallback for CUDA/other platforms where MLX is unavailable.

# AF bridge does D² entirely in C (GEMM + rowSums + broadcast + clamp) with no
# R-level allocation overhead, so it's fast regardless of GPU vs CPU backend.
# MLX path only does the GEMM on GPU; subsequent R ops on 3M-element matrices
# add ~150ms overhead for large matrices — only worthwhile when AF unavailable.
.dist_matrix_sq_gpu <- function(X, Y = NULL) {
  if (isTRUE(.am_af_ok()))
    return(.Call("am_af_dist_sq_bridge", X, Y, PACKAGE = "amatrix.arrayfire"))
  Y_eff <- if (is.null(Y)) X else Y
  if (isTRUE(.am_mlx_ok())) {
    G  <- .Call("amatrix_mlx_tcrossprod_bridge", X, Y, PACKAGE = "amatrix.mlx")
    nx <- rowSums(X^2); ny <- if (is.null(Y)) nx else rowSums(Y_eff^2)
    return(pmax(outer(nx, ny, "+") - 2 * G, 0))
  }
  G  <- am_tcrossprod(X, Y)
  nx <- rowSums(X^2); ny <- if (is.null(Y)) nx else rowSums(Y_eff^2)
  pmax(outer(nx, ny, "+") - 2 * G, 0)
}

.am_kernel_gpu <- function(X, Y = NULL, kernel, sigma, degree, coef) {
  # AF bridge computes entire kernel in C — no R allocation overhead.
  if (isTRUE(.am_af_ok()))
    return(.Call("am_af_kernel_bridge", X, Y, kernel,
                 as.double(sigma), as.integer(degree), as.double(coef),
                 PACKAGE = "amatrix.arrayfire"))
  # MLX: GPU GEMM + R-level transforms (cheap for linear/poly/cosine; heavier for rbf/lap)
  if (isTRUE(.am_mlx_ok())) {
    G     <- .Call("amatrix_mlx_tcrossprod_bridge", X, Y, PACKAGE = "amatrix.mlx")
    Y_eff <- if (is.null(Y)) X else Y
    return(switch(kernel,
      linear     = G,
      polynomial = (coef + G)^degree,
      cosine     = {
        nx <- sqrt(rowSums(X^2)); ny <- if (is.null(Y)) nx else sqrt(rowSums(Y_eff^2))
        G / pmax(outer(nx, ny), .Machine$double.eps)
      },
      rbf        = {
        nx <- rowSums(X^2); ny <- if (is.null(Y)) nx else rowSums(Y_eff^2)
        D_sq <- pmax(outer(nx, ny, "+") - 2 * G, 0)
        if (is.null(Y)) diag(D_sq) <- 0
        exp(-D_sq / (2 * sigma^2))
      },
      laplacian  = {
        nx <- rowSums(X^2); ny <- if (is.null(Y)) nx else rowSums(Y_eff^2)
        D_sq <- pmax(outer(nx, ny, "+") - 2 * G, 0)
        if (is.null(Y)) diag(D_sq) <- 0
        exp(-sqrt(D_sq) / sigma)
      }
    ))
  }
  # CPU fallback
  Y_eff <- if (is.null(Y)) X else Y
  G <- am_tcrossprod(X, Y)

  switch(kernel,
    linear     = G,
    polynomial = (coef + G)^degree,
    cosine     = {
      nx <- sqrt(rowSums(X^2))
      ny <- if (is.null(Y)) nx else sqrt(rowSums(Y_eff^2))
      G / pmax(outer(nx, ny), .Machine$double.eps)
    },
    rbf        = {
      nx   <- rowSums(X^2); ny <- if (is.null(Y)) nx else rowSums(Y_eff^2)
      D_sq <- pmax(outer(nx, ny, "+") - 2 * G, 0)
      if (is.null(Y)) diag(D_sq) <- 0
      exp(-D_sq / (2 * sigma^2))
    },
    laplacian  = {
      nx   <- rowSums(X^2); ny <- if (is.null(Y)) nx else rowSums(Y_eff^2)
      D_sq <- pmax(outer(nx, ny, "+") - 2 * G, 0)
      if (is.null(Y)) diag(D_sq) <- 0
      exp(-sqrt(D_sq) / sigma)
    }
  )
}

# Tiled pairwise distance for large n (avoids GPU OOM on n > 50k).
# Processes row-blocks of X (and Y) independently, assembling the host result
# block by block.  Exploits symmetry when Y = NULL to halve the GEMM count.
.dist_matrix_tiled <- function(X, Y, method, tile_size) {
  m <- nrow(X)
  symmetric <- is.null(Y)
  Y_eff <- if (symmetric) X else Y
  n <- nrow(Y_eff)

  result <- matrix(0.0, nrow = m, ncol = n)

  i_breaks <- c(seq(1L, m, by = tile_size), m + 1L)
  j_breaks <- c(seq(1L, n, by = tile_size), n + 1L)
  ni <- length(i_breaks) - 1L
  nj <- length(j_breaks) - 1L

  for (ii in seq_len(ni)) {
    i0 <- i_breaks[ii]; i1 <- i_breaks[ii + 1L] - 1L
    Xi <- X[i0:i1, , drop = FALSE]

    # Symmetric: compute only lower-triangular blocks, mirror to upper.
    j_end <- if (symmetric) ii else nj
    for (jj in seq_len(j_end)) {
      j0 <- j_breaks[jj]; j1 <- j_breaks[jj + 1L] - 1L
      Xj <- Y_eff[j0:j1, , drop = FALSE]

      D_sq <- .dist_matrix_sq_gpu(Xi, Xj)
      D_block <- if (method == "euclidean") sqrt(D_sq) else D_sq

      result[i0:i1, j0:j1] <- D_block
      if (symmetric && ii != jj) result[j0:j1, i0:i1] <- t(D_block)
    }
  }

  if (symmetric) diag(result) <- 0
  result
}

#' GPU-accelerated pairwise distance matrix
#'
#' Computes the pairwise distance matrix between rows of \code{X} and \code{Y}.
#' The dominant cost (row inner-products via am_tcrossprod) is dispatched to the
#' active GPU backend (ArrayFire or MLX); norm computation and final transforms
#' run on CPU where they are O(mp + np) — negligible versus the O(mnp) GEMM.
#'
#' @param X Numeric matrix or \code{adgeMatrix}, shape [m, p].
#' @param Y Numeric matrix or \code{adgeMatrix}, shape [n, p], or \code{NULL}
#'   to compute pairwise distances within \code{X} (returns [m, m] matrix).
#' @param method One of \code{"euclidean"} (default), \code{"sqeuclidean"},
#'   or \code{"cosine"}.
#' @param tile_size Integer row-block size for tiled computation, or \code{NULL}
#'   (default) to auto-tile when \code{nrow(X) > 50000} (self-distance only).
#'   Set explicitly to process any size in row-blocks; useful when GPU memory
#'   is limited.  Not supported for \code{method = "cosine"}.
#' @return Numeric matrix [m, n] of pairwise distances.
#' @seealso \code{\link{kernel_matrix}}
#' @export
dist_matrix <- function(X, Y = NULL,
                    method = c("euclidean", "sqeuclidean", "cosine"),
                    tile_size = NULL) {
  method <- match.arg(method)
  X_mat <- .am_as_double_matrix(X)
  Y_mat <- if (!is.null(Y)) .am_as_double_matrix(Y) else NULL

  # Auto-tile for large self-distance to prevent GPU OOM
  if (is.null(tile_size) && is.null(Y_mat) && nrow(X_mat) > 50000L)
    tile_size <- 10000L

  if (!is.null(tile_size) && method != "cosine") {
    return(.dist_matrix_tiled(X_mat, Y_mat, method, as.integer(tile_size)))
  }

  if (method == "cosine")
    return(.am_kernel_gpu(X_mat, Y_mat, "cosine", 1.0, 2L, 0.0))

  D_sq <- .dist_matrix_sq_gpu(X_mat, Y_mat)
  if (is.null(Y_mat)) diag(D_sq) <- 0   # fix float32/float64 diagonal mismatch
  if (method == "sqeuclidean") return(D_sq)
  sqrt(D_sq)
}

#' GPU-accelerated pairwise kernel matrix
#'
#' Computes the pairwise kernel matrix between rows of \code{X} and \code{Y}.
#' The expensive am_tcrossprod is GPU-dispatched; element-wise transforms (exp,
#' sqrt, pow) run on CPU.
#'
#' Kernels:
#' \describe{
#'   \item{linear}{k(x,y) = x·y}
#'   \item{rbf}{k(x,y) = exp(-||x-y||² / (2σ²))}
#'   \item{polynomial}{k(x,y) = (coef + x·y)^degree}
#'   \item{cosine}{k(x,y) = x·y / (||x|| ||y||)}
#'   \item{laplacian}{k(x,y) = exp(-||x-y|| / σ)}
#' }
#'
#' @param X Numeric matrix or \code{adgeMatrix}, shape [m, p].
#' @param Y Numeric matrix or \code{adgeMatrix}, shape [n, p], or \code{NULL}.
#' @param kernel Kernel type string (see Details).
#' @param sigma Bandwidth for \code{"rbf"} and \code{"laplacian"}.
#' @param degree Polynomial degree for \code{"polynomial"}.
#' @param coef Constant term for \code{"polynomial"}: (coef + x·y)^degree.
#' @return Numeric matrix [m, n] of kernel values.
#' @seealso \code{\link{dist_matrix}}
#' @export
kernel_matrix <- function(X, Y = NULL,
                      kernel = c("linear", "rbf", "polynomial",
                                 "cosine", "laplacian"),
                      sigma = 1.0, degree = 2L, coef = 0.0) {
  kernel <- match.arg(kernel)
  X_mat  <- .am_as_double_matrix(X)
  Y_mat  <- if (!is.null(Y)) .am_as_double_matrix(Y) else NULL
  Y_eff  <- if (is.null(Y_mat)) X_mat else Y_mat

  .am_kernel_gpu(X_mat, Y_mat, kernel,
                 sigma = sigma, degree = degree, coef = coef)
}

# ---------------------------------------------------------------------------
# Kronecker product
# ---------------------------------------------------------------------------


# norm: matrix and vector norms.
#
# Supported types (matching base::norm):
#   "1" - max absolute column sum (1-norm)
#   "I" - max absolute row sum (infinity-norm)
#   "F" - Frobenius norm  (default for matrices; Euclidean for vectors)
#   "M" - max absolute entry
#   "2" - spectral norm (largest singular value); uses rsvd for large matrices
#
# For numeric vector inputs, "2" and "F" return Euclidean norm,
# "1" returns sum(|x|), "I"/"M" return max(|x|).
#
# Implementation uses S4 methods to override Matrix package's ANY dispatch.

.norm_type <- function(type) {
  toupper(match.arg(type, c("1", "I", "F", "f", "M", "m", "2")))
}

# S4 method for adgeMatrix (primary dispatch)
setMethod("norm", "adgeMatrix", function(x, type = "F", ...) {
  type <- .norm_type(type)
  if (type == "2") {
    sv <- rsvd(x, k = 1L)
    return(sv$d[[1L]])
  }
  base::norm(.am_as_double_matrix(x), type = type)
})

# S4 method for plain base R matrix
setMethod("norm", "matrix", function(x, type = "F", ...) {
  type <- .norm_type(type)
  if (type == "2") {
    sv <- rsvd(as_adgeMatrix(x), k = 1L)
    return(sv$d[[1L]])
  }
  base::norm(x, type = type)
})

# S4 method for numeric vector: vector norms
setMethod("norm", "numeric", function(x, type = "F", ...) {
  type <- .norm_type(type)
  switch(type,
    "2" = ,
    "F" = sqrt(sum(x^2)),
    "1" = sum(abs(x)),
    "I" = ,
    "M" = max(abs(x))
  )
})

# ---------------------------------------------------------------------------
# Batch factorization helpers
# ---------------------------------------------------------------------------

# Normalise a batch argument to a list of plain R matrices.
# Accepts: list of matrices, or a 3-D array dim c(n, n, B).
.am_batch_to_list <- function(A, arg_name = "A") {
  if (is.list(A)) {
    lapply(A, function(a) {
      m <- as.matrix(a)
      if (!is.double(m)) storage.mode(m) <- "double"
      m
    })
  } else if (is.array(A) && length(dim(A)) == 3L) {
    B <- dim(A)[[3L]]
    lapply(seq_len(B), function(b) {
      m <- A[, , b, drop = FALSE]
      dim(m) <- dim(m)[1:2]
      if (!is.double(m)) storage.mode(m) <- "double"
      m
    })
  } else {
    stop(arg_name, " must be a list of matrices or a 3-D array [n, n, B]",
         call. = FALSE)
  }
}

#' Batch Cholesky factorization
#'
#' Factorize B symmetric positive-definite matrices in parallel.  Each matrix
#' is dispatched through the same backend as \code{\link{chol_factor}}, so MLX
#' GPU acceleration applies to every element when available.
#'
#' @param A A list of square numeric matrices, or a 3-D array \code{[n, n, B]}.
#' @return A list of \code{amChol} objects, one per input matrix.
#' @seealso \code{\link{chol_factor}}, \code{\link{batch_solve}}
#' @export
batch_chol <- function(A) {
  mats <- .am_batch_to_list(A, "A")
  lapply(mats, function(m) {
    X <- adgeMatrix(m)
    chol_factor(X)
  })
}

#' Batch triangular solve
#'
#' Solve B linear systems \code{A_b x_b = B_b} where each \code{A_b} is
#' represented by its Cholesky factor from \code{\link{batch_chol}}.
#'
#' @param Ls A list of \code{amChol} objects (output of \code{batch_chol}).
#' @param B A list of right-hand-side matrices/vectors, or a 3-D array
#'   \code{[n, k, B]}.  Length / third dimension must match \code{Ls}.
#' @return A list of solution matrices (or vectors when each rhs is a vector).
#' @seealso \code{\link{batch_chol}}, \code{\link{chol_solve}}
#' @export
batch_solve <- function(Ls, B) {
  if (!is.list(Ls) || !all(vapply(Ls, inherits, logical(1L), "amChol"))) {
    stop("Ls must be a list of amChol objects from batch_chol()", call. = FALSE)
  }

  rhs_list <- if (is.list(B)) {
    B
  } else if (is.array(B) && length(dim(B)) == 3L) {
    nb <- dim(B)[[3L]]
    lapply(seq_len(nb), function(b) B[, , b, drop = TRUE])
  } else {
    stop("B must be a list or a 3-D array [n, k, B]", call. = FALSE)
  }

  if (length(Ls) != length(rhs_list)) {
    stop("Ls and B must have the same batch size", call. = FALSE)
  }

  Map(chol_solve, Ls, rhs_list)
}

#' Eager Kronecker product
#'
#' Computes \code{A ⊗ B} and returns the result as an \code{adgeMatrix}.
#' Accepts plain matrices or any \code{aMatrix} subclass.
#' For a lazy variant that avoids forming the full product see
#' \code{\link{kron_matrix}}.
#'
#' @param A,B Matrices or \code{aMatrix} objects.
#' @return An \code{adgeMatrix} of dimension \code{(nrow(A)*nrow(B)) x (ncol(A)*ncol(B))}.
#' @seealso \code{\link{kron_matrix}}
#' @export
kron <- function(A, B) {
  A_mat <- .am_as_double_matrix(A)
  B_mat <- .am_as_double_matrix(B)
  adgeMatrix(base::kronecker(A_mat, B_mat))
}

#' Batch crossproduct
#'
#' Compute \code{t(A_b) \%*\% A_b} for each matrix in a batch.
#'
#' @param A A list of numeric matrices, or a 3-D array \code{[n, p, B]}.
#' @return A list of \code{p x p} crossproduct matrices.
#' @seealso \code{\link{batch_chol}}
#' @export
batch_crossprod <- function(A) {
  mats <- .am_batch_to_list(A, "A")
  lapply(mats, function(m) crossprod(m))
}
