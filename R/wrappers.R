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

# solve: output is a matrix; store at out_key and bind resident.
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

am_matmul <- function(x, y) {
  choice <- .amatrix_backend_for(x, "matmul", y = y)
  resident <- .amatrix_try_resident_matmul(x, y, choice$name)
  if (!is.null(resident)) {
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

am_row_sums <- function(x, na.rm = FALSE, dims = 1L) {
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

am_col_sums <- function(x, na.rm = FALSE, dims = 1L) {
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
  choice <- .amatrix_backend_for(a, "solve", y = b_arg)
  resident <- .amatrix_try_resident_solve(a, b_arg, choice$name)
  if (!is.null(resident)) {
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

  .amatrix_rewrap_value(
    a,
    amatrix_dispatch_op(
      x = a,
      op = "solve",
      method = "solve",
      y = b_arg,
      args = list(b = .amatrix_host_arg(b_arg), ...),
      fallback = function() base::solve(as.matrix(amatrix_materialize_host(a)), as.matrix(.amatrix_host_arg(b_arg)), ...)
    )
  )
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

am_eigen <- function(x, symmetric, only.values = FALSE, EISPACK = FALSE) {
  amatrix_dispatch_op(
    x = x,
    op = "eigen",
    method = "eigen",
    args = list(symmetric = symmetric, only.values = only.values, EISPACK = EISPACK),
    fallback = function() base::eigen(as.matrix(amatrix_materialize_host(x)), symmetric = symmetric, only.values = only.values, EISPACK = EISPACK)
  )
}

am_diag <- function(x, nrow, ncol, names = TRUE) {
  amatrix_dispatch_op(
    x = x,
    op = "diag",
    method = "diag",
    args = list(nrow = nrow, ncol = ncol, names = names),
    fallback = function() base::diag(as.matrix(amatrix_materialize_host(x)), nrow = nrow, ncol = ncol, names = names)
  )
}

am_ewise <- function(op, e1, e2 = NULL) {
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
