.amatrix_prepare_binary_operand <- function(x, preferred_backend = "auto", precision = NULL, policy = NULL) {
  if (inherits(x, "adgeMatrix") || inherits(x, "adgCMatrix")) {
    return(x)
  }

  precision <- if (is.null(precision)) amatrix_default_precision() else precision
  policy <- if (is.null(policy)) amatrix_default_policy() else policy

  if (inherits(x, "sparseMatrix")) {
    return(as_adgCMatrix(
      x,
      preferred_backend = preferred_backend,
      policy = policy,
      precision = precision
    ))
  }

  if (is.matrix(x) || inherits(x, "denseMatrix") || inherits(x, "dgeMatrix")) {
    return(as_adgeMatrix(
      as.matrix(x),
      preferred_backend = preferred_backend,
      policy = policy,
      precision = precision
    ))
  }

  x
}

.amatrix_select_binary_resident_backend <- function(x, y, op, backend = "auto") {
  if (!is.null(backend) && !identical(backend, "auto")) {
    return(backend)
  }

  if (inherits(x, "aMatrix")) {
    be <- amatrix_resident_backend_for(x, op = op, y = y)
    if (!is.null(be)) {
      return(be)
    }
  }

  # Sparse-right matmul is implemented through a dedicated resident path.
  if (identical(op, "matmul") &&
      inherits(y, "adgCMatrix") &&
      inherits(x, "adgeMatrix")) {
    candidates <- .amatrix_resident_backend_candidates(y, op = op)
    for (backend_name in candidates) {
      backend_obj <- tryCatch(.amatrix_get_backend(backend_name), error = function(e) NULL)
      if (is.null(backend_obj) ||
          !.amatrix_backend_available_safe(backend_obj) ||
          !.amatrix_backend_residency_capable(backend_obj)) {
        next
      }
      if (!(x@precision %in% unique(backend_obj$precision_modes())) ||
          !(y@precision %in% unique(backend_obj$precision_modes()))) {
        next
      }
      if (is.function(backend_obj$dense_sparse_matmul_resident_key)) {
        return(backend_name)
      }
    }
  }

  if (inherits(y, "aMatrix")) {
    amatrix_resident_backend_for(y, op = op, y = x)
  } else {
    NULL
  }
}

#' Prepare operands for a repeated matrix product
#'
#' Converts inputs to \code{amatrix} wrappers when needed, chooses a
#' residency-capable accelerator backend in automatic mode, and binds the
#' operands so repeated products reuse the resident fast path.
#'
#' @param x Left operand.
#' @param y Right operand.
#' @param op Product primitive: \code{"matmul"}, \code{"crossprod"}, or
#'   \code{"tcrossprod"}.
#' @param backend Backend name or \code{"auto"}.
#' @param precision Precision to use when wrapping base matrices.
#' @param policy Policy to use when wrapping base matrices.
#' @return A list with elements \code{x}, \code{y}, and \code{backend}.
#' @export
amatrix_prepare_operands <- function(
  x,
  y,
  op = c("matmul", "crossprod", "tcrossprod"),
  backend = "auto",
  precision = amatrix_default_precision(),
  policy = amatrix_default_policy()
) {
  precision_missing <- missing(precision)
  op <- match.arg(op)
  precision <- .amatrix_resolve_backend_precision(
    backend,
    precision,
    precision_missing = precision_missing
  )

  x_prepared <- .amatrix_prepare_binary_operand(
    x,
    preferred_backend = backend,
    precision = precision,
    policy = policy
  )
  y_prepared <- .amatrix_prepare_binary_operand(
    y,
    preferred_backend = backend,
    precision = precision,
    policy = policy
  )

  backend_name <- .amatrix_select_binary_resident_backend(x_prepared, y_prepared, op = op, backend = backend)

  if (is.null(backend_name) || identical(backend_name, "cpu")) {
    return(list(x = x_prepared, y = y_prepared, backend = NULL))
  }
  if (inherits(x_prepared, "aMatrix")) {
    .amatrix_check_backend_precision(backend_name, x_prepared@precision)
  }
  if (inherits(y_prepared, "aMatrix")) {
    .amatrix_check_backend_precision(backend_name, y_prepared@precision)
  }

  x_bound <- if (inherits(x_prepared, "aMatrix")) {
    amatrix_bind_resident(x_prepared, backend = backend_name, op = op, y = y_prepared)
  } else {
    x_prepared
  }

  y_bound <- if (inherits(y_prepared, "aMatrix")) {
    amatrix_bind_resident(y_prepared, backend = backend_name)
  } else {
    y_prepared
  }

  list(x = x_bound, y = y_bound, backend = backend_name)
}
