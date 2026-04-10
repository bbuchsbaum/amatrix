#' Bind an amatrix object to resident backend storage
#'
#' Upload a dense or sparse matrix to a residency-capable backend and return the
#' corresponding \code{aMatrix} object with a live resident binding. This is
#' primarily useful for repeated GPU work where paying the upload cost once is
#' preferable to relying on cold-path dispatch.
#'
#' @param x An \code{adgeMatrix}, \code{adgCMatrix}, base matrix, or sparse
#'   Matrix object.
#' @param backend Backend name, \code{"auto"}, or \code{NULL}. When left
#'   \code{NULL} or set to \code{"auto"}, \code{amatrix} picks the first
#'   residency-capable accelerator backend that supports the requested resident
#'   operation.
#' @param op Optional operation name such as \code{"matmul"} used when
#'   selecting an automatic resident backend.
#' @param y Optional rhs object used when checking resident-op support for
#'   automatic backend selection.
#' @return An \code{adgeMatrix} or \code{adgCMatrix} with a live resident
#'   binding on \code{backend}. When no suitable accelerator backend is
#'   available in automatic mode, returns \code{x} unchanged.
#' @export
amatrix_bind_resident <- function(x, backend = NULL, op = NULL, y = NULL) {
  if (inherits(x, "adgeMatrix") || inherits(x, "adgCMatrix")) {
    obj <- x
  } else if (is.matrix(x) || inherits(x, "dgeMatrix")) {
    preferred_backend <- if (!is.null(backend)) backend else amatrix_default_policy()
    precision <- .amatrix_resolve_backend_precision(
      preferred_backend,
      amatrix_default_precision(),
      precision_missing = TRUE
    )
    obj <- as_adgeMatrix(x, preferred_backend = preferred_backend, precision = precision)
  } else if (inherits(x, "sparseMatrix")) {
    preferred_backend <- if (!is.null(backend)) backend else amatrix_default_policy()
    precision <- .amatrix_resolve_backend_precision(
      preferred_backend,
      amatrix_default_precision(),
      precision_missing = TRUE
    )
    obj <- as_adgCMatrix(x, preferred_backend = preferred_backend, precision = precision)
  } else {
    stop("x must be an adgeMatrix, adgCMatrix, matrix, or sparse Matrix", call. = FALSE)
  }

  backend_name <- backend
  auto_backend <- is.null(backend_name) || identical(backend_name, "auto")
  if (auto_backend) {
    backend_name <- amatrix_resident_backend_for(obj, op = op, y = y)
    if (is.null(backend_name)) {
      return(obj)
    }
  }
  if (is.null(backend_name) && inherits(obj, "aMatrix") && nzchar(obj@preferred_backend)) {
    backend_name <- obj@preferred_backend
  }
  if (is.null(backend_name) || !nzchar(backend_name)) {
    backend_name <- amatrix_default_policy()
  }

  backend_obj <- .amatrix_get_backend(backend_name)
  if (!isTRUE(backend_obj$available())) {
    stop(sprintf("backend '%s' is not available", backend_name), call. = FALSE)
  }
  if (!.amatrix_backend_residency_capable(backend_obj)) {
    stop(sprintf("backend '%s' does not support residency", backend_name), call. = FALSE)
  }
  .amatrix_check_backend_precision(backend_name, obj@precision)

  existing <- .amatrix_resident_entry(obj)
  if (!is.null(existing) &&
      identical(existing$backend, backend_name) &&
      .amatrix_backend_has_resident_key(backend_obj, existing$resident_key, sparse = isTRUE(existing$sparse))) {
    return(obj)
  }

  if (!is.null(existing)) {
    old_backend <- tryCatch(.amatrix_get_backend(existing$backend), error = function(e) NULL)
    if (!is.null(old_backend) && .amatrix_backend_residency_capable(old_backend)) {
      if (isTRUE(existing$sparse) && is.function(old_backend$sparse_resident_drop)) {
        try(old_backend$sparse_resident_drop(existing$resident_key), silent = TRUE)
      } else if (is.function(old_backend$resident_drop)) {
        try(old_backend$resident_drop(existing$resident_key), silent = TRUE)
      }
    }
    .amatrix_drop_resident_binding(obj)
  }

  resident_key <- .amatrix_next_resident_key(backend_name)

  if (inherits(obj, "adgCMatrix")) {
    if (!is.function(backend_obj$sparse_resident_store)) {
      stop(sprintf("backend '%s' does not support sparse residency", backend_name), call. = FALSE)
    }
    backend_obj$sparse_resident_store(resident_key, amatrix_materialize_host(obj))
    return(.amatrix_bind_resident(obj, backend_name, resident_key, sparse = TRUE))
  }

  backend_obj$resident_store(resident_key, amatrix_materialize_host(obj))
  .amatrix_bind_resident(obj, backend_name, resident_key)
}
