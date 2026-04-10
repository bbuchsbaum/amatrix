.amatrix_auto_resident_backend_order <- function(op = NULL, x = NULL) {
  if (inherits(x, "adgCMatrix") && op %in% c("matmul", "crossprod", "tcrossprod")) {
    c("metal", "arrayfire", "mlx", "opencl", "torch")
  } else {
    c("mlx", "metal", "arrayfire", "opencl", "torch")
  }
}

.amatrix_resident_backend_candidates <- function(x, op = NULL) {
  preferred <- .amatrix_backend_preference(x, op = op)
  explicit <- preferred[!(preferred %in% c("", "auto", "cpu"))]

  registered <- setdiff(amatrix_backend_names(), "cpu")
  auto_priority <- .amatrix_auto_resident_backend_order(op = op, x = x)
  auto_candidates <- unique(c(
    intersect(auto_priority, registered),
    setdiff(registered, auto_priority)
  ))

  unique(c(explicit, auto_candidates))
}

#' Choose a residency-capable accelerator backend for a hot path
#'
#' Returns the first available non-CPU backend that can keep \code{x} resident
#' for the requested operation. This is intended for package authors who want
#' repeated work to stay on the fastest available accelerator without hardcoding
#' backend names such as \code{"metal"} or \code{"mlx"}.
#'
#' @param x An \code{aMatrix}.
#' @param op Optional operation name such as \code{"matmul"}.
#' @param y Optional rhs object used when checking resident-op support.
#' @return A backend name, or \code{NULL} when no residency-capable accelerator
#'   is available for the requested operation.
#' @export
amatrix_resident_backend_for <- function(x, op = NULL, y = NULL) {
  stopifnot(inherits(x, "aMatrix"))

  pinned <- .amatrix_live_resident_backend(x)
  if (!is.null(pinned)) {
    backend <- tryCatch(.amatrix_get_backend(pinned), error = function(e) NULL)
    if (!is.null(backend) &&
        .amatrix_backend_residency_capable(backend) &&
        (is.null(op) || .amatrix_backend_supports_resident_op(backend, op, x = x, y = y))) {
      return(pinned)
    }
  }

  candidates <- .amatrix_resident_backend_candidates(x, op = op)
  if (length(candidates) == 0L) {
    return(NULL)
  }

  for (backend_name in candidates) {
    backend <- tryCatch(.amatrix_get_backend(backend_name), error = function(e) NULL)
    if (is.null(backend) ||
        !isTRUE(backend$available()) ||
        !.amatrix_backend_residency_capable(backend) ||
        !(x@precision %in% unique(backend$precision_modes()))) {
      next
    }
    if (is.null(op) || .amatrix_backend_supports_resident_op(backend, op, x = x, y = y)) {
      return(backend_name)
    }
  }

  NULL
}
