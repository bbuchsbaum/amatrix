# Public wrapper around the internal .amatrix_release_resident() helper
# (defined in R/residency.R). Kept in its own file so the exported user-facing
# API is separate from the residency internals.

#' Release GPU-resident data held by an amatrix object
#'
#' Frees any device-resident buffer associated with \code{x} and drops its
#' residency-registry binding, leaving the host copy as the authoritative
#' storage. This gives long-lived GPU pipelines explicit control over device
#' memory instead of waiting for garbage collection to reclaim resident
#' handles.
#'
#' The object remains fully usable afterwards: its data is served from the host
#' copy and is re-uploaded to the device on the next GPU operation if needed. On
#' CPU-only sessions, or for any object that currently holds no device buffer,
#' this is a safe no-op.
#'
#' @param x An \code{\linkS4class{aMatrix}} object (for example an
#'   \code{\linkS4class{adgeMatrix}} or \code{\linkS4class{adgCMatrix}}).
#'   Non-amatrix inputs are ignored.
#'
#' @return Invisibly, \code{TRUE} if a resident binding was released, and
#'   \code{FALSE} otherwise (including the CPU-only no-op case).
#'
#' @examples
#' A <- adgeMatrix(matrix(1:6, 2, 3))
#' # On a CPU-only session there is no device buffer, so this is a no-op:
#' released <- amatrix_release_resident(A)
#' released
#'
#' @export
amatrix_release_resident <- function(x) {
  invisible(.amatrix_release_resident(x))
}
