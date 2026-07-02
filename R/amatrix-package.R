#' amatrix package
#'
#' Package-level roxygen directives for namespace generation.
#'
#' @section Startup message:
#' On attach, amatrix prints a one-line note listing any installed GPU
#' backend packages (see \code{\link{amatrix_use_gpu}}). Pure-CPU sessions
#' with no backend packages installed are silent. To suppress the note
#' explicitly, set \code{options(amatrix.quiet_startup = TRUE)} or the
#' environment variable \code{AMATRIX_QUIET=1} (either \code{"1"} or
#' \code{"true"}) before loading the package. The note is also skipped
#' when \code{options(amatrix.optional_backends = FALSE)} disables all
#' optional backends, and any backend disabled via
#' \code{options(amatrix.disable_mlx = TRUE)} (or the analogous
#' \code{amatrix.disable_metal}, \code{amatrix.disable_opencl},
#' \code{amatrix.disable_arrayfire}) is omitted from it.
#'
#' @keywords internal
#' @useDynLib amatrix, .registration = TRUE
#' @importFrom Matrix Matrix
#' @importClassesFrom Matrix dgCMatrix dgeMatrix
#' @importFrom methods as callGeneric callNextMethod is new setAs setClass setGeneric setMethod setOldClass setValidity show slotNames
#' @importFrom stats coef fitted median residuals rnorm
"_PACKAGE"
