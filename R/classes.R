.validate_amatrix_slots <- function(object) {
  if (!identical(length(object@preferred_backend), 1L) || is.na(object@preferred_backend)) {
    return("preferred_backend must be a single non-missing string")
  }

  if (!identical(length(object@policy), 1L) || is.na(object@policy)) {
    return("policy must be a single non-missing string")
  }

  if (!(object@policy %in% .amatrix_valid_policies)) {
    return(sprintf("policy must be one of: %s", paste(.amatrix_valid_policies, collapse = ", ")))
  }

  if (!identical(length(object@precision), 1L) || is.na(object@precision)) {
    return("precision must be a single non-missing string")
  }

  if (!(object@precision %in% .amatrix_valid_precisions)) {
    return(sprintf("precision must be one of: %s", paste(.amatrix_valid_precisions, collapse = ", ")))
  }

  if (!identical(length(object@object_id), 1L) || is.na(object@object_id) || !nzchar(object@object_id)) {
    return("object_id must be a single non-missing string")
  }

  TRUE
}

#' Virtual base class for backend-aware matrices
#'
#' \code{aMatrix} is the abstract base from which all concrete amatrix
#' classes inherit. It carries backend-dispatch metadata that controls
#' which compute backend (CPU, GPU, etc.) is used for operations on the
#' matrix.
#'
#' @slot preferred_backend Single string naming the preferred compute
#'   backend; one of \code{"cpu"}, \code{"mlx"}, \code{"metal"},
#'   or \code{"arrayfire"}.
#' @slot policy Single string controlling dispatch policy; one of
#'   \code{"auto"}, \code{"cpu"}, \code{"mlx"}, \code{"metal"},
#'   or \code{"arrayfire"}.
#' @slot precision Single string; either \code{"strict"} (double
#'   precision, exact results) or \code{"fast"} (backend may use
#'   lower precision).
#' @slot object_id Non-empty string uniquely identifying this object
#'   within the session; used for caching and residency tracking.
#' @slot src_id String recording the \code{object_id} of the object
#'   this was derived from, or \code{""} for originals.
#' @slot finalizer_env Environment used to manage GPU-resident memory
#'   and deferred host-copy state.
#'
#' @exportClass aMatrix
setClass(
  "aMatrix",
  contains = "VIRTUAL",
  slots = c(
    preferred_backend = "character",
    policy = "character",
    precision = "character",
    object_id = "character",
    src_id = "character",
    finalizer_env = "environment"
  ),
  prototype = list(
    preferred_backend = "cpu",
    policy = "auto",
    precision = "strict",
    object_id = "",
    src_id = "",
    finalizer_env = new.env(parent = emptyenv())
  ),
  validity = .validate_amatrix_slots
)

#' Dense general matrix with backend-dispatch metadata
#'
#' \code{adgeMatrix} extends both \code{aMatrix} and
#' \code{Matrix::dgeMatrix}, adding backend-dispatch slots to a
#' column-major dense double-precision matrix. All arithmetic
#' generics dispatch through the amatrix backend system rather than
#' directly to BLAS.
#'
#' @exportClass adgeMatrix
#' @seealso \code{\link{adgeMatrix}} for the user-facing constructor,
#'   \code{\link{adgCMatrix}} for the sparse counterpart
setClass(
  "adgeMatrix",
  contains = c("aMatrix", "dgeMatrix")
)

#' Sparse column-compressed matrix with backend-dispatch metadata
#'
#' \code{adgCMatrix} extends both \code{aMatrix} and
#' \code{Matrix::dgCMatrix}, adding backend-dispatch slots to a
#' compressed-column sparse double-precision matrix.
#'
#' @exportClass adgCMatrix
#' @seealso \code{\link{adgCMatrix}} for the user-facing constructor,
#'   \code{\link{adgeMatrix-class}} for the dense counterpart
setClass(
  "adgCMatrix",
  contains = c("aMatrix", "dgCMatrix")
)

#' Dense logical matrix with backend-dispatch metadata
#'
#' \code{adlgeMatrix} extends both \code{aMatrix} and
#' \code{Matrix::lgeMatrix}, adding backend-dispatch slots to a
#' column-major dense logical matrix.
#'
#' @exportClass adlgeMatrix
setClass(
  "adlgeMatrix",
  contains = c("aMatrix", "lgeMatrix")
)

#' Sparse logical matrix with backend-dispatch metadata
#'
#' \code{adlgCMatrix} extends both \code{aMatrix} and
#' \code{Matrix::lgCMatrix}, adding backend-dispatch slots to a
#' compressed-column sparse logical matrix.
#'
#' @exportClass adlgCMatrix
setClass(
  "adlgCMatrix",
  contains = c("aMatrix", "lgCMatrix")
)

#' Lazy transpose view of an adgeMatrix
#'
#' \code{aTransposeView} is a zero-copy structural view representing
#' the transpose of an \code{adgeMatrix}. It carries no independent
#' dense host storage; the underlying data lives in \code{source}.
#' The transposed matrix is materialized on demand via
#' \code{as.matrix()} or \code{amatrix_materialize_host()}.
#'
#' @slot source The originating \code{adgeMatrix}; kept alive by this
#'   reference.
#' @slot Dim Integer vector of length 2 giving the transposed
#'   dimensions \code{c(ncol_src, nrow_src)}.
#' @slot Dimnames List of length 2 with transposed dimnames.
#'
#' @exportClass aTransposeView
setClass(
  "aTransposeView",
  contains = "aMatrix",
  slots = c(
    source   = "adgeMatrix",  # reference to the original; keeps it alive
    Dim      = "integer",     # transposed dims [ncol_src, nrow_src]
    Dimnames = "list"         # transposed dimnames
  )
)

#' @noRd
setMethod("show", "adgeMatrix", function(object) {
  if (.amatrix_is_dead_deferred(object)) {
    .amatrix_dead_deferred_error()
  }

  cat(sprintf(
    "An amatrix dense matrix [%s|policy=%s|precision=%s]\n",
    object@preferred_backend,
    object@policy,
    object@precision
  ))
  # Deferred objects: materialize before callNextMethod() reads @x
  if (isTRUE(object@finalizer_env$host_deferred)) {
    mat <- amatrix_materialize_dense(object)
    show(mat)
  } else {
    callNextMethod()
  }
})

#' @noRd
setMethod("show", "adgCMatrix", function(object) {
  cat(sprintf(
    "An amatrix sparse matrix [%s|policy=%s|precision=%s]\n",
    object@preferred_backend,
    object@policy,
    object@precision
  ))
  callNextMethod()
})

#' @noRd
setMethod("show", "adlgeMatrix", function(object) {
  cat(sprintf(
    "An amatrix dense logical matrix [%s|policy=%s|precision=%s]\n",
    object@preferred_backend,
    object@policy,
    object@precision
  ))
  callNextMethod()
})

#' @noRd
setMethod("show", "adlgCMatrix", function(object) {
  cat(sprintf(
    "An amatrix sparse logical matrix [%s|policy=%s|precision=%s]\n",
    object@preferred_backend,
    object@policy,
    object@precision
  ))
  callNextMethod()
})
