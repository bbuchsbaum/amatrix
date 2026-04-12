.amatrix_build_dgeMatrix <- function(x, dn = base::dimnames(x)) {
  stopifnot(is.matrix(x))
  storage.mode(x) <- "double"
  if (is.null(dn)) {
    dn <- vector("list", 2L)
  }
  new(
    "dgeMatrix",
    x = as.double(x),
    Dim = as.integer(dim(x)),
    Dimnames = dn,
    factors = list()
  )
}

.amatrix_dense_slot_matrix <- function(x) {
  stopifnot(inherits(x, "denseMatrix") || inherits(x, "adgeMatrix"))
  out <- x@x
  dim(out) <- as.integer(x@Dim)
  dimnames(out) <- x@Dimnames
  storage.mode(out) <- "double"
  out
}

.amatrix_new_dense <- function(
  data,
  dim,
  dimnames = NULL,
  factors = list(),
  preferred_backend = "cpu",
  policy = amatrix_default_policy(),
  precision = amatrix_default_precision(),
  src_id = ""
) {
  if (!is.double(data)) {
    storage.mode(data) <- "double"
  }
  if (is.null(dimnames)) {
    dimnames <- vector("list", 2L)
  }
  object_id <- .amatrix_next_object_id()
  new(
    "adgeMatrix",
    x = as.double(data),
    Dim = as.integer(dim),
    Dimnames = dimnames,
    factors = factors,
    preferred_backend = preferred_backend,
    policy = policy,
    precision = precision,
    object_id = object_id,
    src_id = src_id,
    finalizer_env = .amatrix_make_finalizer_env(object_id)
  )
}

.amatrix_dense_base <- function(x) {
  if (inherits(x, "dgeMatrix")) {
    return(x)
  }
  if (inherits(x, "denseMatrix")) {
    return(.amatrix_build_dgeMatrix(as.matrix(x), dn = base::dimnames(x)))
  }
  if (is.matrix(x)) {
    return(.amatrix_build_dgeMatrix(x))
  }
  stop("x must be a base matrix or dgeMatrix")
}

.amatrix_sparse_base <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    return(x)
  }
  if (inherits(x, "sparseMatrix")) {
    base <- x
    if (!inherits(base, "generalMatrix")) {
      base <- as(base, "generalMatrix")
    }
    return(as(base, "dgCMatrix"))
  }
  if (is.matrix(x)) {
    base <- Matrix(x, sparse = TRUE)
    if (!inherits(base, "generalMatrix")) {
      base <- as(base, "generalMatrix")
    }
    return(as(base, "dgCMatrix"))
  }
  stop("x must be a base matrix or dgCMatrix")
}

#' Construct an adgeMatrix from a matrix or dgeMatrix
#'
#' Wraps a base R matrix or \code{Matrix::dgeMatrix} in an
#' \code{adgeMatrix}, attaching backend-dispatch metadata.
#'
#' @param x A base R \code{matrix}, \code{dgeMatrix}, or any
#'   \code{denseMatrix} coercible to \code{dgeMatrix}.
#' @param preferred_backend Single string; the preferred compute
#'   backend. Defaults to \code{"cpu"}.
#' @param policy Single string; backend dispatch policy. Defaults to
#'   \code{amatrix_default_policy()}.
#' @param precision Single string; either \code{"strict"} or
#'   \code{"fast"}. Defaults to \code{amatrix_default_precision()}.
#' @param src_id String recording the source object identifier.
#'   Pass \code{""} (default) for new objects.
#'
#' @return An \code{adgeMatrix} object with the same data as \code{x}.
#'
#' @keywords internal
new_adgeMatrix <- function(
  x,
  preferred_backend = "cpu",
  policy = amatrix_default_policy(),
  precision = amatrix_default_precision(),
  src_id = ""
) {
  if (inherits(x, "dgeMatrix")) {
    return(.amatrix_new_dense(
      data = x@x,
      dim = x@Dim,
      dimnames = x@Dimnames,
      factors = x@factors,
      preferred_backend = preferred_backend,
      policy = policy,
      precision = precision,
      src_id = src_id
    ))
  }

  if (is.matrix(x)) {
    return(.amatrix_new_dense(
      data = x,
      dim = dim(x),
      dimnames = base::dimnames(x),
      factors = list(),
      preferred_backend = preferred_backend,
      policy = policy,
      precision = precision,
      src_id = src_id
    ))
  }

  base <- .amatrix_dense_base(x)
  .amatrix_new_dense(
    data = base@x,
    dim = base@Dim,
    dimnames = base@Dimnames,
    factors = base@factors,
    preferred_backend = preferred_backend,
    policy = policy,
    precision = precision,
    src_id = src_id
  )
}

#' Construct a deferred adgeMatrix with GPU-only storage
#'
#' Creates an \code{adgeMatrix} whose host \code{@x} slot holds a
#' \code{NaN} sentinel vector. The true data lives only on the device
#' until the first host access, which triggers a transparent download.
#'
#' @param dim Integer vector of length 2 giving \code{c(nrow, ncol)}.
#' @param dimnames List of length 2 with row and column names, or
#'   \code{list(NULL, NULL)}.
#' @param preferred_backend Single string naming the preferred backend.
#' @param policy Single string; backend dispatch policy.
#' @param precision Single string; \code{"strict"} or \code{"fast"}.
#' @param src_id String recording the source object identifier.
#'
#' @return An \code{adgeMatrix} with \code{finalizer_env$host_deferred}
#'   set to \code{TRUE}.
#'
#' @keywords internal
new_adgeMatrix_deferred <- function(
  dim,
  dimnames = list(NULL, NULL),
  preferred_backend = "cpu",
  policy = amatrix_default_policy(),
  precision = amatrix_default_precision(),
  src_id = ""
) {
  n <- as.integer(dim[1L]) * as.integer(dim[2L])
  object_id <- .amatrix_next_object_id()
  fenv <- .amatrix_make_finalizer_env(object_id)
  fenv$host_deferred <- TRUE
  fenv$host_x <- NULL

  new(
    "adgeMatrix",
    x = rep(NaN, n),
    Dim = as.integer(dim),
    Dimnames = if (is.null(dimnames)) list(NULL, NULL) else dimnames,
    factors = list(),
    preferred_backend = preferred_backend,
    policy = policy,
    precision = precision,
    object_id = object_id,
    src_id = src_id,
    finalizer_env = fenv
  )
}

#' Construct an adgCMatrix from a sparse or dense matrix
#'
#' Wraps a \code{Matrix::dgCMatrix} or any sparse or dense matrix
#' in an \code{adgCMatrix}, attaching backend-dispatch metadata.
#'
#' @param x A \code{dgCMatrix}, other \code{sparseMatrix}, or base R
#'   \code{matrix} to convert.
#' @param preferred_backend Single string; the preferred compute
#'   backend. Defaults to \code{"cpu"}.
#' @param policy Single string; backend dispatch policy.
#' @param precision Single string; \code{"strict"} or \code{"fast"}.
#'
#' @return An \code{adgCMatrix} object.
#'
#' @keywords internal
new_adgCMatrix <- function(
  x,
  preferred_backend = "cpu",
  policy = amatrix_default_policy(),
  precision = amatrix_default_precision()
) {
  base <- .amatrix_sparse_base(x)
  object_id <- .amatrix_next_object_id()
  new(
    "adgCMatrix",
    i = base@i,
    p = base@p,
    Dim = base@Dim,
    Dimnames = base@Dimnames,
    x = base@x,
    factors = base@factors,
    preferred_backend = preferred_backend,
    policy = policy,
    precision = precision,
    object_id = object_id,
    finalizer_env = .amatrix_make_finalizer_env(object_id)
  )
}

#' Create a backend-aware dense matrix
#'
#' Converts a base R matrix or \code{Matrix::dgeMatrix} to an
#' \code{adgeMatrix} with the specified backend, policy, and precision.
#' This is the primary user-facing constructor for dense amatrix
#' objects.
#'
#' @param x A base R \code{matrix}, \code{dgeMatrix}, or any
#'   \code{denseMatrix} coercible to \code{dgeMatrix}.
#' @param mode Single string shortcut accepted by
#'   \code{.amatrix_resolve_mode()}; used to set \code{backend},
#'   \code{policy}, and \code{precision} together. Pass \code{NULL}
#'   to use the individual arguments instead.
#' @param backend Alias for \code{preferred_backend}; ignored when
#'   \code{preferred_backend} is non-\code{NULL}.
#' @param preferred_backend Single string naming the preferred compute
#'   backend, e.g. \code{"cpu"}, \code{"mlx"}, or \code{"metal"}.
#' @param policy Single string; one of \code{"auto"}, \code{"cpu"},
#'   \code{"mlx"}, \code{"metal"}, \code{"arrayfire"},
#'   \code{"torch"}.
#' @param precision Single string; \code{"strict"} for full
#'   double-precision accuracy or \code{"fast"} to allow reduced
#'   precision on GPU backends.
#'
#' @return An \code{adgeMatrix} with the data from \code{x} and the
#'   requested backend metadata.
#'
#' @examples
#' m <- matrix(1:6, nrow = 2)
#' A <- adgeMatrix(m)
#' A
#'
#' @export
adgeMatrix <- function(
  x,
  mode             = NULL,
  backend          = NULL,
  preferred_backend = NULL,
  policy           = NULL,
  precision        = NULL
) {
  params <- .amatrix_resolve_mode(mode, backend, preferred_backend, policy, precision)
  new_adgeMatrix(x, preferred_backend = params$preferred_backend,
                 policy = params$policy, precision = params$precision)
}

#' Create a backend-aware sparse matrix
#'
#' Converts a sparse or dense matrix to an \code{adgCMatrix} with the
#' specified backend, policy, and precision. This is the primary
#' user-facing constructor for sparse amatrix objects.
#'
#' @param x A \code{dgCMatrix}, other \code{sparseMatrix}, or base R
#'   \code{matrix}.
#' @param mode Single string shortcut passed to
#'   \code{.amatrix_resolve_mode()}. Pass \code{NULL} to use the
#'   individual arguments.
#' @param backend Alias for \code{preferred_backend}; ignored when
#'   \code{preferred_backend} is non-\code{NULL}.
#' @param preferred_backend Single string naming the preferred compute
#'   backend.
#' @param policy Single string; one of \code{"auto"}, \code{"cpu"},
#'   \code{"mlx"}, \code{"metal"}, \code{"arrayfire"},
#'   \code{"torch"}.
#' @param precision Single string; \code{"strict"} or \code{"fast"}.
#'
#' @return An \code{adgCMatrix} with the data from \code{x} and the
#'   requested backend metadata.
#'
#' @examples
#' m <- matrix(c(1, 0, 0, 2), nrow = 2)
#' S <- adgCMatrix(m)
#' S
#'
#' @export
adgCMatrix <- function(
  x,
  mode             = NULL,
  backend          = NULL,
  preferred_backend = NULL,
  policy           = NULL,
  precision        = NULL
) {
  params <- .amatrix_resolve_mode(mode, backend, preferred_backend, policy, precision)
  new_adgCMatrix(x, preferred_backend = params$preferred_backend,
                 policy = params$policy, precision = params$precision)
}

#' Coerce an object to adgeMatrix
#'
#' Converts a matrix-like object or a \code{resident_handle} to an
#' \code{adgeMatrix}. When \code{x} is a \code{resident_handle},
#' ownership of the GPU-resident buffer is transferred to the new
#' \code{adgeMatrix} with no host-side copy.
#'
#' @param x A \code{resident_handle}, base R \code{matrix},
#'   \code{dgeMatrix}, or any \code{denseMatrix}.
#' @param mode Single string shortcut; see \code{\link{adgeMatrix}}.
#' @param backend Alias for \code{preferred_backend}.
#' @param preferred_backend Single string; preferred compute backend.
#' @param policy Single string dispatch policy.
#' @param precision Single string; \code{"strict"} or \code{"fast"}.
#'
#' @return An \code{adgeMatrix}. When \code{x} is a
#'   \code{resident_handle} the host copy is deferred.
#'
#' @examples
#' m <- matrix(1:6, nrow = 2)
#' A <- as_adgeMatrix(m)
#' dim(A)
#'
#' @export
as_adgeMatrix <- function(
  x,
  mode             = NULL,
  backend          = NULL,
  preferred_backend = NULL,
  policy           = NULL,
  precision        = NULL
) {
  # resident_handle → adgeMatrix with ownership transfer
  if (inherits(x, "resident_handle")) {
    return(as_adgeMatrix.resident_handle(x))
  }
  params <- .amatrix_resolve_mode(mode, backend, preferred_backend, policy, precision)
  new_adgeMatrix(x, preferred_backend = params$preferred_backend,
                 policy = params$policy, precision = params$precision)
}

#' Coerce an object to adgCMatrix
#'
#' Converts a sparse or dense matrix-like object to an
#' \code{adgCMatrix} with the requested backend metadata.
#'
#' @param x A \code{dgCMatrix}, other \code{sparseMatrix}, or base R
#'   \code{matrix}.
#' @param mode Single string shortcut; see \code{\link{adgCMatrix}}.
#' @param backend Alias for \code{preferred_backend}.
#' @param preferred_backend Single string; preferred compute backend.
#' @param policy Single string dispatch policy.
#' @param precision Single string; \code{"strict"} or \code{"fast"}.
#'
#' @return An \code{adgCMatrix}.
#'
#' @export
as_adgCMatrix <- function(
  x,
  mode             = NULL,
  backend          = NULL,
  preferred_backend = NULL,
  policy           = NULL,
  precision        = NULL
) {
  params <- .amatrix_resolve_mode(mode, backend, preferred_backend, policy, precision)
  new_adgCMatrix(x, preferred_backend = params$preferred_backend,
                 policy = params$policy, precision = params$precision)
}

# Internal constructor for zero-copy transpose views.
# Always creates a new object_id distinct from the source so caches and
# the residency registry remain independent.
.new_aTransposeView <- function(source) {
  oid <- .amatrix_next_object_id()
  dn <- source@Dimnames
  new(
    "aTransposeView",
    source           = source,
    Dim              = rev(source@Dim),
    Dimnames         = if (length(dn) == 2L) rev(dn) else list(NULL, NULL),
    preferred_backend = source@preferred_backend,
    policy           = source@policy,
    precision        = source@precision,
    object_id        = oid,
    src_id           = source@object_id,
    finalizer_env    = .amatrix_make_finalizer_env(oid)
  )
}

setAs("matrix", "adgeMatrix", function(from) new_adgeMatrix(from))
setAs("dgeMatrix", "adgeMatrix", function(from) new_adgeMatrix(from))
setAs("matrix", "adgCMatrix", function(from) new_adgCMatrix(from))
setAs("dgCMatrix", "adgCMatrix", function(from) new_adgCMatrix(from))
setAs("adgeMatrix", "dgeMatrix", function(from) new("dgeMatrix", x = from@x, Dim = from@Dim, Dimnames = from@Dimnames, factors = from@factors))
setAs("adgCMatrix", "dgCMatrix", function(from) new("dgCMatrix", i = from@i, p = from@p, Dim = from@Dim, Dimnames = from@Dimnames, x = from@x, factors = from@factors))
