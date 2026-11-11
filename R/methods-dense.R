#' Matrix multiplication for adgeMatrix
#'
#' Dispatches \code{\%*\%} through the amatrix backend for dense
#' \code{adgeMatrix} objects on the left-hand side, preserving GPU
#' residency across the operation.
#'
#' @param x An \code{adgeMatrix}, \code{numeric} vector, or \code{matrix}.
#' @param y A matrix-like object: \code{matrix}, \code{Matrix},
#'   \code{adgeMatrix}, \code{adgCMatrix}, or \code{ANY}.
#'
#' @return An \code{adgeMatrix} (or \code{numeric} vector when \code{y}
#'   is a vector and \code{x} is \code{adgeMatrix}), with the same backend
#'   metadata as \code{x}.
#'
#' @examples
#' A <- adgeMatrix(matrix(1:6, 2, 3))
#' B <- matrix(1:3, 3, 1)
#' A %*% B
#'
#' @name matmul-methods
#' @rdname matmul-methods
#' @aliases %*%,adgeMatrix,ANY-method %*%,adgeMatrix,matrix-method
#'   %*%,adgeMatrix,Matrix-method %*%,adgeMatrix,dgeMatrix-method
#'   %*%,adgeMatrix,dgCMatrix-method %*%,adgeMatrix,adgeMatrix-method
#'   %*%,adgeMatrix,adgCMatrix-method %*%,numeric,adgeMatrix-method
#'   %*%,matrix,adgeMatrix-method %*%,aTransposeView,ANY-method
#'   %*%,aTransposeView,adgeMatrix-method %*%,aTransposeView,matrix-method
#'   %*%,aTransposeView,aTransposeView-method
#'   %*%,adgeMatrix,aTransposeView-method %*%,matrix,aTransposeView-method
#'   %*%,KronMatrix,matrix-method %*%,KronMatrix,numeric-method
#'   %*%,matrix,KronMatrix-method %*%,numeric,KronMatrix-method
#'   %*%,dgCMatrix,adgCMatrix-method %*%,dgeMatrix,adgCMatrix-method
#'   %*%,matrix,adgCMatrix-method %*%,numeric,adgCMatrix-method
NULL
setMethod("%*%", signature(x = "adgeMatrix", y = "ANY"),       function(x, y) matmul(x, y))
#' @rdname matmul-methods
#' @aliases %*%,adgeMatrix,matrix-method
setMethod("%*%", signature(x = "adgeMatrix", y = "matrix"),    function(x, y) matmul(x, y))
#' @rdname matmul-methods
#' @aliases %*%,adgeMatrix,Matrix-method
setMethod("%*%", signature(x = "adgeMatrix", y = "Matrix"),    function(x, y) matmul(x, y))
#' @rdname matmul-methods
#' @aliases %*%,adgeMatrix,dgeMatrix-method
setMethod("%*%", signature(x = "adgeMatrix", y = "dgeMatrix"), function(x, y) matmul(x, y))
#' @rdname matmul-methods
#' @aliases %*%,adgeMatrix,dgCMatrix-method
setMethod("%*%", signature(x = "adgeMatrix", y = "dgCMatrix"), function(x, y) matmul(x, y))
#' @rdname matmul-methods
#' @aliases %*%,adgeMatrix,adgeMatrix-method
setMethod("%*%", signature(x = "adgeMatrix", y = "adgeMatrix"),function(x, y) matmul(x, y))
#' @rdname matmul-methods
#' @aliases %*%,adgeMatrix,adgCMatrix-method
setMethod("%*%", signature(x = "adgeMatrix", y = "adgCMatrix"),function(x, y) matmul(x, y))

# Left-hand numeric/matrix — required for irlba's `mult(v, A)` pattern and any
# row-vector times adgeMatrix call. Without these, S4 falls through to base::%*%
# which coerces A to a plain matrix, silently destroying GPU residency.
#
# Route through y (the adgeMatrix) for GPU dispatch. matmul cannot be used
# here because it calls .amatrix_backend_for(x, ...) which expects x to be an
# aMatrix with @preferred_backend slot.
#
# numeric %*% adgeMatrix follows base R: a vector of length nrow(A) is a row
# vector (result 1 x n); if A has a single row the vector is a column instead
# (outer product). The row-vector case is computed as t(A) %*% x_col through
# the backend crossprod path and reshaped to 1 x n (only n values move).
#' @rdname matmul-methods
#' @aliases %*%,numeric,adgeMatrix-method
setMethod("%*%", signature(x = "numeric", y = "adgeMatrix"), function(x, y) {
  wrap <- function(m) {
    new_adgeMatrix(m,
      preferred_backend = y@preferred_backend,
      policy            = y@policy,
      precision         = y@precision)
  }
  if (length(x) == nrow(y)) {
    col <- am_crossprod(y, matrix(as.double(x), ncol = 1L))
    row <- matrix(as.double(as.matrix(col)), nrow = 1L,
                  dimnames = list(NULL, colnames(y)))
    return(wrap(row))
  }
  if (nrow(y) == 1L) {
    return(matmul(wrap(matrix(as.double(x), ncol = 1L)), y))
  }
  stop("non-conformable arguments", call. = FALSE)
})

#' @rdname matmul-methods
#' @aliases %*%,matrix,adgeMatrix-method
# matrix %*% adgeMatrix: x(k×m) %*% A(m×n) = t(am_crossprod(A, t(x))) = k×n.
#' @noRd
setMethod("%*%", signature(x = "matrix",  y = "adgeMatrix"), function(x, y) {
  # x(k×m) %*% y(m×n): wrap x as adgeMatrix and use the standard matmul path.
  matmul(new_adgeMatrix(x,
    preferred_backend = y@preferred_backend,
    policy            = y@policy,
    precision         = y@precision), y)
})

if (!isGeneric("t")) setGeneric("t", function(x) base::t(x))
#' @noRd
setMethod("t", "adgeMatrix",     function(x) am_transpose(x))
#' @noRd
setMethod("t", "aTransposeView", function(x) x@source)

# --- aTransposeView dispatch ------------------------------------------------

#' @noRd
setMethod("show", "aTransposeView", function(object) {
  cat(sprintf(
    "An amatrix transpose view [%s|policy=%s|precision=%s] %d x %d\n",
    object@preferred_backend, object@policy, object@precision,
    object@Dim[1L], object@Dim[2L]
  ))
})

#' @noRd
setMethod("dim",      "aTransposeView", function(x) x@Dim)
#' @noRd
setMethod("nrow",     "aTransposeView", function(x) x@Dim[1L])
#' @noRd
setMethod("ncol",     "aTransposeView", function(x) x@Dim[2L])
#' @noRd
setMethod("dimnames", "aTransposeView", function(x) x@Dimnames)

.amatrix_subset_transpose_view <- function(x, i, j, ..., drop = TRUE, .linear = FALSE) {
  host <- t(as.matrix(amatrix_materialize_dense(x@source)))
  if (missing(i) && missing(j)) {
    value <- host[, , ..., drop = drop]
  } else if (missing(j)) {
    value <- if (isTRUE(.linear)) host[i] else host[i, , ..., drop = drop]
  } else if (missing(i)) {
    value <- host[, j, ..., drop = drop]
  } else {
    value <- host[i, j, ..., drop = drop]
  }

  if (is.matrix(value)) {
    return(new_adgeMatrix(
      value,
      preferred_backend = x@preferred_backend,
      policy = x@policy,
      precision = x@precision,
      src_id = x@source@object_id
    ))
  }

  value
}

#' @noRd
setMethod("[", signature(x = "aTransposeView", i = "ANY", j = "ANY", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  .amatrix_subset_transpose_view(x, i, j, ..., drop = drop, .linear = missing(j) && (nargs() - !missing(drop)) <= 2L)
})
#' @noRd
setMethod("[", signature(x = "aTransposeView", i = "index", j = "index", drop = "logical"), function(x, i, j, ..., drop = TRUE) {
  .amatrix_subset_transpose_view(x, i, j, ..., drop = drop, .linear = missing(j) && (nargs() - !missing(drop)) <= 2L)
})
#' @noRd
setMethod("[", signature(x = "aTransposeView", i = "index", j = "index", drop = "missing"), function(x, i, j, ..., drop) {
  .amatrix_subset_transpose_view(x, i, j, ..., drop = TRUE, .linear = missing(j) && (nargs() - !missing(drop)) <= 2L)
})
#' @noRd
setMethod("[", signature(x = "aTransposeView", i = "missing", j = "index", drop = "logical"), function(x, i, j, ..., drop = TRUE) {
  .amatrix_subset_transpose_view(x, i, j, ..., drop = drop, .linear = missing(j) && (nargs() - !missing(drop)) <= 2L)
})
#' @noRd
setMethod("[", signature(x = "aTransposeView", i = "missing", j = "index", drop = "missing"), function(x, i, j, ..., drop) {
  .amatrix_subset_transpose_view(x, i, j, ..., drop = TRUE, .linear = missing(j) && (nargs() - !missing(drop)) <= 2L)
})
#' @noRd
setMethod("[", signature(x = "aTransposeView", i = "index", j = "missing", drop = "logical"), function(x, i, j, ..., drop = TRUE) {
  .amatrix_subset_transpose_view(x, i, j, ..., drop = drop, .linear = missing(j) && (nargs() - !missing(drop)) <= 2L)
})
#' @noRd
setMethod("[", signature(x = "aTransposeView", i = "index", j = "missing", drop = "missing"), function(x, i, j, ..., drop) {
  .amatrix_subset_transpose_view(x, i, j, ..., drop = TRUE, .linear = missing(j) && (nargs() - !missing(drop)) <= 2L)
})

# t(A) %*% B — route to crossprod(A, B) using the source resident key
#' @noRd
setMethod("%*%", signature(x = "aTransposeView", y = "adgeMatrix"),
  function(x, y) am_crossprod(x@source, y))
#' @noRd
setMethod("%*%", signature(x = "aTransposeView", y = "matrix"),
  function(x, y) am_crossprod(x@source, y))
#' @noRd
setMethod("%*%", signature(x = "aTransposeView", y = "ANY"),
  function(x, y) am_crossprod(x@source, y))

# A %*% t(B) — route to tcrossprod(A, B) using the source resident key
#' @noRd
setMethod("%*%", signature(x = "adgeMatrix",     y = "aTransposeView"),
  function(x, y) am_tcrossprod(x, y@source))
#' @noRd
setMethod("%*%", signature(x = "matrix",         y = "aTransposeView"),
  function(x, y) am_tcrossprod(new_adgeMatrix(x,
    preferred_backend = y@source@preferred_backend,
    policy            = y@source@policy,
    precision         = y@source@precision), y@source))

# t(A) %*% t(B) = tcrossprod(B, A)
#' @noRd
setMethod("%*%", signature(x = "aTransposeView", y = "aTransposeView"),
  function(x, y) am_tcrossprod(y@source, x@source))

# Arithmetic: materialize to adgeMatrix then delegate
#' @noRd
setMethod("Ops", signature(e1 = "aTransposeView", e2 = "ANY"), function(e1, e2) {
  callGeneric(as(e1, "adgeMatrix"), e2)
})
#' @noRd
setMethod("Ops", signature(e1 = "ANY", e2 = "aTransposeView"), function(e1, e2) {
  callGeneric(e1, as(e2, "adgeMatrix"))
})

setAs("aTransposeView", "adgeMatrix", function(from) {
  new_adgeMatrix(t(as.matrix(amatrix_materialize_dense(from@source))),
    preferred_backend = from@preferred_backend,
    policy            = from@policy,
    precision         = from@precision)
})

# ── aTransposeView as an ordinary matrix ─────────────────────────────────────
#
# t(<adgeMatrix>) returns a lazy view so that t(A) %*% B can be lowered to a
# crossprod without materialising the transpose. The view is not a Matrix, so
# everything else a caller may reasonably do with t(A) is forwarded to the
# concrete transposed adgeMatrix here. Without these methods length(t(A)) was
# 1, mean(t(A)) was NA and sum/rowSums/exp/solve/... raised errors.

.amatrix_view_concrete <- function(x) methods::as(x, "adgeMatrix")

#' @noRd
setMethod("length", "aTransposeView", function(x) prod(x@Dim))
#' @noRd
setMethod("Summary", "aTransposeView", function(x, ..., na.rm = FALSE) {
  callGeneric(.amatrix_view_concrete(x), ..., na.rm = na.rm)
})
#' @noRd
setMethod("Math", "aTransposeView", function(x) callGeneric(.amatrix_view_concrete(x)))
#' @noRd
setMethod("Math2", "aTransposeView", function(x, digits) {
  if (missing(digits)) callGeneric(.amatrix_view_concrete(x)) else callGeneric(.amatrix_view_concrete(x), digits)
})
#' @noRd
setMethod("log", "aTransposeView", function(x, ...) log(.amatrix_view_concrete(x), ...))
#' @noRd
setMethod("Ops", signature(e1 = "aTransposeView", e2 = "missing"), function(e1, e2) {
  callGeneric(.amatrix_view_concrete(e1))
})
#' @noRd
setMethod("as.vector", "aTransposeView", function(x, mode = "any") as.vector(as.matrix(x), mode = mode))
#' @noRd
setMethod("as.numeric", "aTransposeView", function(x, ...) as.numeric(as.matrix(x)))
#' @noRd
setMethod("mean", "aTransposeView", function(x, ...) mean(as.matrix(x), ...))
#' @noRd
setMethod("rowSums", "aTransposeView", function(x, na.rm = FALSE, dims = 1L, ...) {
  rowSums(.amatrix_view_concrete(x), na.rm = na.rm, dims = dims, ...)
})
#' @noRd
setMethod("colSums", "aTransposeView", function(x, na.rm = FALSE, dims = 1L, ...) {
  colSums(.amatrix_view_concrete(x), na.rm = na.rm, dims = dims, ...)
})
#' @noRd
setMethod("rowMeans", "aTransposeView", function(x, na.rm = FALSE, dims = 1L, ...) {
  rowMeans(.amatrix_view_concrete(x), na.rm = na.rm, dims = dims, ...)
})
#' @noRd
setMethod("colMeans", "aTransposeView", function(x, na.rm = FALSE, dims = 1L, ...) {
  colMeans(.amatrix_view_concrete(x), na.rm = na.rm, dims = dims, ...)
})
#' @noRd
setMethod("diag", "aTransposeView", function(x = 1, nrow, ncol, names = TRUE) {
  diag(.amatrix_view_concrete(x), names = names)
})
#' @noRd
setMethod("crossprod", signature(x = "aTransposeView", y = "ANY"), function(x, y = NULL, ...) {
  if (inherits(y, "aTransposeView")) y <- .amatrix_view_concrete(y)
  crossprod(.amatrix_view_concrete(x), y, ...)
})
#' @noRd
setMethod("crossprod", signature(x = "aTransposeView", y = "missing"), function(x, y = NULL, ...) {
  # crossprod(t(A)) is A %*% t(A).
  tcrossprod(x@source)
})
#' @noRd
setMethod("tcrossprod", signature(x = "aTransposeView", y = "ANY"), function(x, y = NULL, ...) {
  if (inherits(y, "aTransposeView")) y <- .amatrix_view_concrete(y)
  tcrossprod(.amatrix_view_concrete(x), y, ...)
})
#' @noRd
setMethod("tcrossprod", signature(x = "aTransposeView", y = "missing"), function(x, y = NULL, ...) {
  # tcrossprod(t(A)) is t(A) %*% A.
  crossprod(x@source)
})
#' @noRd
setMethod("solve", signature(a = "aTransposeView", b = "ANY"), function(a, b, ...) {
  if (missing(b)) solve(.amatrix_view_concrete(a), ...) else solve(.amatrix_view_concrete(a), b, ...)
})
#' @noRd
setMethod("solve", signature(a = "aTransposeView", b = "missing"), function(a, b, ...) {
  solve(.amatrix_view_concrete(a), ...)
})
#' @noRd
setMethod("norm", signature(x = "aTransposeView", type = "missing"), function(x, type, ...) {
  norm(.amatrix_view_concrete(x))
})
#' @noRd
setMethod("norm", signature(x = "aTransposeView", type = "character"), function(x, type, ...) {
  norm(.amatrix_view_concrete(x), type)
})
#' @noRd
setMethod("%*%", signature(x = "numeric", y = "aTransposeView"), function(x, y) {
  x %*% .amatrix_view_concrete(y)
})

# Matrix's isSymmetric() method for dgeMatrix calls t(), which for an
# adgeMatrix yields the lazy view above and then fails on a missing `@x`.
#' @noRd
setMethod("isSymmetric", "adgeMatrix", function(object, ...) {
  host <- amatrix_materialize_dense(object)
  isSymmetric(as.matrix(host), ...)
})

#' Cross-product methods for adgeMatrix
#'
#' Compute \eqn{t(x) \%*\% y} (\code{crossprod}) or
#' \eqn{x \%*\% t(y)} (\code{tcrossprod}) for \code{adgeMatrix} objects,
#' dispatching through the amatrix backend to preserve GPU residency.
#'
#' @param x An \code{adgeMatrix}.
#' @param y A matrix-like object, or \code{NULL} for the symmetric form
#'   \eqn{t(x) \%*\% x} or \eqn{x \%*\% t(x)}.
#' @param ... Further arguments passed to the underlying backend operation.
#'
#' @return An \code{adgeMatrix} containing the result.
#'
#' @examples
#' A <- adgeMatrix(matrix(rnorm(12), 4, 3))
#' crossprod(A)
#' tcrossprod(A)
#'
#' @rdname crossprod-methods
#' @aliases crossprod,adgeMatrix,ANY-method
setMethod("crossprod", signature(x = "adgeMatrix", y = "ANY"), function(x, y = NULL, ...) am_crossprod(x, y = y, ...))
#' @rdname crossprod-methods
#' @aliases crossprod,adgeMatrix,missing-method
setMethod("crossprod", signature(x = "adgeMatrix", y = "missing"), function(x, y, ...) am_crossprod(x, y = NULL, ...))
#' @rdname crossprod-methods
#' @aliases tcrossprod,adgeMatrix,ANY-method
setMethod("tcrossprod", signature(x = "adgeMatrix", y = "ANY"), function(x, y = NULL, ...) am_tcrossprod(x, y = y, ...))
#' @rdname crossprod-methods
#' @aliases tcrossprod,adgeMatrix,missing-method
setMethod("tcrossprod", signature(x = "adgeMatrix", y = "missing"), function(x, y, ...) am_tcrossprod(x, y = NULL, ...))

#' Row and column summary methods for adgeMatrix
#'
#' Compute row or column sums and means for an \code{adgeMatrix}, dispatching
#' through the amatrix backend when GPU acceleration is available.
#'
#' @param x An \code{adgeMatrix}.
#' @param na.rm Logical; if \code{TRUE}, \code{NA} values are ignored.
#' @param dims Integer; dimensions to sum over (passed to the backend).
#'
#' @return A numeric vector of length equal to the number of rows or columns.
#'
#' @examples
#' A <- adgeMatrix(matrix(1:12, 3, 4))
#' rowSums(A)
#' colMeans(A)
#'
#' @rdname rowcol-summary-methods
#' @aliases rowSums,adgeMatrix-method
setMethod("rowSums", "adgeMatrix", function(x, na.rm = FALSE, dims = 1L) rowsums(x, na.rm = na.rm, dims = dims))
#' @rdname rowcol-summary-methods
#' @aliases colSums,adgeMatrix-method
setMethod("colSums", "adgeMatrix", function(x, na.rm = FALSE, dims = 1L) colsums(x, na.rm = na.rm, dims = dims))
#' @rdname rowcol-summary-methods
#' @aliases rowMeans,adgeMatrix-method
setMethod("rowMeans", "adgeMatrix", function(x, na.rm = FALSE, dims = 1L) rowmeans(x, na.rm = na.rm))
#' @rdname rowcol-summary-methods
#' @aliases colMeans,adgeMatrix-method
setMethod("colMeans", "adgeMatrix", function(x, na.rm = FALSE, dims = 1L) colmeans(x, na.rm = na.rm))

#' @noRd
setMethod("cbind2", signature(x = "aMatrix", y = "aMatrix"), function(x, y, ...) {
  .amatrix_bind2("cbind2", x, y)
})

# Explicit dense signatures prevent Matrix's inherited dgeMatrix methods from
# winning before amatrix has a chance to rewrap the result.
#' @noRd
setMethod("cbind2", signature(x = "adgeMatrix", y = "matrix"), function(x, y, ...) {
  .amatrix_bind2("cbind2", x, y)
})

#' @noRd
setMethod("cbind2", signature(x = "adgeMatrix", y = "Matrix"), function(x, y, ...) {
  .amatrix_bind2("cbind2", x, y)
})

#' @noRd
setMethod("cbind2", signature(x = "adgeMatrix", y = "numeric"), function(x, y, ...) {
  .amatrix_bind2("cbind2", x, y)
})

#' @noRd
setMethod("cbind2", signature(x = "adgeMatrix", y = "integer"), function(x, y, ...) {
  .amatrix_bind2("cbind2", x, y)
})

#' @noRd
setMethod("cbind2", signature(x = "adgeMatrix", y = "logical"), function(x, y, ...) {
  .amatrix_bind2("cbind2", x, y)
})

#' @noRd
setMethod("cbind2", signature(x = "aMatrix", y = "ANY"), function(x, y, ...) {
  .amatrix_bind2("cbind2", x, y)
})

#' @noRd
setMethod("cbind2", signature(x = "matrix", y = "adgeMatrix"), function(x, y, ...) {
  .amatrix_bind2("cbind2", x, y)
})

#' @noRd
setMethod("cbind2", signature(x = "Matrix", y = "adgeMatrix"), function(x, y, ...) {
  .amatrix_bind2("cbind2", x, y)
})

#' @noRd
setMethod("cbind2", signature(x = "numeric", y = "adgeMatrix"), function(x, y, ...) {
  .amatrix_bind2("cbind2", x, y)
})

#' @noRd
setMethod("cbind2", signature(x = "integer", y = "adgeMatrix"), function(x, y, ...) {
  .amatrix_bind2("cbind2", x, y)
})

#' @noRd
setMethod("cbind2", signature(x = "logical", y = "adgeMatrix"), function(x, y, ...) {
  .amatrix_bind2("cbind2", x, y)
})

#' @noRd
setMethod("cbind2", signature(x = "ANY", y = "aMatrix"), function(x, y, ...) {
  .amatrix_bind2("cbind2", x, y)
})

#' @noRd
setMethod("cbind2", signature(x = "matrix", y = "aMatrix"), function(x, y, ...) {
  .amatrix_bind2("cbind2", x, y)
})

#' @noRd
setMethod("cbind2", signature(x = "Matrix", y = "aMatrix"), function(x, y, ...) {
  .amatrix_bind2("cbind2", x, y)
})

#' @noRd
setMethod("rbind2", signature(x = "aMatrix", y = "aMatrix"), function(x, y, ...) {
  .amatrix_bind2("rbind2", x, y)
})

#' @noRd
setMethod("rbind2", signature(x = "adgeMatrix", y = "matrix"), function(x, y, ...) {
  .amatrix_bind2("rbind2", x, y)
})

#' @noRd
setMethod("rbind2", signature(x = "adgeMatrix", y = "Matrix"), function(x, y, ...) {
  .amatrix_bind2("rbind2", x, y)
})

#' @noRd
setMethod("rbind2", signature(x = "adgeMatrix", y = "numeric"), function(x, y, ...) {
  .amatrix_bind2("rbind2", x, y)
})

#' @noRd
setMethod("rbind2", signature(x = "adgeMatrix", y = "integer"), function(x, y, ...) {
  .amatrix_bind2("rbind2", x, y)
})

#' @noRd
setMethod("rbind2", signature(x = "adgeMatrix", y = "logical"), function(x, y, ...) {
  .amatrix_bind2("rbind2", x, y)
})

#' @noRd
setMethod("rbind2", signature(x = "aMatrix", y = "ANY"), function(x, y, ...) {
  .amatrix_bind2("rbind2", x, y)
})

#' @noRd
setMethod("rbind2", signature(x = "matrix", y = "adgeMatrix"), function(x, y, ...) {
  .amatrix_bind2("rbind2", x, y)
})

#' @noRd
setMethod("rbind2", signature(x = "Matrix", y = "adgeMatrix"), function(x, y, ...) {
  .amatrix_bind2("rbind2", x, y)
})

#' @noRd
setMethod("rbind2", signature(x = "numeric", y = "adgeMatrix"), function(x, y, ...) {
  .amatrix_bind2("rbind2", x, y)
})

#' @noRd
setMethod("rbind2", signature(x = "integer", y = "adgeMatrix"), function(x, y, ...) {
  .amatrix_bind2("rbind2", x, y)
})

#' @noRd
setMethod("rbind2", signature(x = "logical", y = "adgeMatrix"), function(x, y, ...) {
  .amatrix_bind2("rbind2", x, y)
})

#' @noRd
setMethod("rbind2", signature(x = "ANY", y = "aMatrix"), function(x, y, ...) {
  .amatrix_bind2("rbind2", x, y)
})

#' @noRd
setMethod("rbind2", signature(x = "matrix", y = "aMatrix"), function(x, y, ...) {
  .amatrix_bind2("rbind2", x, y)
})

#' @noRd
setMethod("rbind2", signature(x = "Matrix", y = "aMatrix"), function(x, y, ...) {
  .amatrix_bind2("rbind2", x, y)
})

#' @noRd
setMethod("[", signature(x = "adgeMatrix", i = "ANY", j = "ANY", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  am_subset(x, i, j, ..., drop = drop, .linear = missing(j) && (nargs() - !missing(drop)) <= 2L)
})

#' @noRd
setMethod("[", signature(x = "adgeMatrix", i = "index", j = "index", drop = "logical"), function(x, i, j, ..., drop = TRUE) {
  am_subset(x, i, j, ..., drop = drop, .linear = missing(j) && (nargs() - !missing(drop)) <= 2L)
})

#' @noRd
setMethod("[", signature(x = "adgeMatrix", i = "index", j = "index", drop = "missing"), function(x, i, j, ..., drop) {
  am_subset(x, i, j, ..., drop = TRUE, .linear = missing(j) && (nargs() - !missing(drop)) <= 2L)
})

#' @noRd
setMethod("[", signature(x = "adgeMatrix", i = "missing", j = "index", drop = "logical"), function(x, i, j, ..., drop = TRUE) {
  am_subset(x, i, j, ..., drop = drop, .linear = missing(j) && (nargs() - !missing(drop)) <= 2L)
})

#' @noRd
setMethod("[", signature(x = "adgeMatrix", i = "missing", j = "index", drop = "missing"), function(x, i, j, ..., drop) {
  am_subset(x, i, j, ..., drop = TRUE, .linear = missing(j) && (nargs() - !missing(drop)) <= 2L)
})

#' @noRd
setMethod("[", signature(x = "adgeMatrix", i = "index", j = "missing", drop = "logical"), function(x, i, j, ..., drop = TRUE) {
  am_subset(x, i, j, ..., drop = drop, .linear = missing(j) && (nargs() - !missing(drop)) <= 2L)
})

#' @noRd
setMethod("[", signature(x = "adgeMatrix", i = "index", j = "missing", drop = "missing"), function(x, i, j, ..., drop) {
  am_subset(x, i, j, ..., drop = TRUE, .linear = missing(j) && (nargs() - !missing(drop)) <= 2L)
})

#' Scalar extraction from an adgeMatrix
#'
#' \code{[[} extracts a single element of an \code{\linkS4class{adgeMatrix}},
#' matching base matrix semantics: \code{A[[i]]} selects the \code{i}-th element
#' in column-major order, and \code{A[[i, j]]} selects the element in row
#' \code{i}, column \code{j}. Without this method \code{A[[1]]} fails with
#' \dQuote{this S4 class is not subsettable} (as it does for a bare
#' \code{dgeMatrix}).
#'
#' The host copy is materialized before extraction (transparently downloading
#' from the device for GPU-resident objects), so the returned value always
#' reflects the current contents.
#'
#' @param x An \code{\linkS4class{adgeMatrix}}.
#' @param i,j Element subscripts, following base matrix \code{[[} semantics.
#' @param ... Unused.
#'
#' @return A length-one numeric value.
#'
#' @examples
#' A <- adgeMatrix(matrix(1:6, 2, 3))
#' A[[1]]
#' A[[2, 3]]
#'
#' @name adgeMatrix-extract2
#' @aliases [[,adgeMatrix-method
#' @keywords internal
setMethod("[[", signature(x = "adgeMatrix"), function(x, i, j, ...) {
  m <- as.matrix(amatrix_materialize_host(x))
  if (missing(j)) m[[i, ...]] else m[[i, j, ...]]
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "ANY", j = "ANY", value = "ANY"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value, .linear = missing(j) && nargs() <= 3L)
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "Matrix", j = "missing", value = "ANY"), function(x, i, j, ..., value) {
  am_subassign(x, i, value = value, .linear = TRUE)
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "ndenseMatrix", j = "missing", value = "numeric"), function(x, i, j, ..., value) {
  am_subassign(x, i, value = value, .linear = TRUE)
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "ngeMatrix", j = "missing", value = "numeric"), function(x, i, j, ..., value) {
  am_subassign(x, i, value = value, .linear = TRUE)
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "matrix", j = "missing", value = "ANY"), function(x, i, j, ..., value) {
  am_subassign(x, i, value = value, .linear = TRUE)
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "index", j = "index", value = "numeric"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value, .linear = missing(j) && nargs() <= 3L)
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "index", j = "index", value = "integer"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value, .linear = missing(j) && nargs() <= 3L)
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "index", j = "index", value = "logical"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value, .linear = missing(j) && nargs() <= 3L)
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "index", j = "index", value = "matrix"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value, .linear = missing(j) && nargs() <= 3L)
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "index", j = "index", value = "Matrix"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value, .linear = missing(j) && nargs() <= 3L)
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "missing", j = "index", value = "ANY"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value, .linear = missing(j) && nargs() <= 3L)
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "index", j = "missing", value = "ANY"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value, .linear = missing(j) && nargs() <= 3L)
})

#' Solve a linear system for adgeMatrix
#'
#' Compute the solution to \eqn{a x = b} or the matrix inverse of \code{a}
#' when \code{b} is missing, dispatching through the amatrix backend.
#'
#' @param a An \code{adgeMatrix} coefficient matrix.
#' @param b A matrix or vector right-hand side, or missing for matrix
#'   inversion.
#' @param ... Further arguments passed to the backend.
#'
#' @return An \code{adgeMatrix} (or numeric vector when \code{b} is a
#'   plain vector) containing the solution.
#'
#' @examples
#' A <- adgeMatrix(crossprod(matrix(rnorm(9), 3, 3)) + 3 * diag(3))
#' solve(A)
#'
#' @rdname solve-methods
#' @aliases solve,adgeMatrix,missing-method
setMethod("solve", signature(a = "adgeMatrix", b = "missing"), function(a, b, ...) am_solve(a, ...))
#' @rdname solve-methods
#' @aliases solve,adgeMatrix,ANY-method
setMethod("solve", signature(a = "adgeMatrix", b = "ANY"), function(a, b, ...) am_solve(a, b = b, ...))

#' Determinant for adgeMatrix
#'
#' Computes the determinant of an \code{adgeMatrix} with base matrix semantics.
#' The method materializes the current host value before delegating to
#' \code{\link[base]{determinant}} so \code{base::det()} works even when
#' Matrix is not attached on the search path.
#'
#' @param x An \code{adgeMatrix}.
#' @param logarithm Logical; if \code{TRUE}, return the logarithm of the
#'   modulus.
#' @param ... Further arguments passed to \code{\link[base]{determinant}}.
#'
#' @return A \code{"det"} object with \code{modulus} and \code{sign}.
#'
#' @examples
#' A <- adgeMatrix(crossprod(matrix(rnorm(9), 3, 3)) + 3 * diag(3))
#' det(A)
#'
#' @keywords internal
#' @method determinant adgeMatrix
#' @export
determinant.adgeMatrix <- function(x, logarithm = TRUE, ...) {
  base::determinant(as.matrix(amatrix_materialize_host(x)),
    logarithm = logarithm, ...)
}

#' @noRd
setMethod("determinant", signature(x = "adgeMatrix", logarithm = "logical"),
  determinant.adgeMatrix)

#' @noRd
setMethod("determinant", signature(x = "adgeMatrix", logarithm = "missing"),
  function(x, logarithm, ...) determinant.adgeMatrix(x, logarithm = TRUE, ...))

#' Cholesky factorization for adgeMatrix
#'
#' Compute the Cholesky factor of a symmetric positive-definite
#' \code{adgeMatrix}, dispatching through the amatrix backend.
#'
#' @param x A symmetric positive-definite \code{adgeMatrix}.
#' @param ... Further arguments passed to the backend.
#'
#' @return An \code{adgeMatrix} containing the upper triangular Cholesky
#'   factor.
#'
#' @examples
#' S <- adgeMatrix(crossprod(matrix(rnorm(9), 3, 3)) + 3 * diag(3))
#' R <- chol(S)
#'
#' @rdname chol-methods
#' @aliases chol,adgeMatrix-method
setMethod("chol", "adgeMatrix", function(x, ...) am_chol(x, ...))
#' @noRd
setMethod("qr", "adgeMatrix", function(x, ...) am_qr(x, ...))

#' Singular value decomposition for adgeMatrix
#'
#' Compute the singular value decomposition of an \code{adgeMatrix},
#' dispatching through the amatrix backend. A fallback method for plain
#' \code{matrix} objects is also provided to preserve base R behaviour
#' after the generic is promoted to S4.
#'
#' @param x An \code{adgeMatrix} or plain \code{matrix}.
#' @param nu Number of left singular vectors to compute.
#' @param nv Number of right singular vectors to compute.
#' @param LINPACK Ignored; retained for signature compatibility.
#' @param ... Further arguments passed to the backend.
#'
#' @return A list with components \code{d} (singular values), \code{u}
#'   (left singular vectors, \code{nrow(x)} by \code{nu}), and \code{v}
#'   (right singular vectors, \code{ncol(x)} by \code{nv}).
#'
#' @examples
#' A <- adgeMatrix(matrix(rnorm(12), 4, 3))
#' s <- svd(A)
#' length(s$d)
#'
#' @rdname svd-methods
#' @name svd-methods
#' @aliases svd svd,adgeMatrix-method
#' @export
setMethod("svd", "adgeMatrix", function(x, nu = min(dim(x)), nv = min(dim(x)), LINPACK = FALSE, ...) {
  am_svd(x, nu = nu, nv = nv, LINPACK = LINPACK, ...)
})

# Plain matrices (and everything else) fall through to base::svd, which is the
# default method of the standard generic declared in R/generics.R.

#' Eigendecomposition for adgeMatrix
#'
#' Compute eigenvalues and eigenvectors of an \code{adgeMatrix},
#' dispatching through the amatrix backend for symmetric matrices when
#' GPU acceleration is available. A fallback method for plain \code{matrix}
#' preserves base R behaviour after the generic is promoted to S4.
#'
#' @param x An \code{adgeMatrix} or plain \code{matrix}.
#' @param symmetric Logical indicating whether \code{x} is symmetric.
#'   When missing, symmetry is detected automatically from the host copy.
#' @param only.values Logical; if \code{TRUE} only eigenvalues are returned.
#' @param EISPACK Ignored; retained for signature compatibility.
#'
#' @return A list with components \code{values} (numeric vector) and
#'   \code{vectors} (matrix, omitted when \code{only.values = TRUE}).
#'
#' @examples
#' S <- adgeMatrix(crossprod(matrix(rnorm(9), 3, 3)))
#' ev <- eigen(S, symmetric = TRUE)
#' length(ev$values)
#'
#' @rdname eigen-methods
#' @aliases eigen,adgeMatrix-method
setMethod("eigen", "adgeMatrix", function(x, symmetric, only.values = FALSE, EISPACK = FALSE) {
  sym <- if (missing(symmetric)) NULL else symmetric
  am_eigen(x, symmetric = sym, only.values = only.values, EISPACK = EISPACK)
})

#' @rdname eigen-methods
#' @aliases eigen,matrix-method
# Fallback: keep base::eigen working for plain matrix
#' @noRd
setMethod("eigen", "matrix", function(x, symmetric, only.values = FALSE, EISPACK = FALSE) {
  sym <- if (missing(symmetric)) isSymmetric(x) else symmetric
  base::eigen(x, symmetric = sym, only.values = only.values, EISPACK = EISPACK)
})

#' @noRd
setMethod("diag", "adgeMatrix", function(x = 1, nrow, ncol, names = TRUE) {
  if (missing(nrow) && missing(ncol)) am_diag(x, names = names)
  else if (missing(ncol))             am_diag(x, nrow = nrow, names = names)
  else                                am_diag(x, nrow = nrow, ncol = ncol, names = names)
})

#' @noRd
setMethod("Math", "adgeMatrix", function(x) {
  .amatrix_rewrap_value(x, callGeneric(amatrix_materialize_host(x)))
})

# log() is in the Math group, but the group method only receives `x`, so
# log(x, base) silently computed the natural log. An explicit method keeps
# the base argument.
#' @noRd
setMethod("log", "adgeMatrix", function(x, ...) {
  .amatrix_rewrap_value(x, log(amatrix_materialize_host(x), ...))
})

#' @noRd
setReplaceMethod("diag", "adgeMatrix", function(x, value) {
  host_x <- amatrix_materialize_host(x)
  base::diag(host_x) <- value
  .amatrix_cache_invalidate_object(x@object_id)
  .amatrix_rewrap_value(x, host_x)
})

#' @noRd
setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "ANY"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

# Explicit same-class and mixed-amatrix signatures prevent S4 from preferring
# Matrix's Ops(dgeMatrix, dgeMatrix) method (adgeMatrix extends dgeMatrix, so
# the inherited method has distance 1+1=2 which beats (adgeMatrix, ANY) when
# ANY has distance > 2 for the second slot).
#' @noRd
setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "adgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "adgCMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "adgCMatrix", e2 = "adgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "numeric"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "integer"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "logical"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "matrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "Matrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "dgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "ANY", e2 = "adgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "numeric", e2 = "adgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "integer", e2 = "adgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "logical", e2 = "adgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "matrix", e2 = "adgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "Matrix", e2 = "adgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "dgeMatrix", e2 = "adgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "adlgeMatrix", e2 = "ANY"), function(e1, e2) {
  .amatrix_rewrap_value(e1, callGeneric(.amatrix_logical_host_arg(e1), .amatrix_logical_host_arg(e2)))
})

#' @noRd
setMethod("Ops", signature(e1 = "adlgeMatrix", e2 = "missing"), function(e1, e2) {
  .amatrix_rewrap_value(e1, callGeneric(.amatrix_logical_host_arg(e1)))
})

#' @noRd
setMethod("Ops", signature(e1 = "adlgeMatrix", e2 = "adlgeMatrix"), function(e1, e2) {
  .amatrix_rewrap_value(e1, callGeneric(.amatrix_logical_host_arg(e1), .amatrix_logical_host_arg(e2)))
})

#' @noRd
setMethod("Ops", signature(e1 = "ANY", e2 = "adlgeMatrix"), function(e1, e2) {
  .amatrix_rewrap_value(e2, callGeneric(.amatrix_logical_host_arg(e1), .amatrix_logical_host_arg(e2)))
})

#' @noRd
setMethod("!", "adlgeMatrix", function(x) {
  .amatrix_rewrap_value(x, !.amatrix_logical_host_arg(x))
})

# norm() methods live in R/wrappers.R — removed duplicate definitions here
