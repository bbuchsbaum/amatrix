#' Matrix multiplication for adgCMatrix
#'
#' Dispatches \code{\%*\%} through the amatrix backend for sparse
#' \code{adgCMatrix} objects on the left-hand side, preserving GPU
#' residency metadata across the operation.
#'
#' @param x An \code{adgCMatrix}.
#' @param y A matrix-like object: \code{matrix}, \code{Matrix},
#'   \code{adgeMatrix}, \code{adgCMatrix}, or \code{ANY}.
#'
#' @return An \code{adgeMatrix} or \code{adgCMatrix} containing the
#'   product, with backend metadata inherited from \code{x}.
#'
#' @examples
#' sp <- as(matrix(c(1, 0, 0, 2), 2, 2), "dgCMatrix")
#' A  <- adgCMatrix(sp)
#' B  <- matrix(1:4, 2, 2)
#' A %*% B
#'
#' @rdname matmul-sparse-methods
#' @aliases %*%,adgCMatrix,ANY-method
setMethod("%*%", signature(x = "adgCMatrix", y = "ANY"), function(x, y) matmul(x, y))
#' @rdname matmul-sparse-methods
#' @aliases %*%,adgCMatrix,matrix-method
setMethod("%*%", signature(x = "adgCMatrix", y = "matrix"), function(x, y) matmul(x, y))
#' @rdname matmul-sparse-methods
#' @aliases %*%,adgCMatrix,Matrix-method
setMethod("%*%", signature(x = "adgCMatrix", y = "Matrix"), function(x, y) matmul(x, y))
#' @rdname matmul-sparse-methods
#' @aliases %*%,adgCMatrix,dgeMatrix-method
setMethod("%*%", signature(x = "adgCMatrix", y = "dgeMatrix"), function(x, y) matmul(x, y))
#' @rdname matmul-sparse-methods
#' @aliases %*%,adgCMatrix,dgCMatrix-method
setMethod("%*%", signature(x = "adgCMatrix", y = "dgCMatrix"), function(x, y) matmul(x, y))
#' @rdname matmul-sparse-methods
#' @aliases %*%,adgCMatrix,adgeMatrix-method
setMethod("%*%", signature(x = "adgCMatrix", y = "adgeMatrix"), function(x, y) matmul(x, y))
#' @rdname matmul-sparse-methods
#' @aliases %*%,adgCMatrix,adgCMatrix-method
setMethod("%*%", signature(x = "adgCMatrix", y = "adgCMatrix"), function(x, y) matmul(x, y))

#' @noRd
setMethod("t", "adgCMatrix", function(x) am_transpose(x))

#' Cross-product methods for adgCMatrix
#'
#' Compute \eqn{t(x) \%*\% y} (\code{crossprod}) or
#' \eqn{x \%*\% t(y)} (\code{tcrossprod}) for sparse \code{adgCMatrix}
#' objects, dispatching through the amatrix backend.
#'
#' @param x An \code{adgCMatrix}.
#' @param y A matrix-like object or \code{NULL}/missing for the symmetric
#'   form.
#' @param ... Further arguments passed to the backend.
#'
#' @return An \code{adgeMatrix} or \code{adgCMatrix} containing the result.
#'
#' @examples
#' sp <- as(matrix(c(1, 0, 0, 2, 1, 0), 3, 2), "dgCMatrix")
#' A  <- adgCMatrix(sp)
#' crossprod(A)
#'
#' @rdname crossprod-sparse-methods
#' @aliases crossprod,adgCMatrix,missing-method
setMethod("crossprod",  signature(x = "adgCMatrix", y = "missing"), function(x, y, ...) am_crossprod(x,  y = NULL, ...))
#' @rdname crossprod-sparse-methods
#' @aliases crossprod,adgCMatrix,ANY-method
setMethod("crossprod",  signature(x = "adgCMatrix", y = "ANY"),     function(x, y = NULL, ...) am_crossprod(x,  y = y, ...))
#' @rdname crossprod-sparse-methods
#' @aliases crossprod,adgCMatrix,matrix-method
setMethod("crossprod",  signature(x = "adgCMatrix", y = "matrix"),  function(x, y = NULL, ...) am_crossprod(x,  y = y, ...))
#' @rdname crossprod-sparse-methods
#' @aliases crossprod,adgCMatrix,Matrix-method
setMethod("crossprod",  signature(x = "adgCMatrix", y = "Matrix"),  function(x, y = NULL, ...) am_crossprod(x,  y = y, ...))
#' @rdname crossprod-sparse-methods
#' @aliases crossprod,adgCMatrix,dgeMatrix-method
setMethod("crossprod",  signature(x = "adgCMatrix", y = "dgeMatrix"), function(x, y = NULL, ...) am_crossprod(x, y = y, ...))
#' @rdname crossprod-sparse-methods
#' @aliases crossprod,adgCMatrix,dgCMatrix-method
setMethod("crossprod",  signature(x = "adgCMatrix", y = "dgCMatrix"), function(x, y = NULL, ...) am_crossprod(x, y = y, ...))
#' @rdname crossprod-sparse-methods
#' @aliases crossprod,adgCMatrix,adgeMatrix-method
setMethod("crossprod",  signature(x = "adgCMatrix", y = "adgeMatrix"), function(x, y = NULL, ...) am_crossprod(x, y = y, ...))
#' @rdname crossprod-sparse-methods
#' @aliases crossprod,adgCMatrix,adgCMatrix-method
setMethod("crossprod",  signature(x = "adgCMatrix", y = "adgCMatrix"), function(x, y = NULL, ...) am_crossprod(x, y = y, ...))
#' @rdname crossprod-sparse-methods
#' @aliases tcrossprod,adgCMatrix,missing-method
setMethod("tcrossprod", signature(x = "adgCMatrix", y = "missing"), function(x, y, ...) am_tcrossprod(x, y = NULL, ...))
#' @rdname crossprod-sparse-methods
#' @aliases tcrossprod,adgCMatrix,ANY-method
setMethod("tcrossprod", signature(x = "adgCMatrix", y = "ANY"),     function(x, y = NULL, ...) am_tcrossprod(x, y = y, ...))
#' @rdname crossprod-sparse-methods
#' @aliases tcrossprod,adgCMatrix,matrix-method
setMethod("tcrossprod", signature(x = "adgCMatrix", y = "matrix"),  function(x, y = NULL, ...) am_tcrossprod(x, y = y, ...))
#' @rdname crossprod-sparse-methods
#' @aliases tcrossprod,adgCMatrix,Matrix-method
setMethod("tcrossprod", signature(x = "adgCMatrix", y = "Matrix"),  function(x, y = NULL, ...) am_tcrossprod(x, y = y, ...))
#' @rdname crossprod-sparse-methods
#' @aliases tcrossprod,adgCMatrix,dgeMatrix-method
setMethod("tcrossprod", signature(x = "adgCMatrix", y = "dgeMatrix"), function(x, y = NULL, ...) am_tcrossprod(x, y = y, ...))
#' @rdname crossprod-sparse-methods
#' @aliases tcrossprod,adgCMatrix,dgCMatrix-method
setMethod("tcrossprod", signature(x = "adgCMatrix", y = "dgCMatrix"), function(x, y = NULL, ...) am_tcrossprod(x, y = y, ...))
#' @rdname crossprod-sparse-methods
#' @aliases tcrossprod,adgCMatrix,adgeMatrix-method
setMethod("tcrossprod", signature(x = "adgCMatrix", y = "adgeMatrix"), function(x, y = NULL, ...) am_tcrossprod(x, y = y, ...))
#' @rdname crossprod-sparse-methods
#' @aliases tcrossprod,adgCMatrix,adgCMatrix-method
setMethod("tcrossprod", signature(x = "adgCMatrix", y = "adgCMatrix"), function(x, y = NULL, ...) am_tcrossprod(x, y = y, ...))

#' Row and column summary methods for adgCMatrix
#'
#' Compute row or column sums and means for a sparse \code{adgCMatrix},
#' dispatching through the amatrix backend when GPU acceleration is
#' available.
#'
#' @param x An \code{adgCMatrix}.
#' @param na.rm Logical; if \code{TRUE}, \code{NA} values are ignored.
#' @param dims Integer; dimensions to sum over (passed to the backend).
#'
#' @return A numeric vector of length equal to the number of rows or
#'   columns.
#'
#' @examples
#' sp <- as(matrix(c(1, 0, 2, 0, 3, 0), 2, 3), "dgCMatrix")
#' A  <- adgCMatrix(sp)
#' rowSums(A)
#' colMeans(A)
#'
#' @rdname rowcol-summary-sparse-methods
#' @aliases rowSums,adgCMatrix-method
setMethod("rowSums", "adgCMatrix", function(x, na.rm = FALSE, dims = 1L) rowsums(x, na.rm = na.rm, dims = dims))
#' @rdname rowcol-summary-sparse-methods
#' @aliases colSums,adgCMatrix-method
setMethod("colSums", "adgCMatrix", function(x, na.rm = FALSE, dims = 1L) colsums(x, na.rm = na.rm, dims = dims))
#' @rdname rowcol-summary-sparse-methods
#' @aliases rowMeans,adgCMatrix-method
setMethod("rowMeans", "adgCMatrix", function(x, na.rm = FALSE, dims = 1L) rowmeans(x, na.rm = na.rm))
#' @rdname rowcol-summary-sparse-methods
#' @aliases colMeans,adgCMatrix-method
setMethod("colMeans", "adgCMatrix", function(x, na.rm = FALSE, dims = 1L) colmeans(x, na.rm = na.rm))

#' @noRd
setMethod("[", signature(x = "adgCMatrix", i = "ANY", j = "ANY", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  am_subset(x, i, j, ..., drop = drop)
})

#' @noRd
setMethod("[", signature(x = "adgCMatrix", i = "index", j = "index", drop = "logical"), function(x, i, j, ..., drop = TRUE) {
  am_subset(x, i, j, ..., drop = drop)
})

#' @noRd
setMethod("[", signature(x = "adgCMatrix", i = "missing", j = "index", drop = "logical"), function(x, i, j, ..., drop = TRUE) {
  am_subset(x, i, j, ..., drop = drop)
})

#' @noRd
setMethod("[", signature(x = "adgCMatrix", i = "index", j = "missing", drop = "logical"), function(x, i, j, ..., drop = TRUE) {
  am_subset(x, i, j, ..., drop = drop)
})

setReplaceMethod("[", signature(x = "adgCMatrix", i = "ANY", j = "ANY", value = "ANY"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value)
})

setReplaceMethod("[", signature(x = "adgCMatrix", i = "index", j = "index", value = "numeric"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value)
})

setReplaceMethod("[", signature(x = "adgCMatrix", i = "index", j = "index", value = "integer"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value)
})

setReplaceMethod("[", signature(x = "adgCMatrix", i = "index", j = "index", value = "logical"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value)
})

setReplaceMethod("[", signature(x = "adgCMatrix", i = "index", j = "index", value = "matrix"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value)
})

setReplaceMethod("[", signature(x = "adgCMatrix", i = "index", j = "index", value = "Matrix"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value)
})

#' @noRd
setMethod("norm", "adgCMatrix", function(x, type = "1", ...) {
  Matrix::norm(amatrix_materialize_host(x), type = type)
})

#' Solve a linear system for adgCMatrix
#'
#' Compute the solution to \eqn{a x = b} or the inverse of \code{a} when
#' \code{b} is missing, for a sparse \code{adgCMatrix} coefficient matrix.
#'
#' @param a An \code{adgCMatrix} coefficient matrix.
#' @param b A matrix or vector right-hand side, or missing for inversion.
#' @param ... Further arguments passed to the backend.
#'
#' @return An \code{adgeMatrix} or \code{adgCMatrix} containing the
#'   solution.
#'
#' @rdname solve-sparse-methods
#' @aliases solve,adgCMatrix,missing-method
setMethod("solve", signature(a = "adgCMatrix", b = "missing"), function(a, b, ...) am_solve(a, ...))
#' @rdname solve-sparse-methods
#' @aliases solve,adgCMatrix,ANY-method
setMethod("solve", signature(a = "adgCMatrix", b = "ANY"), function(a, b, ...) am_solve(a, b = b, ...))

#' Cholesky factorization for adgCMatrix
#'
#' Compute the Cholesky factorization of a sparse symmetric
#' positive-definite \code{adgCMatrix}.
#'
#' @param x A symmetric positive-definite \code{adgCMatrix}.
#' @param ... Further arguments passed to \code{Matrix::chol}.
#'
#' @return An \code{adgCMatrix} or sparse Cholesky factor object.
#'
#' @rdname chol-sparse-methods
#' @aliases chol,adgCMatrix-method
setMethod("chol", "adgCMatrix", function(x, ...) am_chol(x, ...))
#' @noRd
setMethod("qr", "adgCMatrix", function(x, ...) am_qr(x, ...))

#' Singular value decomposition for adgCMatrix
#'
#' Compute the singular value decomposition of a sparse
#' \code{adgCMatrix}, dispatching through the amatrix backend.
#'
#' @param x An \code{adgCMatrix}.
#' @param nu Number of left singular vectors to compute.
#' @param nv Number of right singular vectors to compute.
#' @param LINPACK Ignored; retained for signature compatibility.
#' @param ... Further arguments passed to the backend.
#'
#' @return A list with components \code{d}, \code{u}, and \code{v}.
#'
#' @rdname svd-sparse-methods
#' @aliases svd,adgCMatrix-method
setMethod("svd", "adgCMatrix", function(x, nu = min(dim(x)), nv = min(dim(x)), LINPACK = FALSE, ...) {
  am_svd(x, nu = nu, nv = nv, LINPACK = LINPACK, ...)
})

#' Eigendecomposition for adgCMatrix
#'
#' Compute eigenvalues and eigenvectors of a sparse \code{adgCMatrix},
#' dispatching through the amatrix backend.
#'
#' @param x An \code{adgCMatrix}.
#' @param symmetric Logical; whether \code{x} is symmetric. Auto-detected
#'   when missing.
#' @param only.values Logical; if \code{TRUE} only eigenvalues are returned.
#' @param EISPACK Ignored; retained for signature compatibility.
#'
#' @return A list with components \code{values} and \code{vectors}.
#'
#' @rdname eigen-sparse-methods
#' @aliases eigen,adgCMatrix-method
setMethod("eigen", "adgCMatrix", function(x, symmetric, only.values = FALSE, EISPACK = FALSE) {
  sym <- if (missing(symmetric)) NULL else symmetric
  am_eigen(x, symmetric = sym, only.values = only.values, EISPACK = EISPACK)
})

#' @noRd
setMethod("diag", "adgCMatrix", function(x = 1, nrow, ncol, names = TRUE) {
  if (missing(nrow) && missing(ncol)) am_diag(x, names = names)
  else if (missing(ncol))             am_diag(x, nrow = nrow, names = names)
  else                                am_diag(x, nrow = nrow, ncol = ncol, names = names)
})

#' @noRd
setMethod("Ops", signature(e1 = "adgCMatrix", e2 = "ANY"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

# Explicit same-class signature: adgCMatrix extends dgCMatrix, so without this
# Matrix's Ops(dgCMatrix, dgCMatrix) would win by distance 1+1 vs (adgCMatrix, ANY).
#' @noRd
setMethod("Ops", signature(e1 = "adgCMatrix", e2 = "adgCMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "adgCMatrix", e2 = "numeric"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "adgCMatrix", e2 = "integer"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "adgCMatrix", e2 = "logical"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "adgCMatrix", e2 = "matrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "adgCMatrix", e2 = "Matrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "ANY", e2 = "adgCMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "numeric", e2 = "adgCMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "integer", e2 = "adgCMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "logical", e2 = "adgCMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "matrix", e2 = "adgCMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

#' @noRd
setMethod("Ops", signature(e1 = "Matrix", e2 = "adgCMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})
