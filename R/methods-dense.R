setMethod("%*%", signature(x = "adgeMatrix", y = "ANY"),       function(x, y) matmul(x, y))
setMethod("%*%", signature(x = "adgeMatrix", y = "matrix"),    function(x, y) matmul(x, y))
setMethod("%*%", signature(x = "adgeMatrix", y = "Matrix"),    function(x, y) matmul(x, y))
setMethod("%*%", signature(x = "adgeMatrix", y = "dgeMatrix"), function(x, y) matmul(x, y))
setMethod("%*%", signature(x = "adgeMatrix", y = "dgCMatrix"), function(x, y) matmul(x, y))
setMethod("%*%", signature(x = "adgeMatrix", y = "adgeMatrix"),function(x, y) matmul(x, y))
setMethod("%*%", signature(x = "adgeMatrix", y = "adgCMatrix"),function(x, y) matmul(x, y))

# Left-hand numeric/matrix — required for irlba's `mult(v, A)` pattern and any
# row-vector times adgeMatrix call. Without these, S4 falls through to base::%*%
# which coerces A to a plain matrix, silently destroying GPU residency.
#
# Route through y (the adgeMatrix) for GPU dispatch. matmul cannot be used
# here because it calls .amatrix_backend_for(x, ...) which expects x to be an
# aMatrix with @preferred_backend slot.
#
# numeric %*% adgeMatrix: am_crossprod(A, x_col) = t(A) %*% x_col = n×1.
# Values identical to x_row %*% A (1×n); drop()/t() in irlba normalise both.
setMethod("%*%", signature(x = "numeric", y = "adgeMatrix"), function(x, y) {
  am_crossprod(y, matrix(x, ncol = 1L))
})

# matrix %*% adgeMatrix: x(k×m) %*% A(m×n) = t(am_crossprod(A, t(x))) = k×n.
setMethod("%*%", signature(x = "matrix",  y = "adgeMatrix"), function(x, y) {
  # x(k×m) %*% y(m×n): wrap x as adgeMatrix and use the standard matmul path.
  matmul(new_adgeMatrix(x,
    preferred_backend = y@preferred_backend,
    policy            = y@policy,
    precision         = y@precision), y)
})

if (!isGeneric("t")) setGeneric("t", function(x) base::t(x))
setMethod("t", "adgeMatrix",     function(x) am_transpose(x))
setMethod("t", "aTransposeView", function(x) x@source)

# --- aTransposeView dispatch ------------------------------------------------

setMethod("show", "aTransposeView", function(object) {
  cat(sprintf(
    "An amatrix transpose view [%s|policy=%s|precision=%s] %d x %d\n",
    object@preferred_backend, object@policy, object@precision,
    object@Dim[1L], object@Dim[2L]
  ))
})

setMethod("dim",      "aTransposeView", function(x) x@Dim)
setMethod("nrow",     "aTransposeView", function(x) x@Dim[1L])
setMethod("ncol",     "aTransposeView", function(x) x@Dim[2L])
setMethod("dimnames", "aTransposeView", function(x) x@Dimnames)

# t(A) %*% B — route to crossprod(A, B) using the source resident key
setMethod("%*%", signature(x = "aTransposeView", y = "adgeMatrix"),
  function(x, y) am_crossprod(x@source, y))
setMethod("%*%", signature(x = "aTransposeView", y = "matrix"),
  function(x, y) am_crossprod(x@source, y))
setMethod("%*%", signature(x = "aTransposeView", y = "ANY"),
  function(x, y) am_crossprod(x@source, y))

# A %*% t(B) — route to tcrossprod(A, B) using the source resident key
setMethod("%*%", signature(x = "adgeMatrix",     y = "aTransposeView"),
  function(x, y) am_tcrossprod(x, y@source))
setMethod("%*%", signature(x = "matrix",         y = "aTransposeView"),
  function(x, y) am_tcrossprod(new_adgeMatrix(x,
    preferred_backend = y@source@preferred_backend,
    policy            = y@source@policy,
    precision         = y@source@precision), y@source))

# t(A) %*% t(B) = tcrossprod(B, A)
setMethod("%*%", signature(x = "aTransposeView", y = "aTransposeView"),
  function(x, y) am_tcrossprod(y@source, x@source))

# Arithmetic: materialize to adgeMatrix then delegate
setMethod("Ops", signature(e1 = "aTransposeView", e2 = "ANY"), function(e1, e2) {
  callGeneric(as(e1, "adgeMatrix"), e2)
})
setMethod("Ops", signature(e1 = "ANY", e2 = "aTransposeView"), function(e1, e2) {
  callGeneric(e1, as(e2, "adgeMatrix"))
})

setAs("aTransposeView", "adgeMatrix", function(from) {
  new_adgeMatrix(t(as.matrix(amatrix_materialize_dense(from@source))),
    preferred_backend = from@preferred_backend,
    policy            = from@policy,
    precision         = from@precision)
})

setMethod("crossprod", signature(x = "adgeMatrix", y = "ANY"), function(x, y = NULL, ...) am_crossprod(x, y = y, ...))
setMethod("crossprod", signature(x = "adgeMatrix", y = "missing"), function(x, y, ...) am_crossprod(x, y = NULL, ...))
setMethod("tcrossprod", signature(x = "adgeMatrix", y = "ANY"), function(x, y = NULL, ...) am_tcrossprod(x, y = y, ...))
setMethod("tcrossprod", signature(x = "adgeMatrix", y = "missing"), function(x, y, ...) am_tcrossprod(x, y = NULL, ...))

setMethod("rowSums", "adgeMatrix", function(x, na.rm = FALSE, dims = 1L) rowsums(x, na.rm = na.rm, dims = dims))
setMethod("colSums", "adgeMatrix", function(x, na.rm = FALSE, dims = 1L) colsums(x, na.rm = na.rm, dims = dims))
setMethod("rowMeans", "adgeMatrix", function(x, na.rm = FALSE, dims = 1L) rowmeans(x, na.rm = na.rm))
setMethod("colMeans", "adgeMatrix", function(x, na.rm = FALSE, dims = 1L) colmeans(x, na.rm = na.rm))

setMethod("[", signature(x = "adgeMatrix", i = "ANY", j = "ANY", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  am_subset(x, i, j, ..., drop = drop)
})

setMethod("[", signature(x = "adgeMatrix", i = "index", j = "index", drop = "logical"), function(x, i, j, ..., drop = TRUE) {
  am_subset(x, i, j, ..., drop = drop)
})

setMethod("[", signature(x = "adgeMatrix", i = "missing", j = "index", drop = "logical"), function(x, i, j, ..., drop = TRUE) {
  am_subset(x, i, j, ..., drop = drop)
})

setMethod("[", signature(x = "adgeMatrix", i = "index", j = "missing", drop = "logical"), function(x, i, j, ..., drop = TRUE) {
  am_subset(x, i, j, ..., drop = drop)
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "ANY", j = "ANY", value = "ANY"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value)
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "index", j = "index", value = "numeric"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value)
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "index", j = "index", value = "integer"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value)
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "index", j = "index", value = "logical"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value)
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "index", j = "index", value = "matrix"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value)
})

setReplaceMethod("[", signature(x = "adgeMatrix", i = "index", j = "index", value = "Matrix"), function(x, i, j, ..., value) {
  am_subassign(x, i, j, ..., value = value)
})

setMethod("solve", signature(a = "adgeMatrix", b = "missing"), function(a, b, ...) am_solve(a, ...))
setMethod("solve", signature(a = "adgeMatrix", b = "ANY"), function(a, b, ...) am_solve(a, b = b, ...))
setMethod("chol", "adgeMatrix", function(x, ...) am_chol(x, ...))
setMethod("qr", "adgeMatrix", function(x, ...) am_qr(x, ...))

setMethod("svd", "adgeMatrix", function(x, nu = min(dim(x)), nv = min(dim(x)), LINPACK = FALSE, ...) {
  am_svd(x, nu = nu, nv = nv, LINPACK = LINPACK, ...)
})

setMethod("eigen", "adgeMatrix", function(x, symmetric, only.values = FALSE, EISPACK = FALSE) {
  sym <- if (missing(symmetric)) NULL else symmetric
  am_eigen(x, symmetric = sym, only.values = only.values, EISPACK = EISPACK)
})

setMethod("diag", "adgeMatrix", function(x = 1, nrow, ncol, names = TRUE) {
  if (missing(nrow) && missing(ncol)) am_diag(x, names = names)
  else if (missing(ncol))             am_diag(x, nrow = nrow, names = names)
  else                                am_diag(x, nrow = nrow, ncol = ncol, names = names)
})

setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "ANY"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

# Explicit same-class and mixed-amatrix signatures prevent S4 from preferring
# Matrix's Ops(dgeMatrix, dgeMatrix) method (adgeMatrix extends dgeMatrix, so
# the inherited method has distance 1+1=2 which beats (adgeMatrix, ANY) when
# ANY has distance > 2 for the second slot).
setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "adgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "adgCMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "adgCMatrix", e2 = "adgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "numeric"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "integer"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "logical"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "matrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "Matrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "ANY", e2 = "adgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "numeric", e2 = "adgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "integer", e2 = "adgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "logical", e2 = "adgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "matrix", e2 = "adgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "Matrix", e2 = "adgeMatrix"), function(e1, e2) {
  ewise(.Generic, e1, e2)
})

# ---------------------------------------------------------------------------
# norm()  — matrix / vector norms dispatched on adgeMatrix
# ---------------------------------------------------------------------------
if (!isGeneric("norm")) {
  setGeneric("norm", function(x, type = "F", ...) standardGeneric("norm"))
}

setMethod("norm", "adgeMatrix", function(x, type = "F", ...) {
  type <- match.arg(toupper(substr(type, 1L, 1L)),
                    c("F", "1", "O", "I", "M", "2"))
  X_mat <- as.matrix(amatrix_materialize_host(x))
  switch(type,
    "F" = sqrt(sum(X_mat * X_mat)),
    "1" = ,
    "O" = max(colSums(abs(X_mat))),
    "I" = max(rowSums(abs(X_mat))),
    "M" = max(abs(X_mat)),
    "2" = {
      sv <- rsvd(x, k = 1L)
      sv$d[[1L]]
    }
  )
})

# Fallback for plain matrix — delegate to base::norm
setMethod("norm", "matrix", function(x, type = "F", ...) {
  base::norm(x, type = type, ...)
})

# Vector norms (no base::norm equivalent for numeric vectors)
setMethod("norm", "numeric", function(x, type = "2", ...) {
  type <- match.arg(toupper(substr(type, 1L, 1L)), c("2", "1", "I"))
  switch(type,
    "2" = sqrt(sum(x * x)),
    "1" = sum(abs(x)),
    "I" = max(abs(x))
  )
})
