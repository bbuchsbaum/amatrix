setMethod("%*%", signature(x = "adgeMatrix", y = "ANY"), function(x, y) am_matmul(x, y))
setMethod("%*%", signature(x = "adgeMatrix", y = "matrix"), function(x, y) am_matmul(x, y))
setMethod("%*%", signature(x = "adgeMatrix", y = "Matrix"), function(x, y) am_matmul(x, y))
setMethod("%*%", signature(x = "adgeMatrix", y = "dgeMatrix"), function(x, y) am_matmul(x, y))
setMethod("%*%", signature(x = "adgeMatrix", y = "dgCMatrix"), function(x, y) am_matmul(x, y))
setMethod("%*%", signature(x = "adgeMatrix", y = "adgeMatrix"), function(x, y) am_matmul(x, y))
setMethod("%*%", signature(x = "adgeMatrix", y = "adgCMatrix"), function(x, y) am_matmul(x, y))

if (!isGeneric("t")) setGeneric("t", function(x) base::t(x))
setMethod("t", "adgeMatrix", function(x) am_transpose(x))

setMethod("crossprod", signature(x = "adgeMatrix", y = "ANY"), function(x, y = NULL, ...) am_crossprod(x, y = y, ...))
setMethod("crossprod", signature(x = "adgeMatrix", y = "missing"), function(x, y, ...) am_crossprod(x, y = NULL, ...))
setMethod("tcrossprod", signature(x = "adgeMatrix", y = "ANY"), function(x, y = NULL, ...) am_tcrossprod(x, y = y, ...))
setMethod("tcrossprod", signature(x = "adgeMatrix", y = "missing"), function(x, y, ...) am_tcrossprod(x, y = NULL, ...))

setMethod("rowSums", "adgeMatrix", function(x, na.rm = FALSE, dims = 1L) am_row_sums(x, na.rm = na.rm, dims = dims))
setMethod("colSums", "adgeMatrix", function(x, na.rm = FALSE, dims = 1L) am_col_sums(x, na.rm = na.rm, dims = dims))

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
  am_eigen(x, symmetric = symmetric, only.values = only.values, EISPACK = EISPACK)
})

setMethod("diag", "adgeMatrix", function(x = 1, nrow, ncol, names = TRUE) {
  am_diag(x, nrow = nrow, ncol = ncol, names = names)
})

setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "ANY"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "numeric"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "integer"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "logical"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "matrix"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "adgeMatrix", e2 = "Matrix"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "ANY", e2 = "adgeMatrix"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "numeric", e2 = "adgeMatrix"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "integer", e2 = "adgeMatrix"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "logical", e2 = "adgeMatrix"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "matrix", e2 = "adgeMatrix"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "Matrix", e2 = "adgeMatrix"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})
