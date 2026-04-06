setMethod("%*%", signature(x = "adgCMatrix", y = "ANY"), function(x, y) am_matmul(x, y))
setMethod("%*%", signature(x = "adgCMatrix", y = "matrix"), function(x, y) am_matmul(x, y))
setMethod("%*%", signature(x = "adgCMatrix", y = "Matrix"), function(x, y) am_matmul(x, y))
setMethod("%*%", signature(x = "adgCMatrix", y = "dgeMatrix"), function(x, y) am_matmul(x, y))
setMethod("%*%", signature(x = "adgCMatrix", y = "dgCMatrix"), function(x, y) am_matmul(x, y))
setMethod("%*%", signature(x = "adgCMatrix", y = "adgeMatrix"), function(x, y) am_matmul(x, y))
setMethod("%*%", signature(x = "adgCMatrix", y = "adgCMatrix"), function(x, y) am_matmul(x, y))

setMethod("t", "adgCMatrix", function(x) am_transpose(x))

setMethod("crossprod", signature(x = "adgCMatrix", y = "ANY"), function(x, y = NULL, ...) am_crossprod(x, y = y, ...))
setMethod("tcrossprod", signature(x = "adgCMatrix", y = "ANY"), function(x, y = NULL, ...) am_tcrossprod(x, y = y, ...))

setMethod("rowSums", "adgCMatrix", function(x, na.rm = FALSE, dims = 1L) am_rowsums(x, na.rm = na.rm, dims = dims))
setMethod("colSums", "adgCMatrix", function(x, na.rm = FALSE, dims = 1L) am_colsums(x, na.rm = na.rm, dims = dims))

setMethod("[", signature(x = "adgCMatrix", i = "ANY", j = "ANY", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  am_subset(x, i, j, ..., drop = drop)
})

setMethod("[", signature(x = "adgCMatrix", i = "index", j = "index", drop = "logical"), function(x, i, j, ..., drop = TRUE) {
  am_subset(x, i, j, ..., drop = drop)
})

setMethod("[", signature(x = "adgCMatrix", i = "missing", j = "index", drop = "logical"), function(x, i, j, ..., drop = TRUE) {
  am_subset(x, i, j, ..., drop = drop)
})

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

setMethod("solve", signature(a = "adgCMatrix", b = "missing"), function(a, b, ...) am_solve(a, ...))
setMethod("solve", signature(a = "adgCMatrix", b = "ANY"), function(a, b, ...) am_solve(a, b = b, ...))
setMethod("chol", "adgCMatrix", function(x, ...) am_chol(x, ...))
setMethod("qr", "adgCMatrix", function(x, ...) am_qr(x, ...))

setMethod("svd", "adgCMatrix", function(x, nu = min(dim(x)), nv = min(dim(x)), LINPACK = FALSE, ...) {
  am_svd(x, nu = nu, nv = nv, LINPACK = LINPACK, ...)
})

setMethod("eigen", "adgCMatrix", function(x, symmetric, only.values = FALSE, EISPACK = FALSE) {
  am_eigen(x, symmetric = symmetric, only.values = only.values, EISPACK = EISPACK)
})

setMethod("diag", "adgCMatrix", function(x = 1, nrow, ncol, names = TRUE) {
  am_diag(x, nrow = nrow, ncol = ncol, names = names)
})

setMethod("Ops", signature(e1 = "adgCMatrix", e2 = "ANY"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "adgCMatrix", e2 = "numeric"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "adgCMatrix", e2 = "integer"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "adgCMatrix", e2 = "logical"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "adgCMatrix", e2 = "matrix"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "adgCMatrix", e2 = "Matrix"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "ANY", e2 = "adgCMatrix"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "numeric", e2 = "adgCMatrix"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "integer", e2 = "adgCMatrix"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "logical", e2 = "adgCMatrix"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "matrix", e2 = "adgCMatrix"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})

setMethod("Ops", signature(e1 = "Matrix", e2 = "adgCMatrix"), function(e1, e2) {
  am_ewise(.Generic, e1, e2)
})
