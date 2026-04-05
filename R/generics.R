if (!isGeneric("rowSums")) {
  setGeneric("rowSums")
}

if (!isGeneric("colSums")) {
  setGeneric("colSums")
}

if (!isGeneric("solve")) {
  setGeneric("solve", function(a, b, ...) standardGeneric("solve"))
}

if (!isGeneric("chol")) {
  setGeneric("chol", function(x, ...) standardGeneric("chol"))
}

if (!isGeneric("qr")) {
  setGeneric("qr", function(x, ...) standardGeneric("qr"))
}

if (!isGeneric("qr.Q")) {
  setGeneric("qr.Q", function(qr, complete = FALSE, ...) standardGeneric("qr.Q"))
}

if (!isGeneric("qr.R")) {
  setGeneric("qr.R", function(qr, complete = FALSE, ...) standardGeneric("qr.R"))
}

if (!isGeneric("qr.solve")) {
  setGeneric("qr.solve", function(a, b, tol = 1e-07) standardGeneric("qr.solve"))
}

if (!isGeneric("qr.coef")) {
  setGeneric("qr.coef", function(qr, y) standardGeneric("qr.coef"))
}

if (!isGeneric("qr.qty")) {
  setGeneric("qr.qty", function(qr, y) standardGeneric("qr.qty"))
}

if (!isGeneric("qr.qy")) {
  setGeneric("qr.qy", function(qr, y) standardGeneric("qr.qy"))
}

if (!isGeneric("qr.fitted")) {
  setGeneric("qr.fitted", function(qr, y, k = NULL) standardGeneric("qr.fitted"))
}

if (!isGeneric("qr.resid")) {
  setGeneric("qr.resid", function(qr, y) standardGeneric("qr.resid"))
}

if (!isGeneric("svd")) {
  setGeneric("svd", function(x, nu = min(n, p), nv = min(n, p), LINPACK = FALSE, ...) standardGeneric("svd"))
}

if (!isGeneric("eigen")) {
  setGeneric("eigen", function(x, symmetric, only.values = FALSE, EISPACK = FALSE) standardGeneric("eigen"))
}

if (!isGeneric("diag")) {
  setGeneric("diag", function(x = 1, nrow, ncol, names = TRUE) standardGeneric("diag"))
}

if (!isGeneric("as.matrix")) {
  setGeneric("as.matrix", function(x, ...) standardGeneric("as.matrix"))
}

if (!isGeneric("as.array")) {
  setGeneric("as.array", function(x, ...) standardGeneric("as.array"))
}

if (!isGeneric("dimnames<-")) {
  setGeneric("dimnames<-", function(x, value) standardGeneric("dimnames<-"))
}
