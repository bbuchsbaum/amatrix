#' Internal S4 method dispatch entries
#'
#' These S4 methods are dispatch implementations of standard generics
#' (`%*%`, `[`, `crossprod`, `t`, etc.) for amatrix classes. They are
#' not part of the public API -- use the generics directly. This help
#' page exists only to satisfy R CMD check.
#'
#' @aliases
#'   !,adlgCMatrix-method
#'   !,adlgeMatrix-method
#'   %*%,dgCMatrix,adgeMatrix-method
#'   %*%,dgeMatrix,adgeMatrix-method
#'   [,adgCMatrix,ANY,ANY,ANY-method
#'   [,adgCMatrix,index,index,logical-method
#'   [,adgCMatrix,index,index,missing-method
#'   [,adgCMatrix,index,missing,logical-method
#'   [,adgCMatrix,index,missing,missing-method
#'   [,adgCMatrix,missing,index,logical-method
#'   [,adgCMatrix,missing,index,missing-method
#'   [,adgeMatrix,ANY,ANY,ANY-method
#'   [,adgeMatrix,index,index,logical-method
#'   [,adgeMatrix,index,index,missing-method
#'   [,adgeMatrix,index,missing,logical-method
#'   [,adgeMatrix,index,missing,missing-method
#'   [,adgeMatrix,missing,index,logical-method
#'   [,adgeMatrix,missing,index,missing-method
#'   [,aTransposeView,ANY,ANY,ANY-method
#'   [,aTransposeView,index,index,logical-method
#'   [,aTransposeView,index,index,missing-method
#'   [,aTransposeView,index,missing,logical-method
#'   [,aTransposeView,index,missing,missing-method
#'   [,aTransposeView,missing,index,logical-method
#'   [,aTransposeView,missing,index,missing-method
#'   [<-,adgCMatrix,ANY,ANY,ANY-method
#'   [<-,adgCMatrix,index,index,Matrix-method
#'   [<-,adgCMatrix,index,index,integer-method
#'   [<-,adgCMatrix,index,index,logical-method
#'   [<-,adgCMatrix,index,index,matrix-method
#'   [<-,adgCMatrix,index,index,numeric-method
#'   [<-,adgCMatrix,index,missing,ANY-method
#'   [<-,adgCMatrix,index,missing,Matrix-method
#'   [<-,adgCMatrix,index,missing,integer-method
#'   [<-,adgCMatrix,index,missing,logical-method
#'   [<-,adgCMatrix,index,missing,matrix-method
#'   [<-,adgCMatrix,index,missing,numeric-method
#'   [<-,adgCMatrix,missing,index,ANY-method
#'   [<-,adgCMatrix,missing,index,Matrix-method
#'   [<-,adgCMatrix,missing,index,integer-method
#'   [<-,adgCMatrix,missing,index,logical-method
#'   [<-,adgCMatrix,missing,index,matrix-method
#'   [<-,adgCMatrix,missing,index,numeric-method
#'   [<-,adgeMatrix,ANY,ANY,ANY-method
#'   [<-,adgeMatrix,index,index,Matrix-method
#'   [<-,adgeMatrix,index,index,integer-method
#'   [<-,adgeMatrix,index,index,logical-method
#'   [<-,adgeMatrix,index,index,matrix-method
#'   [<-,adgeMatrix,index,index,numeric-method
#'   [<-,adgeMatrix,index,missing,ANY-method
#'   [<-,adgeMatrix,matrix,missing,ANY-method
#'   [<-,adgeMatrix,Matrix,missing,ANY-method
#'   [<-,adgeMatrix,missing,index,ANY-method
#'   [<-,adgeMatrix,ndenseMatrix,missing,numeric-method
#'   [<-,adgeMatrix,ngeMatrix,missing,numeric-method
#'   crossprod,KronMatrix,matrix-method
#'   crossprod,KronMatrix,missing-method
#'   crossprod,KronMatrix,numeric-method
#'   crossprod,dgCMatrix,adgCMatrix-method
#'   crossprod,dgCMatrix,adgeMatrix-method
#'   crossprod,dgeMatrix,adgCMatrix-method
#'   crossprod,dgeMatrix,adgeMatrix-method
#'   crossprod,matrix,adgCMatrix-method
#'   crossprod,matrix,adgeMatrix-method
#'   crossprod,numeric,adgCMatrix-method
#'   crossprod,numeric,adgeMatrix-method
#'   determinant,adgeMatrix,logical-method
#'   determinant,adgeMatrix,missing-method
#'   diag,adgCMatrix-method
#'   diag,adgeMatrix-method
#'   dim,KronMatrix-method
#'   dim,aTransposeView-method
#'   dimnames,aTransposeView-method
#'   dimnames<-,adgCMatrix,ANY-method
#'   dimnames<-,adgeMatrix,ANY-method
#'   eigen,matrix-method
#'   %*%,adgCMatrix,vector-method
#'   %*%,adgeMatrix,vector-method
#'   %*%,numeric,aTransposeView-method
#'   [<-,adgCMatrix,Matrix,missing,replValue-method
#'   [<-,adgCMatrix,index,index,replValue-method
#'   [<-,adgCMatrix,index,missing,replValue-method
#'   [<-,adgCMatrix,matrix,missing,replValue-method
#'   [<-,adgCMatrix,missing,index,replValue-method
#'   [<-,adgeMatrix,Matrix,missing,replValue-method
#'   [<-,adgeMatrix,index,index,replValue-method
#'   [<-,adgeMatrix,index,missing,replValue-method
#'   [<-,adgeMatrix,matrix,missing,replValue-method
#'   [<-,adgeMatrix,missing,index,replValue-method
#'   as.numeric,aTransposeView-method
#'   colMeans,aTransposeView-method
#'   colSums,aTransposeView-method
#'   crossprod,aTransposeView,ANY-method
#'   crossprod,aTransposeView,missing-method
#'   crossprod,adgCMatrix,vector-method
#'   crossprod,adgeMatrix,CsparseMatrix-method
#'   crossprod,adgeMatrix,denseMatrix-method
#'   crossprod,adgeMatrix,matrix-method
#'   crossprod,adgeMatrix,vector-method
#'   diag,aTransposeView-method
#'   kronecker,CsparseMatrix,adgCMatrix-method
#'   kronecker,CsparseMatrix,adgeMatrix-method
#'   kronecker,adgCMatrix,CsparseMatrix-method
#'   kronecker,adgCMatrix,denseMatrix-method
#'   kronecker,adgeMatrix,CsparseMatrix-method
#'   kronecker,adgeMatrix,denseMatrix-method
#'   kronecker,denseMatrix,adgCMatrix-method
#'   kronecker,denseMatrix,adgeMatrix-method
#'   length,aTransposeView-method
#'   log,aTransposeView-method
#'   norm,aTransposeView,character-method
#'   norm,aTransposeView,missing-method
#'   rowMeans,aTransposeView-method
#'   rowSums,aTransposeView-method
#'   solve,aTransposeView,ANY-method
#'   solve,aTransposeView,missing-method
#'   tcrossprod,aTransposeView,ANY-method
#'   tcrossprod,aTransposeView,missing-method
#'   tcrossprod,adgCMatrix,vector-method
#'   tcrossprod,adgeMatrix,CsparseMatrix-method
#'   tcrossprod,adgeMatrix,denseMatrix-method
#'   tcrossprod,adgeMatrix,matrix-method
#'   tcrossprod,adgeMatrix,vector-method
#'   log,adgCMatrix-method
#'   log,adgeMatrix-method
#'   norm,adgCMatrix,character-method
#'   norm,adgCMatrix,missing-method
#'   norm,adgeMatrix,character-method
#'   norm,adgeMatrix,missing-method
#'   norm,numeric,character-method
#'   norm,numeric,missing-method
#'   qr,adgCMatrix-method
#'   qr,adgeMatrix-method
#'   qr.Q,amQR-method
#'   qr.R,amQR-method
#'   qr.coef,amQR,ANY-method
#'   qr.fitted,amQR,ANY-method
#'   qr.qty,amQR,ANY-method
#'   qr.qy,amQR,ANY-method
#'   qr.resid,amQR,ANY-method
#'   qr.solve,amQR,ANY-method
#'   qr.solve,amQR,missing-method
#'   solve,KronMatrix,matrix-method
#'   solve,KronMatrix,missing-method
#'   solve,KronMatrix,numeric-method
#'   t,KronMatrix-method
#'   t,aTransposeView-method
#'   t,adgCMatrix-method
#'   t,adgeMatrix-method
#'   tcrossprod,dgCMatrix,adgCMatrix-method
#'   tcrossprod,dgCMatrix,adgeMatrix-method
#'   tcrossprod,dgeMatrix,adgCMatrix-method
#'   tcrossprod,dgeMatrix,adgeMatrix-method
#'   tcrossprod,matrix,adgCMatrix-method
#'   tcrossprod,matrix,adgeMatrix-method
#'   tcrossprod,numeric,adgCMatrix-method
#'   tcrossprod,numeric,adgeMatrix-method
#' @keywords internal
#' @name amatrix-methods-internal
NULL

# ── Identity-safe unary / rounding methods ───────────────────────────────────
#
# Matrix implements unary minus, round(), signif() and zapsmall() as
# `x@x <- f(x@x); x`, which preserves every amatrix slot — including
# `object_id` and `finalizer_env`. The result would then be a *different*
# matrix sharing the identity of its source, so identity-keyed caches
# (Cholesky / QR / XtX) and the residency registry would hand back state that
# belongs to the original data. These methods route through the host value and
# mint a fresh wrapper instead.

.amatrix_fresh_like <- function(x, host_value) {
  if (inherits(x, "adgCMatrix")) {
    if (!inherits(host_value, "dgCMatrix")) {
      host_value <- methods::as(methods::as(methods::as(host_value, "dMatrix"), "generalMatrix"), "CsparseMatrix")
    }
    return(new_adgCMatrix(
      host_value,
      preferred_backend = x@preferred_backend,
      policy = x@policy,
      precision = x@precision
    ))
  }
  .amatrix_rewrap_like(x, host_value)
}

.amatrix_plain_host <- function(x) {
  host <- amatrix_materialize_host(x)
  if (inherits(x, "adgCMatrix")) {
    methods::as(host, "dgCMatrix")
  } else {
    out <- as.matrix(host)
    dimnames(out) <- dimnames(x)
    out
  }
}

#' @noRd
setMethod("Arith", signature(e1 = "adgeMatrix", e2 = "missing"), function(e1, e2) {
  if (identical(.Generic, "+")) return(e1)
  .amatrix_fresh_like(e1, methods::callGeneric(.amatrix_plain_host(e1)))
})

#' @noRd
setMethod("Arith", signature(e1 = "adgCMatrix", e2 = "missing"), function(e1, e2) {
  if (identical(.Generic, "+")) return(e1)
  .amatrix_fresh_like(e1, methods::callGeneric(.amatrix_plain_host(e1)))
})

#' @noRd
setMethod("Math2", signature(x = "adgeMatrix"), function(x, digits) {
  host <- .amatrix_plain_host(x)
  value <- if (missing(digits)) methods::callGeneric(host) else methods::callGeneric(host, digits)
  .amatrix_fresh_like(x, value)
})

#' @noRd
setMethod("Math2", signature(x = "adgCMatrix"), function(x, digits) {
  host <- .amatrix_plain_host(x)
  value <- if (missing(digits)) methods::callGeneric(host) else methods::callGeneric(host, digits)
  .amatrix_fresh_like(x, value)
})

# base::zapsmall() gained arguments in R 4.4, so the method definition takes
# its formals from the running R rather than hard-coding either signature.
.amatrix_zapsmall_method <- local({
  f <- function() NULL
  formals(f) <- formals(base::zapsmall)
  body(f) <- quote({
    cl <- match.call()
    cl[[1L]] <- quote(zapsmall)
    cl$x <- .amatrix_plain_host(x)
    .amatrix_fresh_like(x, eval(cl, parent.frame()))
  })
  f
})

#' @noRd
setMethod("zapsmall", signature(x = "adgeMatrix"), .amatrix_zapsmall_method)

#' @noRd
setMethod("zapsmall", signature(x = "adgCMatrix"), .amatrix_zapsmall_method)
