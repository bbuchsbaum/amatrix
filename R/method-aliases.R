#' Internal S4 method dispatch entries
#'
#' These S4 methods are dispatch implementations of standard generics
#' (`%*%`, `[`, `crossprod`, `t`, etc.) for amatrix classes. They are
#' not part of the public API -- use the generics directly. This help
#' page exists only to satisfy R CMD check.
#'
#' @aliases
#'   [,adgCMatrix,ANY,ANY,ANY-method
#'   [,adgCMatrix,index,index,logical-method
#'   [,adgCMatrix,index,missing,logical-method
#'   [,adgCMatrix,missing,index,logical-method
#'   [,adgeMatrix,ANY,ANY,ANY-method
#'   [,adgeMatrix,index,index,logical-method
#'   [,adgeMatrix,index,missing,logical-method
#'   [,adgeMatrix,missing,index,logical-method
#'   [<-,adgCMatrix,ANY,ANY,ANY-method
#'   [<-,adgCMatrix,index,index,Matrix-method
#'   [<-,adgCMatrix,index,index,integer-method
#'   [<-,adgCMatrix,index,index,logical-method
#'   [<-,adgCMatrix,index,index,matrix-method
#'   [<-,adgCMatrix,index,index,numeric-method
#'   [<-,adgeMatrix,ANY,ANY,ANY-method
#'   [<-,adgeMatrix,index,index,Matrix-method
#'   [<-,adgeMatrix,index,index,integer-method
#'   [<-,adgeMatrix,index,index,logical-method
#'   [<-,adgeMatrix,index,index,matrix-method
#'   [<-,adgeMatrix,index,index,numeric-method
#'   crossprod,KronMatrix,matrix-method
#'   crossprod,KronMatrix,missing-method
#'   crossprod,KronMatrix,numeric-method
#'   crossprod,dgCMatrix,adgCMatrix-method
#'   crossprod,dgeMatrix,adgCMatrix-method
#'   crossprod,matrix,adgCMatrix-method
#'   crossprod,matrix,adgeMatrix-method
#'   crossprod,numeric,adgCMatrix-method
#'   crossprod,numeric,adgeMatrix-method
#'   diag,adgCMatrix-method
#'   diag,adgeMatrix-method
#'   dim,KronMatrix-method
#'   dim,aTransposeView-method
#'   dimnames,aTransposeView-method
#'   dimnames<-,adgCMatrix,ANY-method
#'   dimnames<-,adgeMatrix,ANY-method
#'   eigen,matrix-method
#'   norm,adgCMatrix,ANY-method
#'   norm,adgeMatrix,ANY-method
#'   norm,matrix,ANY-method
#'   norm,numeric,ANY-method
#'   qr.Q,amQR-method
#'   qr.R,amQR-method
#'   qr.coef,amQR,ANY-method
#'   qr.fitted,amQR,ANY-method
#'   qr.qty,amQR,ANY-method
#'   qr.qy,amQR,ANY-method
#'   qr.resid,amQR,ANY-method
#'   qr.solve,amQR,ANY-method
#'   qr.solve,amQR,missing-method
#'   qr,adgCMatrix-method
#'   qr,adgeMatrix-method
#'   solve,KronMatrix,matrix-method
#'   solve,KronMatrix,missing-method
#'   solve,KronMatrix,numeric-method
#'   svd,matrix-method
#'   t,KronMatrix-method
#'   t,aTransposeView-method
#'   t,adgCMatrix-method
#'   t,adgeMatrix-method
#'   tcrossprod,dgCMatrix,adgCMatrix-method
#'   tcrossprod,dgeMatrix,adgCMatrix-method
#'   tcrossprod,matrix,adgCMatrix-method
#'   tcrossprod,matrix,adgeMatrix-method
#'   tcrossprod,numeric,adgCMatrix-method
#' @keywords internal
#' @name amatrix-methods-internal
NULL
