#' Coerce amatrix objects to base R types
#'
#' Convert \code{adgeMatrix}, \code{adgCMatrix}, or \code{aTransposeView}
#' objects to base R \code{matrix}, numeric vector, or array by
#' materializing the host copy.
#'
#' @param x An \code{adgeMatrix}, \code{adgCMatrix}, or
#'   \code{aTransposeView}.
#' @param ... Further arguments passed to the corresponding base R
#'   coercion function.
#' @param mode Storage mode string passed to \code{as.vector}.
#'
#' @return A plain R \code{matrix}, numeric vector, or \code{array}
#'   containing the materialized host data.
#'
#' @examples
#' A <- adgeMatrix(matrix(1:6, 2, 3))
#' as.matrix(A)
#'
#' @rdname coerce-methods
#' @aliases as.matrix,adgeMatrix-method
setMethod("as.matrix", "adgeMatrix",      function(x, ...) as.matrix(amatrix_materialize_host(x), ...))
#' @rdname coerce-methods
#' @aliases as.matrix,adgCMatrix-method
setMethod("as.matrix", "adgCMatrix",      function(x, ...) as.matrix(amatrix_materialize_host(x), ...))
#' @rdname coerce-methods
#' @aliases as.matrix,aTransposeView-method
setMethod("as.matrix", "aTransposeView",  function(x, ...) t(as.matrix(amatrix_materialize_dense(x@source), ...)))

#' @rdname coerce-methods
#' @aliases as.numeric,adgeMatrix-method
setMethod("as.numeric", "adgeMatrix", function(x, ...) as.numeric(as.matrix(amatrix_materialize_host(x)), ...))
#' @rdname coerce-methods
#' @aliases as.vector,adgeMatrix-method
setMethod("as.vector",  "adgeMatrix", function(x, mode = "any") as.vector(as.matrix(amatrix_materialize_host(x)), mode))

#' @rdname coerce-methods
#' @aliases as.array,adgeMatrix-method
setMethod("as.array", "adgeMatrix", function(x, ...) as.array(as.matrix(amatrix_materialize_host(x)), ...))
#' @rdname coerce-methods
#' @aliases as.array,adgCMatrix-method
setMethod("as.array", "adgCMatrix", function(x, ...) as.array(as.matrix(amatrix_materialize_host(x)), ...))

setReplaceMethod("dimnames", "adgeMatrix", function(x, value) {
  am_set_dimnames(x, value)
})

setReplaceMethod("dimnames", "adgCMatrix", function(x, value) {
  am_set_dimnames(x, value)
})
