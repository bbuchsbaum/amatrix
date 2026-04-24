#' Covariance-to-correlation methods for amatrix objects
#'
#' Bridge \code{cov2cor()} through Matrix's covariance-to-correlation methods
#' so standard workflows such as \code{cov2cor(crossprod(X))} keep working
#' when \code{crossprod()} preserves an amatrix class.
#'
#' @param V A square \code{adgeMatrix} or \code{adgCMatrix}.
#'
#' @return A base R correlation matrix, matching \code{stats::cov2cor()} on
#'   the corresponding host matrix.
#'
#' @examples
#' X <- adgeMatrix(matrix(1:9 + 0, 3, 3))
#' cov2cor(crossprod(X))
#'
#' @rdname cov2cor-methods
#' @aliases cov2cor,adgeMatrix-method
#' @exportMethod cov2cor
setMethod("cov2cor", "adgeMatrix", function(V) {
  stats::cov2cor(as.matrix(V))
})

#' @rdname cov2cor-methods
#' @aliases cov2cor,adgCMatrix-method
setMethod("cov2cor", "adgCMatrix", function(V) {
  stats::cov2cor(as.matrix(V))
})
