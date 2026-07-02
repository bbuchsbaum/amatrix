# Kronecker product methods for amatrix classes
#
# base::kronecker() and the %x% operator dispatch through the S4 generic that
# Matrix installs. For the parent dgeMatrix / dgCMatrix classes Matrix's own
# methods fire, returning plain (non-amatrix) results and dropping the
# backend-dispatch metadata. adgeMatrix / adgCMatrix therefore silently demote
# to dgeMatrix / dgCMatrix. These methods intercept the amatrix subclasses,
# compute the product on the underlying Matrix contents (via Matrix's own
# kronecker methods, so values match base::kronecker exactly), and re-wrap the
# result as an amatrix, preserving the backend metadata of the first amatrix
# operand. The result is sparse (adgCMatrix) when Matrix returns a sparse
# result and dense (adgeMatrix) otherwise.

#' @importFrom Matrix kronecker
NULL

# meta: preserve backend metadata from the first amatrix operand.
.am_kronecker_result <- function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
  meta <- if (inherits(X, "aMatrix")) X else Y
  xc <- if (inherits(X, "aMatrix")) amatrix_materialize_host(X) else X
  yc <- if (inherits(Y, "aMatrix")) amatrix_materialize_host(Y) else Y
  res <- kronecker(xc, yc, FUN = FUN, make.dimnames = make.dimnames, ...)
  if (inherits(res, "sparseMatrix")) {
    new_adgCMatrix(
      res,
      preferred_backend = meta@preferred_backend,
      policy = meta@policy,
      precision = meta@precision
    )
  } else {
    new_adgeMatrix(
      res,
      preferred_backend = meta@preferred_backend,
      policy = meta@policy,
      precision = meta@precision
    )
  }
}

#' Kronecker product of backend-aware matrices
#'
#' S4 methods for \code{\link[base]{kronecker}} and the \code{\%x\%} operator
#' that keep the result as an amatrix. Without them, \code{kronecker(A, B)} and
#' \code{A \%x\% B} dispatch to the \pkg{Matrix} methods for the parent
#' \code{dgeMatrix} / \code{dgCMatrix} classes and silently demote to a plain
#' (non-amatrix) result, discarding backend-dispatch metadata.
#'
#' The product itself is computed by \pkg{Matrix}'s own Kronecker methods on the
#' materialized host contents, so values are identical to
#' \code{base::kronecker()} on the dense contents. The result is re-wrapped as
#' an \code{\linkS4class{adgCMatrix}} when it is sparse and an
#' \code{\linkS4class{adgeMatrix}} otherwise, inheriting the preferred backend,
#' policy, and precision of the first amatrix operand.
#'
#' @param X,Y Kronecker factors. At least one is an \code{\linkS4class{aMatrix}}
#'   subclass; the other may be an amatrix, a base \code{matrix}, or a
#'   \pkg{Matrix} object.
#' @param FUN Function (or its name) applied to the outer products; passed to
#'   the underlying \pkg{Matrix} method. Defaults to \code{"*"}.
#' @param make.dimnames Logical; construct dimnames from the factors. Passed to
#'   the underlying method.
#' @param ... Further arguments passed to the underlying method.
#'
#' @return An \code{\linkS4class{adgeMatrix}} (dense) or
#'   \code{\linkS4class{adgCMatrix}} (sparse).
#'
#' @examples
#' A <- adgeMatrix(matrix(1:4, 2, 2))
#' B <- adgeMatrix(diag(2))
#' kronecker(A, B)
#' A %x% B
#'
#' @name kronecker-methods
#' @aliases kronecker,adgeMatrix,adgeMatrix-method
#'   kronecker,adgeMatrix,adgCMatrix-method
#'   kronecker,adgCMatrix,adgeMatrix-method
#'   kronecker,adgCMatrix,adgCMatrix-method
#'   kronecker,adgeMatrix,matrix-method
#'   kronecker,matrix,adgeMatrix-method
#'   kronecker,adgCMatrix,matrix-method
#'   kronecker,matrix,adgCMatrix-method
#' @exportMethod kronecker
setMethod("kronecker", signature("adgeMatrix", "adgeMatrix"), function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
  .am_kronecker_result(X, Y, FUN, make.dimnames, ...)
})

#' @rdname kronecker-methods
setMethod("kronecker", signature("adgeMatrix", "adgCMatrix"), function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
  .am_kronecker_result(X, Y, FUN, make.dimnames, ...)
})

#' @rdname kronecker-methods
setMethod("kronecker", signature("adgCMatrix", "adgeMatrix"), function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
  .am_kronecker_result(X, Y, FUN, make.dimnames, ...)
})

#' @rdname kronecker-methods
setMethod("kronecker", signature("adgCMatrix", "adgCMatrix"), function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
  .am_kronecker_result(X, Y, FUN, make.dimnames, ...)
})

#' @rdname kronecker-methods
setMethod("kronecker", signature("adgeMatrix", "matrix"), function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
  .am_kronecker_result(X, Y, FUN, make.dimnames, ...)
})

#' @rdname kronecker-methods
setMethod("kronecker", signature("matrix", "adgeMatrix"), function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
  .am_kronecker_result(X, Y, FUN, make.dimnames, ...)
})

#' @rdname kronecker-methods
setMethod("kronecker", signature("adgCMatrix", "matrix"), function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
  .am_kronecker_result(X, Y, FUN, make.dimnames, ...)
})

#' @rdname kronecker-methods
setMethod("kronecker", signature("matrix", "adgCMatrix"), function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
  .am_kronecker_result(X, Y, FUN, make.dimnames, ...)
})
