# dispatch-hardening.R
#
# Additional S4 methods that close dispatch gaps where a plain matrix or numeric
# is passed as the *left-hand* argument and an adgeMatrix appears on the right.
# Without these, S4 falls through to the base R generic which coerces the
# adgeMatrix to a plain matrix, silently destroying GPU residency.
#
# Covered gaps:
#   crossprod(matrix,  adgeMatrix) — x(k×m) → t(x)(m×k), result = m×k %*% m×n = k×n
#   tcrossprod(matrix, adgeMatrix) — x(k×m) %*% t(y)(n×m) = k×n
#   crossprod(numeric, adgeMatrix) — treat numeric as column vector, result = 1×n
#
# Implementation strategy: wrap the plain-matrix argument in a temporary
# adgeMatrix sharing y's preferred_backend so the resident path is available.
# The temporary object is short-lived (no object_id binding needed beyond the
# call) and carries no resident key, so no residency leak occurs.

.amatrix_wrap_for_sparse_rhs <- function(x, y, vector_as = c("row", "col")) {
  stopifnot(inherits(y, "adgCMatrix"))
  vector_as <- match.arg(vector_as)

  if (inherits(x, "aMatrix")) {
    return(x)
  }

  if (is.numeric(x) && is.null(dim(x))) {
    x <- if (identical(vector_as, "row")) {
      matrix(x, nrow = 1L)
    } else {
      matrix(x, ncol = 1L)
    }
  }

  if (inherits(x, "dgCMatrix")) {
    return(new_adgCMatrix(
      x,
      preferred_backend = y@preferred_backend,
      policy = y@policy,
      precision = y@precision
    ))
  }

  new_adgeMatrix(
    as.matrix(x),
    preferred_backend = y@preferred_backend,
    policy = y@policy,
    precision = y@precision
  )
}

# ── crossprod(matrix, adgeMatrix) ─────────────────────────────────────────────
# crossprod(A, B) = t(A) %*% B
# A: k×m plain matrix   B: m×n adgeMatrix   result: k×n
#' @noRd
setMethod("crossprod",
  signature(x = "matrix", y = "adgeMatrix"),
  function(x, y, ...) {
    x_wrapped <- new_adgeMatrix(x, preferred_backend = y@preferred_backend,
                                policy = y@policy, precision = y@precision)
    crossprod(x_wrapped, y, ...)
  }
)

# ── tcrossprod(matrix, adgeMatrix) ────────────────────────────────────────────
# tcrossprod(A, B) = A %*% t(B)
# A: k×m plain matrix   B: n×m adgeMatrix   result: k×n
#' @noRd
setMethod("tcrossprod",
  signature(x = "matrix", y = "adgeMatrix"),
  function(x, y, ...) {
    x_wrapped <- new_adgeMatrix(x, preferred_backend = y@preferred_backend,
                                policy = y@policy, precision = y@precision)
    tcrossprod(x_wrapped, y, ...)
  }
)

# ── crossprod(numeric, adgeMatrix) ────────────────────────────────────────────
# Treats the numeric vector as a column vector (n×1), so t(v) %*% B = 1×p.
# This mirrors what base::crossprod does with a numeric y: crossprod(v, B)
# = t(as.matrix(v)) %*% B.
#' @noRd
setMethod("crossprod",
  signature(x = "numeric", y = "adgeMatrix"),
  function(x, y, ...) {
    x_mat <- matrix(x, ncol = 1L)
    x_wrapped <- new_adgeMatrix(x_mat, preferred_backend = y@preferred_backend,
                                policy = y@policy, precision = y@precision)
    crossprod(x_wrapped, y, ...)
  }
)

# ── Sparse RHS dispatch hardening ────────────────────────────────────────────
# Without these, `matrix %*% adgCMatrix` etc. fall through to Matrix/base,
# silently returning unwrapped dgeMatrix and losing GPU residency.

#' @noRd
setMethod("%*%", signature(x = "matrix",     y = "adgCMatrix"), function(x, y) {
  matmul(.amatrix_wrap_for_sparse_rhs(x, y, vector_as = "row"), y)
})
#' @noRd
setMethod("%*%", signature(x = "numeric",    y = "adgCMatrix"), function(x, y) {
  matmul(.amatrix_wrap_for_sparse_rhs(x, y, vector_as = "row"), y)
})
#' @noRd
setMethod("%*%", signature(x = "dgeMatrix",  y = "adgCMatrix"), function(x, y) {
  matmul(.amatrix_wrap_for_sparse_rhs(x, y, vector_as = "row"), y)
})
#' @noRd
setMethod("%*%", signature(x = "dgCMatrix",  y = "adgCMatrix"), function(x, y) {
  matmul(.amatrix_wrap_for_sparse_rhs(x, y, vector_as = "row"), y)
})
#' @noRd
setMethod("crossprod",  signature(x = "matrix",  y = "adgCMatrix"), function(x, y, ...) {
  am_crossprod(.amatrix_wrap_for_sparse_rhs(x, y, vector_as = "col"), y, ...)
})
#' @noRd
setMethod("crossprod",  signature(x = "numeric", y = "adgCMatrix"), function(x, y, ...) {
  am_crossprod(.amatrix_wrap_for_sparse_rhs(x, y, vector_as = "col"), y, ...)
})
#' @noRd
setMethod("crossprod",  signature(x = "dgeMatrix", y = "adgCMatrix"), function(x, y, ...) {
  am_crossprod(.amatrix_wrap_for_sparse_rhs(x, y, vector_as = "col"), y, ...)
})
#' @noRd
setMethod("crossprod",  signature(x = "dgCMatrix", y = "adgCMatrix"), function(x, y, ...) {
  am_crossprod(.amatrix_wrap_for_sparse_rhs(x, y, vector_as = "col"), y, ...)
})
#' @noRd
setMethod("tcrossprod", signature(x = "matrix",  y = "adgCMatrix"), function(x, y, ...) {
  am_tcrossprod(.amatrix_wrap_for_sparse_rhs(x, y, vector_as = "row"), y, ...)
})
#' @noRd
setMethod("tcrossprod", signature(x = "numeric", y = "adgCMatrix"), function(x, y, ...) {
  am_tcrossprod(.amatrix_wrap_for_sparse_rhs(x, y, vector_as = "row"), y, ...)
})
#' @noRd
setMethod("tcrossprod", signature(x = "dgeMatrix", y = "adgCMatrix"), function(x, y, ...) {
  am_tcrossprod(.amatrix_wrap_for_sparse_rhs(x, y, vector_as = "row"), y, ...)
})
#' @noRd
setMethod("tcrossprod", signature(x = "dgCMatrix", y = "adgCMatrix"), function(x, y, ...) {
  am_tcrossprod(.amatrix_wrap_for_sparse_rhs(x, y, vector_as = "row"), y, ...)
})
