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

# ── Dense RHS dispatch hardening (Matrix-family LHS × adgeMatrix RHS) ────────
# Without these, `dgeMatrix %*% adgeMatrix` and friends fall through to the
# Matrix package's dispatch and silently coerce the adgeMatrix to a base
# matrix, destroying GPU residency. Pattern mirrors the sparse block above.
# Found by tools/audit-dispatch.R — do not remove without updating that script.

.amatrix_wrap_for_dense_rhs <- function(x, y, vector_as = c("row", "col")) {
  stopifnot(inherits(y, "adgeMatrix"))
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

  new_adgeMatrix(
    as.matrix(x),
    preferred_backend = y@preferred_backend,
    policy = y@policy,
    precision = y@precision
  )
}

#' @noRd
setMethod("%*%", signature(x = "dgeMatrix", y = "adgeMatrix"), function(x, y) {
  matmul(.amatrix_wrap_for_dense_rhs(x, y, vector_as = "row"), y)
})
#' @noRd
setMethod("%*%", signature(x = "dgCMatrix", y = "adgeMatrix"), function(x, y) {
  matmul(.amatrix_wrap_for_dense_rhs(x, y, vector_as = "row"), y)
})
#' @noRd
setMethod("crossprod", signature(x = "dgeMatrix", y = "adgeMatrix"), function(x, y, ...) {
  am_crossprod(.amatrix_wrap_for_dense_rhs(x, y, vector_as = "col"), y, ...)
})
#' @noRd
setMethod("crossprod", signature(x = "dgCMatrix", y = "adgeMatrix"), function(x, y, ...) {
  am_crossprod(.amatrix_wrap_for_dense_rhs(x, y, vector_as = "col"), y, ...)
})
#' @noRd
setMethod("tcrossprod", signature(x = "numeric", y = "adgeMatrix"), function(x, y, ...) {
  am_tcrossprod(.amatrix_wrap_for_dense_rhs(x, y, vector_as = "row"), y, ...)
})
#' @noRd
setMethod("tcrossprod", signature(x = "dgeMatrix", y = "adgeMatrix"), function(x, y, ...) {
  am_tcrossprod(.amatrix_wrap_for_dense_rhs(x, y, vector_as = "row"), y, ...)
})
#' @noRd
setMethod("tcrossprod", signature(x = "dgCMatrix", y = "adgeMatrix"), function(x, y, ...) {
  am_tcrossprod(.amatrix_wrap_for_dense_rhs(x, y, vector_as = "row"), y, ...)
})

# ── Signatures that Matrix otherwise wins ────────────────────────────────────
#
# S4 picks the method with the smallest total class distance. Matrix defines
# methods on (denseMatrix, vector), (denseMatrix, matrix),
# (denseMatrix, index, missing, replValue), ... and those beat amatrix's
# (adgeMatrix, ANY) fallbacks because "ANY" is further away than the
# two-step adgeMatrix -> denseMatrix inheritance. The inherited method then
# returns a plain dgeMatrix / dgCMatrix: values are right, but the amatrix
# class, its metadata and backend dispatch are silently lost. Which Matrix
# signatures exist differs between Matrix releases, so the pairs below are
# declared explicitly: each is exact on the amatrix argument and matches
# Matrix's class on the other, which makes the amatrix method strictly closer.

for (.cls in c("adgeMatrix", "adgCMatrix")) {
  setMethod("%*%", signature(x = .cls, y = "vector"), function(x, y) matmul(x, y))
  setMethod("crossprod", signature(x = .cls, y = "vector"),
            function(x, y = NULL, ...) am_crossprod(x, y = y, ...))
  setMethod("tcrossprod", signature(x = .cls, y = "vector"),
            function(x, y = NULL, ...) am_tcrossprod(x, y = y, ...))
}

for (.ycls in c("matrix", "denseMatrix", "CsparseMatrix")) {
  setMethod("crossprod", signature(x = "adgeMatrix", y = .ycls),
            function(x, y = NULL, ...) am_crossprod(x, y = y, ...))
  setMethod("tcrossprod", signature(x = "adgeMatrix", y = .ycls),
            function(x, y = NULL, ...) am_tcrossprod(x, y = y, ...))
}

for (.pair in list(c("adgeMatrix", "denseMatrix"), c("denseMatrix", "adgeMatrix"),
                   c("adgeMatrix", "CsparseMatrix"), c("CsparseMatrix", "adgeMatrix"),
                   c("adgCMatrix", "CsparseMatrix"), c("CsparseMatrix", "adgCMatrix"),
                   c("adgCMatrix", "denseMatrix"), c("denseMatrix", "adgCMatrix"))) {
  setMethod("kronecker", signature(.pair[[1L]], .pair[[2L]]),
            function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
              .am_kronecker_result(X, Y, FUN, make.dimnames, ...)
            })
}

# Sub-assignment. `[<-` dispatches on (x, i, j, value) and Matrix declares its
# methods with value = "replValue".
for (.cls in c("adgeMatrix", "adgCMatrix")) {
  setReplaceMethod("[", signature(x = .cls, i = "index", j = "missing", value = "replValue"),
    function(x, i, j, ..., value) {
      am_subassign(x, i, j, ..., value = value, .linear = missing(j) && nargs() <= 3L)
    })
  setReplaceMethod("[", signature(x = .cls, i = "missing", j = "index", value = "replValue"),
    function(x, i, j, ..., value) {
      am_subassign(x, i, j, ..., value = value, .linear = FALSE)
    })
  setReplaceMethod("[", signature(x = .cls, i = "index", j = "index", value = "replValue"),
    function(x, i, j, ..., value) {
      am_subassign(x, i, j, ..., value = value, .linear = FALSE)
    })
  # Matrix-valued subscripts (logical masks, two-column index matrices) are
  # always the one-subscript form.
  for (.icls in c("matrix", "Matrix")) {
    setReplaceMethod("[", signature(x = .cls, i = .icls, j = "missing", value = "replValue"),
      function(x, i, j, ..., value) {
        am_subassign(x, i, value = value, .linear = TRUE)
      })
  }
}
rm(.cls, .ycls, .pair, .icls)

# Matrix defines `Logic` group methods on its logical classes. A method on the
# more specific group generic is preferred over one on `Ops`, so the amatrix
# `Ops` methods for the logical wrappers never saw `&` / `|` and the result
# came back as a plain lgCMatrix / lgeMatrix.
for (.e1 in c("adlgCMatrix", "adlgeMatrix")) {
  for (.e2 in c("adlgCMatrix", "adlgeMatrix")) {
    setMethod("Logic", signature(e1 = .e1, e2 = .e2), function(e1, e2) {
      .amatrix_rewrap_value(e1, callGeneric(.amatrix_logical_host_arg(e1), .amatrix_logical_host_arg(e2)))
    })
  }
}
rm(.e1, .e2)
