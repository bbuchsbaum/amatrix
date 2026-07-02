# KronMatrix — lazy Kronecker product A⊗B
#
# Stores factors A (m×n) and B (p×q) separately.  The full (mp × nq) matrix
# is never formed; all operations use the vec-permutation identity:
#
#   (A⊗B) vec(X) = vec(B X A^T)
#
# where X = matrix(y, nrow = q, ncol = n) for an input vector y of length nq.

#' Lazy Kronecker product of two matrices
#'
#' \code{KronMatrix} stores the two factor matrices \code{A} (m x n)
#' and \code{B} (p x q) without forming the full (mp x nq) Kronecker
#' product. Matrix-vector and matrix-matrix products are evaluated
#' using the vec-permutation identity
#' \code{(A x B) vec(X) = vec(B X t(A))}, keeping memory use at
#' \code{O(mn + pq)} rather than \code{O(mnpq)}.
#'
#' @slot A Numeric matrix; the left factor of the Kronecker product.
#' @slot B Numeric matrix; the right factor of the Kronecker product.
#'
#' @exportClass KronMatrix
#' @seealso \code{\link{kron_matrix}}
setClass(
  "KronMatrix",
  slots = list(A = "matrix", B = "matrix")
)

# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------

#' Construct a lazy Kronecker product
#'
#' Creates a \code{\linkS4class{KronMatrix}} that stores the factor
#' matrices \code{A} and \code{B} without materializing the full
#' Kronecker product. Standard operations such as \code{\%*\%},
#' \code{crossprod}, \code{solve}, and \code{determinant} are
#' available and exploit the Kronecker structure.
#'
#' @param A Numeric matrix or object coercible via \code{as.matrix()};
#'   the left Kronecker factor.
#' @param B Numeric matrix or object coercible via \code{as.matrix()};
#'   the right Kronecker factor.
#'
#' @return A \code{\linkS4class{KronMatrix}} of implicit dimensions
#'   \code{c(nrow(A) * nrow(B), ncol(A) * ncol(B))}.
#'
#' @examples
#' A <- matrix(1:4, 2, 2)
#' B <- diag(3)
#' K <- kron_matrix(A, B)
#' dim(K)
#' as.matrix(K)
#'
#' @export
kron_matrix <- function(A, B) {
  A_mat <- if (is.matrix(A)) A else as.matrix(A)
  B_mat <- if (is.matrix(B)) B else as.matrix(B)
  storage.mode(A_mat) <- "double"
  storage.mode(B_mat) <- "double"
  new("KronMatrix", A = A_mat, B = B_mat)
}

# ---------------------------------------------------------------------------
# Basic shape
# ---------------------------------------------------------------------------

#' @noRd
setMethod("dim", "KronMatrix", function(x) {
  c(nrow(x@A) * nrow(x@B), ncol(x@A) * ncol(x@B))
})

#' @noRd
setMethod("show", "KronMatrix", function(object) {
  dm <- dim(object)
  cat(sprintf(
    "KronMatrix [%d\u00d7%d] = (%d\u00d7%d) \u2297 (%d\u00d7%d)\n",
    dm[1L], dm[2L],
    nrow(object@A), ncol(object@A),
    nrow(object@B), ncol(object@B)
  ))
  invisible(object)
})

# ---------------------------------------------------------------------------
# Materialization
# ---------------------------------------------------------------------------

#' @rdname amatrix-s3-methods
#' @usage \method{as.matrix}{KronMatrix}(x, ...)
#' @export
as.matrix.KronMatrix <- function(x, ...) base::kronecker(x@A, x@B)

# ---------------------------------------------------------------------------
# Subsetting  (materialize-on-subset)
# ---------------------------------------------------------------------------

#' Subset a lazy Kronecker product
#'
#' Extracts elements of a \code{\linkS4class{KronMatrix}} using standard matrix
#' indexing (\code{K[i, j]}, \code{K[i, ]}, \code{K[, j]}, or linear
#' \code{K[i]}). Without this method \code{K[i, j]} fails with
#' \dQuote{object of type 'S4' is not subsettable}.
#'
#' \strong{Note:} this is a materialize-on-subset implementation. The full
#' \eqn{(mp \times nq)} Kronecker product is formed via
#' \code{\link[base]{kronecker}} before indexing, so subsetting does not
#' preserve the memory advantage of the lazy representation. It is intended for
#' convenient inspection rather than large-scale extraction.
#'
#' @param x A \code{\linkS4class{KronMatrix}}.
#' @param i,j Row and column subscripts, following base matrix semantics.
#' @param drop Logical; drop dimensions when the result has a single row or
#'   column. Defaults to \code{TRUE}, as for base matrices.
#' @param ... Unused.
#'
#' @return A numeric vector or matrix, exactly as base matrix indexing of
#'   \code{as.matrix(x)} would return.
#'
#' @examples
#' K <- kron_matrix(matrix(1:4, 2, 2), diag(2))
#' K[1, ]
#' K[, 2]
#' K[2, 3]
#'
#' @name KronMatrix-subset
#' @aliases [,KronMatrix,ANY,ANY,ANY-method
#' @keywords internal
setMethod("[", signature(x = "KronMatrix", i = "ANY", j = "ANY", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  full <- as.matrix(x)
  ndrop <- if (missing(drop)) 0L else 1L
  # Single-subscript form K[i] / K[]: no comma in the call (nargs counts only
  # x and i). Distinguishes linear indexing K[5] from row indexing K[2, ].
  if (missing(j) && (nargs() - ndrop) <= 2L) {
    if (missing(i)) {
      return(full)
    }
    return(full[i])
  }
  if (missing(i) && missing(j)) {
    return(full[, , drop = drop])
  }
  if (missing(i)) {
    return(full[, j, drop = drop])
  }
  if (missing(j)) {
    return(full[i, , drop = drop])
  }
  full[i, j, drop = drop]
})

# ---------------------------------------------------------------------------
# Transpose  t(A⊗B) = A^T ⊗ B^T
# ---------------------------------------------------------------------------

#' @noRd
setMethod("t", "KronMatrix", function(x) {
  kron_matrix(t(x@A), t(x@B))
})

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# (A⊗B) y  where y is a numeric vector of length ncol(A)*ncol(B)
.kron_matvec <- function(A, B, y) {
  X <- matrix(y, nrow = ncol(B), ncol = ncol(A))
  as.vector(B %*% X %*% t(A))
}

# (A⊗B) Y  where Y is a matrix — applies column by column
.kron_matmat <- function(A, B, Y) {
  nr <- nrow(B) * nrow(A)
  nc <- ncol(Y)
  out <- matrix(0.0, nr, nc)
  for (j in seq_len(nc)) out[, j] <- .kron_matvec(A, B, Y[, j])
  out
}

# ---------------------------------------------------------------------------
# Matrix multiplication  (A⊗B) %*% rhs
# ---------------------------------------------------------------------------

#' @noRd
setMethod("%*%", signature("KronMatrix", "numeric"), function(x, y) {
  .kron_matvec(x@A, x@B, y)
})

#' @noRd
setMethod("%*%", signature("KronMatrix", "matrix"), function(x, y) {
  .kron_matmat(x@A, x@B, y)
})

# lhs %*% (A⊗B)  ≡  t((A⊗B)^T %*% t(lhs))  = t((A^T⊗B^T) %*% t(lhs))
#' @noRd
setMethod("%*%", signature("numeric", "KronMatrix"), function(x, y) {
  as.vector(.kron_matmat(t(y@A), t(y@B), matrix(x, ncol = 1L)))
})

#' @noRd
setMethod("%*%", signature("matrix", "KronMatrix"), function(x, y) {
  t(.kron_matmat(t(y@A), t(y@B), t(x)))
})

# ---------------------------------------------------------------------------
# Crossproduct  crossprod(K) = (A^T A) ⊗ (B^T B);  crossprod(K, Y) = K^T Y
# ---------------------------------------------------------------------------

#' @noRd
setMethod("crossprod", signature("KronMatrix", "missing"), function(x, y = NULL) {
  kron_matrix(crossprod(x@A), crossprod(x@B))
})

#' @noRd
setMethod("crossprod", signature("KronMatrix", "matrix"), function(x, y) {
  .kron_matmat(t(x@A), t(x@B), y)
})

#' @noRd
setMethod("crossprod", signature("KronMatrix", "numeric"), function(x, y) {
  .kron_matvec(t(x@A), t(x@B), y)
})

# ---------------------------------------------------------------------------
# Solve  (A⊗B)^{-1} = A^{-1} ⊗ B^{-1}  (requires square A and B)
# ---------------------------------------------------------------------------

.kron_check_square <- function(K, call = sys.call(-1L)) {
  if (nrow(K@A) != ncol(K@A) || nrow(K@B) != ncol(K@B)) {
    stop("solve() requires square A and B factor matrices", call. = FALSE)
  }
}

#' @noRd
setMethod("solve", signature("KronMatrix", "missing"), function(a, b, ...) {
  .kron_check_square(a)
  kron_matrix(solve(a@A), solve(a@B))
})

#' @noRd
setMethod("solve", signature("KronMatrix", "numeric"), function(a, b, ...) {
  .kron_check_square(a)
  inv <- solve(a)
  .kron_matvec(inv@A, inv@B, b)
})

#' @noRd
setMethod("solve", signature("KronMatrix", "matrix"), function(a, b, ...) {
  .kron_check_square(a)
  inv <- solve(a)
  .kron_matmat(inv@A, inv@B, b)
})

# ---------------------------------------------------------------------------
# Determinant  det(A⊗B) = det(A)^nrow(B) * det(B)^nrow(A)
# ---------------------------------------------------------------------------

#' @noRd
setMethod("determinant", "KronMatrix", function(x, logarithm = TRUE, ...) {
  .kron_check_square(x)
  m <- nrow(x@A)   # dim of A factor
  n <- nrow(x@B)   # dim of B factor
  dA <- determinant(x@A, logarithm = TRUE)
  dB <- determinant(x@B, logarithm = TRUE)
  log_mod <- n * as.double(dA$modulus) + m * as.double(dB$modulus)
  sgn     <- (dA$sign ^ n) * (dB$sign ^ m)
  mod <- if (logarithm) log_mod else exp(log_mod)
  attr(mod, "logarithm") <- logarithm
  structure(list(modulus = mod, sign = sgn), class = "det")
})

# kron(): eager wrapper lives in wrappers.R — handles aMatrix inputs via
# .am_as_double_matrix(). kron_matrix() (above) is the lazy path.
