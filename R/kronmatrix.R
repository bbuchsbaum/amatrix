# KronMatrix — lazy Kronecker product A⊗B
#
# Stores factors A (m×n) and B (p×q) separately.  The full (mp × nq) matrix
# is never formed; all operations use the vec-permutation identity:
#
#   (A⊗B) vec(X) = vec(B X A^T)
#
# where X = matrix(y, nrow = q, ncol = n) for an input vector y of length nq.

setClass(
  "KronMatrix",
  slots = list(A = "matrix", B = "matrix")
)

# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------

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

setMethod("dim", "KronMatrix", function(x) {
  c(nrow(x@A) * nrow(x@B), ncol(x@A) * ncol(x@B))
})

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

#' @export
as.matrix.KronMatrix <- function(x, ...) base::kronecker(x@A, x@B)

# ---------------------------------------------------------------------------
# Transpose  t(A⊗B) = A^T ⊗ B^T
# ---------------------------------------------------------------------------

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

setMethod("%*%", signature("KronMatrix", "numeric"), function(x, y) {
  .kron_matvec(x@A, x@B, y)
})

setMethod("%*%", signature("KronMatrix", "matrix"), function(x, y) {
  .kron_matmat(x@A, x@B, y)
})

# lhs %*% (A⊗B)  ≡  t((A⊗B)^T %*% t(lhs))  = t((A^T⊗B^T) %*% t(lhs))
setMethod("%*%", signature("numeric", "KronMatrix"), function(x, y) {
  as.vector(.kron_matmat(t(y@A), t(y@B), matrix(x, ncol = 1L)))
})

setMethod("%*%", signature("matrix", "KronMatrix"), function(x, y) {
  t(.kron_matmat(t(y@A), t(y@B), t(x)))
})

# ---------------------------------------------------------------------------
# Crossproduct  crossprod(K) = (A^T A) ⊗ (B^T B);  crossprod(K, Y) = K^T Y
# ---------------------------------------------------------------------------

setMethod("crossprod", signature("KronMatrix", "missing"), function(x, y = NULL) {
  kron_matrix(crossprod(x@A), crossprod(x@B))
})

setMethod("crossprod", signature("KronMatrix", "matrix"), function(x, y) {
  .kron_matmat(t(x@A), t(x@B), y)
})

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

setMethod("solve", signature("KronMatrix", "missing"), function(a, b, ...) {
  .kron_check_square(a)
  kron_matrix(solve(a@A), solve(a@B))
})

setMethod("solve", signature("KronMatrix", "numeric"), function(a, b, ...) {
  .kron_check_square(a)
  inv <- solve(a)
  .kron_matvec(inv@A, inv@B, b)
})

setMethod("solve", signature("KronMatrix", "matrix"), function(a, b, ...) {
  .kron_check_square(a)
  inv <- solve(a)
  .kron_matmat(inv@A, inv@B, b)
})

# ---------------------------------------------------------------------------
# Determinant  det(A⊗B) = det(A)^nrow(B) * det(B)^nrow(A)
# ---------------------------------------------------------------------------

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
