# Woodbury matrix identity helpers.
#
# For a large system A and a rank-k update U C V, these functions compute
# (A + UCV)^{-1} b  and  log|A + UCV|  in O(nk^2 + k^3) time, avoiding
# the O(n^3) cost of refactorising the updated matrix.
#
# Reference: Woodbury, M.A. (1950). Inverting modified matrices.
#   Statistical Research Group Memo 42, Princeton.

#' Solve a linear system using the Woodbury matrix identity
#'
#' Computes \eqn{(A + UCV)^{-1} b} in \eqn{O(nk^2 + k^3)} time using the
#' Woodbury matrix identity, avoiding an \eqn{O(n^3)} refactorisation of
#' the updated matrix.
#'
#' @param A_factor An \code{amChol} object from \code{chol_factor()}, or a
#'   square numeric matrix that is automatically Cholesky-factored.
#' @param U Numeric matrix of shape \code{[n, k]}; low-rank left factor.
#' @param b Numeric matrix of shape \code{[n, rhs]}; right-hand side(s).
#' @param V Numeric matrix of shape \code{[k, n]}; low-rank right factor.
#'   Defaults to \code{t(U)} (symmetric update).
#' @param C_inv Numeric matrix of shape \code{[k, k]}; inverse of the
#'   central factor \eqn{C}. Defaults to \code{diag(k)} (pure rank-k
#'   update with \eqn{C = I}).
#'
#' @return Numeric matrix of shape \code{[n, rhs]}: the solution
#'   \eqn{(A + UCV)^{-1} b}.
#'
#' @examples
#' A <- crossprod(matrix(rnorm(25), 5)) + diag(5)
#' U <- matrix(rnorm(10), 5, 2)
#' b <- rnorm(5)
#' x <- woodbury_solve(A, U, b)
#' length(x)
#'
#' @seealso \code{\link{woodbury_logdet}}, \code{\link{chol_factor}}
#' @export
woodbury_solve <- function(A_factor, U, b, V = NULL, C_inv = NULL) {
  # A_factor : amChol (from chol_factor()) or a square matrix (auto-factored)
  # U        : n x k matrix (low-rank left factor)
  # b        : n x rhs right-hand sides
  # V        : k x n matrix; default t(U)  (symmetric case)
  # C_inv    : k x k matrix; default diag(k)  (identity = pure rank-k update)
  #
  # Returns: (A + UCV)^{-1} b  [n x rhs]

  # Auto-factor plain matrix input
  if (!inherits(A_factor, "amChol")) {
    A_factor <- chol_factor(as_adgeMatrix(as.matrix(A_factor)))
  }

  U <- as.matrix(U)
  b <- as.matrix(b)
  k <- ncol(U)

  if (is.null(V))     V     <- t(U)
  if (is.null(C_inv)) C_inv <- diag(k)

  V     <- as.matrix(V)
  C_inv <- as.matrix(C_inv)

  stopifnot(nrow(U) == nrow(b), nrow(V) == k, ncol(V) == nrow(U),
            nrow(C_inv) == k, ncol(C_inv) == k)

  # Step 1-2: A^{-1} b  and  A^{-1} U
  Ainv_b <- chol_solve(A_factor, b)   # n x rhs
  Ainv_U <- chol_solve(A_factor, U)   # n x k

  # Step 3: k x k inner matrix
  M <- C_inv + V %*% Ainv_U           # k x k

  # Step 4: solve M for V A^{-1} b
  inner <- solve(M, V %*% Ainv_b)     # k x rhs

  # Step 5: apply correction
  Ainv_b - Ainv_U %*% inner           # n x rhs
}

#' Log-determinant via the Woodbury matrix determinant lemma
#'
#' Computes \eqn{\log|A + UCV|} using the matrix determinant lemma
#' in \eqn{O(nk^2 + k^3)} time, reusing an existing Cholesky factor of
#' \eqn{A}.
#'
#' @param A_factor An \code{amChol} object from \code{chol_factor()}, or a
#'   square numeric matrix that is automatically Cholesky-factored.
#' @param U Numeric matrix of shape \code{[n, k]}; low-rank left factor.
#' @param V Numeric matrix of shape \code{[k, n]}; low-rank right factor.
#'   Defaults to \code{t(U)} (symmetric update).
#' @param C_inv Numeric matrix of shape \code{[k, k]}; inverse of the
#'   central factor \eqn{C}. Defaults to \code{diag(k)}.
#'
#' @return A length-1 numeric: \eqn{\log|A + UCV|}.
#'
#' @examples
#' A <- crossprod(matrix(rnorm(25), 5)) + diag(5)
#' U <- matrix(rnorm(10), 5, 2)
#' ld <- woodbury_logdet(A, U)
#' is.finite(ld)
#'
#' @seealso \code{\link{woodbury_solve}}, \code{\link{chol_factor}}
#' @export
woodbury_logdet <- function(A_factor, U, V = NULL, C_inv = NULL) {
  # Returns log|A + UCV| using the matrix determinant lemma:
  #   log|A + UCV| = log|A| + log|C^{-1}| + log|C^{-1} + V A^{-1} U|
  # (signs chosen so the result is the log-determinant of the updated matrix)

  if (!inherits(A_factor, "amChol")) {
    A_factor <- chol_factor(as_adgeMatrix(as.matrix(A_factor)))
  }

  U <- as.matrix(U)
  k <- ncol(U)

  if (is.null(V))     V     <- t(U)
  if (is.null(C_inv)) C_inv <- diag(k)

  V     <- as.matrix(V)
  C_inv <- as.matrix(C_inv)

  # log|A| from existing Cholesky
  logdet_A <- chol_logdet(A_factor)

  # A^{-1} U  [n x k]
  Ainv_U <- chol_solve(A_factor, U)

  # Inner matrix: C^{-1} + V A^{-1} U  [k x k]
  M <- C_inv + V %*% Ainv_U

  # log|C^{-1}| = -log|C|; use determinant() on the small k x k matrices
  logdet_Cinv  <- as.numeric(determinant(C_inv, logarithm = TRUE)$modulus)
  logdet_M     <- as.numeric(determinant(M,     logarithm = TRUE)$modulus)

  # Matrix determinant lemma:
  #   det(A + UCV) = det(A) * det(C) * det(C^{-1} + V A^{-1} U)
  # log|A + UCV| = log|A| + log|C| + log|M|
  #              = log|A| - log|C^{-1}| + log|M|
  logdet_A - logdet_Cinv + logdet_M
}
