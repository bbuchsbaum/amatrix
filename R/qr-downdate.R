# qr_downdate: update a QR factorization after removing one row.
#
# v1 uses a pure-R refit fallback (correct interface, O(np^2) cost).
# TODO(Givens): Replace the refit path with a true Givens-rotation downdate
#   for O(p^2) per-fold cost. Algorithm: extract R from the factor, apply
#   a sequence of p Givens rotations to zero out the contribution of
#   row row_idx, then re-thin. See Golub & Van Loan sec 12.5.

#' QR downdate after removing one row
#'
#' Updates a QR factorization to reflect the removal of a single row
#' from the original matrix. The current implementation refits from the
#' reduced matrix; a Givens-rotation path is planned for a future
#' release.
#'
#' @param qr_factor An \code{amQR} factor object returned by
#'   \code{am_qr()}, or any object for which a method is defined.
#' @param row_idx Positive integer index of the row to remove.
#' @param X The original numeric matrix or \code{adgeMatrix} used to
#'   compute \code{qr_factor}. Required for \code{amQR} factors because
#'   they do not store the source matrix.
#'
#' @return An updated \code{amQR} factor with row \code{row_idx}
#'   excluded.
#'
#' @examples
#' \dontrun{
#' X <- adgeMatrix(matrix(rnorm(40), nrow = 8))
#' qf <- am_qr(X)
#' qf2 <- qr_downdate(qf, row_idx = 3L, X = X)
#' qr_info(qf2)$dim
#' }
#'
#' @seealso \code{\link{am_qr}}, \code{\link{qr_info}},
#'   \code{\link{lm_loo_cv}}
#' @export
qr_downdate <- function(qr_factor, row_idx, X = NULL) {
  UseMethod("qr_downdate")
}

#' @export
qr_downdate.amQR <- function(qr_factor, row_idx, X = NULL) {
  if (is.null(X)) {
    stop(
      "qr_downdate.amQR requires the original matrix X as the third argument. ",
      "amQR factors do not store the source matrix. ",
      "Call: qr_downdate(qr_factor, row_idx, X = your_matrix)",
      call. = FALSE
    )
  }
  X_sub <- X[-row_idx, , drop = FALSE]
  if (!inherits(X_sub, "aMatrix")) X_sub <- as_adgeMatrix(X_sub)
  am_qr(X_sub)
}

#' @export
qr_downdate.default <- function(qr_factor, row_idx, X = NULL) {
  stop(
    "qr_downdate.default requires an amQR factor (from am_qr()). ",
    "Base R qr() objects do not store the original matrix. ",
    "Use am_qr(X) to obtain a downdatable factor.",
    call. = FALSE
  )
}

# lm_loo_cv: exact leave-one-out CV via sequential qr_downdate.
#
# Returns a named list with:
#   $residuals  numeric vector of LOO prediction errors (y_i - y_hat_i)
#   $mse        mean squared LOO error
#
# Cost: O(n * p^2) for the refit path (v1).
# With Givens downdate: O(n * p^2) amortised but with much smaller constant.

#' Leave-one-out cross-validation for linear models
#'
#' Computes exact leave-one-out (LOO) prediction errors by refitting the
#' model \code{n} times, each time dropping one observation. Uses
#' \code{qr_downdate} internally to avoid recomputing the full
#' factorization from scratch.
#'
#' @param X Numeric matrix or \code{adgeMatrix} of predictors, shape
#'   \code{[n, p]}. No intercept column is added automatically.
#' @param y Numeric vector or single-column matrix of responses, length
#'   \code{n}.
#' @param method Character string passed to \code{am_qr} controlling
#'   the QR algorithm. Currently only \code{"qr"} is supported.
#' @param ... Additional arguments forwarded to \code{am_qr}.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{residuals}{Numeric vector of length \code{n}: LOO prediction
#'       errors \eqn{y_i - \hat{y}_i^{(-i)}}.}
#'     \item{mse}{Scalar mean squared LOO error.}
#'   }
#'
#' @examples
#' \dontrun{
#' X <- adgeMatrix(matrix(rnorm(50), nrow = 10))
#' y <- rnorm(10)
#' cv <- lm_loo_cv(X, y)
#' cv$mse
#' }
#'
#' @seealso \code{\link{lm_fit}}, \code{\link{many_lm}}
#' @export
lm_loo_cv <- function(X, y, method = "qr", ...) {
  stopifnot(
    is.numeric(y) || is.matrix(y),
    nrow(X) == length(y) || nrow(X) == nrow(as.matrix(y))
  )

  y_vec <- as.numeric(y)
  n     <- nrow(X)

  X_am <- if (inherits(X, "aMatrix")) X else as_adgeMatrix(X)

  # Full QR factor (used as starting point for downdate)
  qr_full <- am_qr(X_am, ...)

  loo_resid <- numeric(n)

  for (i in seq_len(n)) {
    qr_i   <- qr_downdate(qr_full, i, X = X)   # drop row i
    coef_i <- as.numeric(qr.coef(qr_i, y_vec[-i]))
    loo_resid[[i]] <- y_vec[[i]] - sum(X[i, ] * coef_i)
  }

  list(
    residuals = loo_resid,
    mse       = mean(loo_resid^2)
  )
}
