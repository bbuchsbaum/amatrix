#' GPU-native randomized SVD (Halko et al. 2011)
#'
#' Computes a truncated SVD via randomized projection entirely on the GPU.
#' All QR, matmul, and SVD steps stay on device; a single \code{mlx_eval}
#' materializes the results.  Falls back to \code{irlba::svdr} on CPU if no
#' GPU backend with rsvd support is active.
#'
#' @param x   An \code{adgeMatrix} or plain numeric matrix.
#' @param k   Number of singular values/vectors to compute.
#' @param n_oversamples  Extra columns for the random projection (default 10).
#'   Increasing this improves accuracy at modest cost.
#' @param n_iter  Number of power-iteration passes (default 2).
#'   More passes give better accuracy for matrices with slowly decaying spectra.
#' @param ...  Ignored (for forward compatibility).
#'
#' @return A list with components \code{u} (m x k), \code{d} (length-k
#'   singular values, decreasing), and \code{v} (n x k).
#'
#' @references Halko, N., Martinsson, P. G., & Tropp, J. A. (2011).
#'   Finding structure with randomness: Probabilistic algorithms for
#'   constructing approximate matrix decompositions.
#'   \emph{SIAM Review}, 53(2), 217-288.
#'
#' @export
rsvd <- function(x, k, n_oversamples = 10L, n_iter = 2L, ...) {
  k             <- as.integer(k)
  n_oversamples <- as.integer(n_oversamples)
  n_iter        <- as.integer(n_iter)

  if (is(x, "adgeMatrix")) {
    bk_name <- tryCatch(.amatrix_svd_factor_rsvd_backend(x), error = function(e) NULL)
    if (!is.null(bk_name)) {
      bk <- tryCatch(.amatrix_get_backend(bk_name), error = function(e) NULL)
      if (!is.null(bk) && is.function(bk$rsvd)) {
        return(bk$rsvd(x, k = k, n_oversamples = n_oversamples, n_iter = n_iter))
      }
    }
  }

  # CPU fallback: use irlba::svdr if available
  mat <- if (is(x, "adgeMatrix")) as.matrix(x) else x
  if (requireNamespace("irlba", quietly = TRUE)) {
    res <- irlba::svdr(mat, k = k, extra = n_oversamples, it = n_iter)
    return(list(u = res$u, d = res$d, v = res$v))
  }

  # Last resort: base svd (full decomposition — only for small matrices)
  res <- base::svd(mat, nu = k, nv = k)
  list(u = res$u[, seq_len(k), drop = FALSE],
       d = res$d[seq_len(k)],
       v = res$v[, seq_len(k), drop = FALSE])
}
