ridge_path <- function(X, Y, lambdas, k = NULL, ...) {
  stopifnot(is.numeric(lambdas), length(lambdas) >= 1L, all(lambdas > 0))

  Y_mat <- if (is.null(dim(Y))) matrix(Y, ncol = 1L) else as.matrix(Y)
  stopifnot(nrow(Y_mat) == nrow(X))

  if (!inherits(X, "aMatrix")) X <- as_adgeMatrix(X)
  if (is.null(k)) k <- min(nrow(X), ncol(X))

  fac <- svd_factor(X, k = k, ...)

  k_use <- fac@k
  U_k <- fac@u[, seq_len(k_use), drop = FALSE]
  d_k <- fac@d[seq_len(k_use)]
  V_k <- fac@v[, seq_len(k_use), drop = FALSE]

  UtY <- base::crossprod(U_k, Y_mat)

  p  <- ncol(X)
  q  <- ncol(Y_mat)
  nl <- length(lambdas)

  coef_array <- array(NA_real_, dim = c(p, q, nl))

  for (i in seq_len(nl)) {
    d_lam <- d_k / (d_k^2 + lambdas[[i]])
    coef_array[,, i] <- V_k %*% (d_lam * UtY)
  }

  structure(
    list(coef = coef_array, lambdas = lambdas, svd = fac, k = k_use),
    class = "ridge_path"
  )
}

print.ridge_path <- function(x, ...) {
  cat(sprintf("ridge_path: %d lambda values, k=%d singular values\n",
    length(x$lambdas), x$k))
  cat(sprintf("  coef dim: [%s]\n", paste(dim(x$coef), collapse = " x ")))
  cat(sprintf("  lambda range: [%.4g, %.4g]\n",
    min(x$lambdas), max(x$lambdas)))
  invisible(x)
}

coef.ridge_path <- function(object, lambda = NULL, ...) {
  if (is.null(lambda)) return(object$coef)
  idx <- which.min(abs(object$lambdas - lambda))
  object$coef[,, idx, drop = FALSE]
}

predict.ridge_path <- function(object, newdata, lambda = NULL, ...) {
  coefs <- coef(object, lambda = lambda)
  newdata <- as.matrix(newdata)
  if (is.null(lambda)) {
    result <- array(NA_real_, dim = c(nrow(newdata), dim(coefs)[2L], dim(coefs)[3L]))
    for (i in seq_len(dim(coefs)[3L])) result[,, i] <- newdata %*% coefs[,, i]
    result
  } else {
    newdata %*% coefs[,, 1L, drop = FALSE]
  }
}
