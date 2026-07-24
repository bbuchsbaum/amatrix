#' Logistic regression via IRLS, with amatrix-accelerated kernels
#'
#' Fits a binary logistic regression by iteratively reweighted least
#' squares (IRLS) — the same algorithm that powers
#' `glm(family = binomial())`. The per-iteration bottleneck, forming the
#' weighted normal equations, is expressed with `amatrix` kernels:
#'
#' * `amatrix::crossprod_weighted(X, w)` — \eqn{X^\top W X}
#' * `amatrix::xty_weighted(X, w, z)` — \eqn{X^\top W z}
#' * `solve()` on the resulting \eqn{p \times p} system
#'
#' The function contains **no backend-specific code**. Passed a plain base
#' matrix it computes on the CPU reference path. To accelerate it, wrap the
#' design matrix once at the boundary:
#'
#' ```r
#' X_gpu <- amatrix::adgeMatrix(X, mode = "fast")   # picks a healthy GPU backend
#' fit   <- logit_fit(X_gpu, y)
#' ```
#'
#' Everything else — the algorithm, the convergence test, the return value —
#' is identical. That is the adoption story this package demonstrates.
#'
#' @param X design matrix: a base `matrix`, or any `amatrix` dense matrix
#'   (e.g. the result of `amatrix::adgeMatrix()`). Should include an
#'   intercept column if one is wanted.
#' @param y binary response of length `nrow(X)`; values in \{0, 1\}.
#' @param max_iter maximum IRLS iterations.
#' @param tol convergence tolerance on the relative deviance change
#'   (same criterion as `glm`).
#' @param weight_floor lower bound applied to the IRLS weights
#'   \eqn{\mu(1-\mu)} for numerical stability near fitted probabilities of
#'   0 or 1.
#'
#' @return an object of class `"logit_fit"`: a list with elements
#'   `coefficients`, `fitted`, `deviance`, `iter`, `converged`.
#'
#' @section Package-author note:
#' `%*%` and the arithmetic operators dispatch on S4 operands automatically,
#' but `solve()` is a regular function — a package using it on `amatrix`
#' objects must import the method, e.g.
#' `@importMethodsFrom amatrix solve` (done here in this package).
#' @importMethodsFrom amatrix solve
#' @export
logit_fit <- function(X, y, max_iter = 25L, tol = 1e-8, weight_floor = 1e-10) {
  y <- as.numeric(y)
  if (!all(y %in% c(0, 1))) {
    stop("y must be a binary 0/1 response", call. = FALSE)
  }
  n <- nrow(X)
  p <- ncol(X)
  if (length(y) != n) {
    stop("length(y) must equal nrow(X)", call. = FALSE)
  }

  beta <- numeric(p)
  dev  <- Inf
  converged <- FALSE
  eta <- numeric(n)
  iter <- 0L

  for (iter in seq_len(max_iter)) {
    mu <- stats::plogis(eta)
    w  <- pmax(mu * (1 - mu), weight_floor)
    z  <- eta + (y - mu) / w

    # Hot kernels: one weighted crossprod (p x p), one weighted xty (p x 1),
    # one small solve. On an accelerated input the first two run on the GPU.
    xtwx <- amatrix::crossprod_weighted(X, w)
    xtwz <- amatrix::xty_weighted(X, w, z)
    beta <- as.numeric(as.matrix(solve(xtwx, xtwz)))
    if (!all(is.finite(beta))) {
      stop("IRLS produced non-finite coefficients ",
           "(singular or ill-conditioned X'WX)", call. = FALSE)
    }

    eta <- as.numeric(as.matrix(X %*% matrix(beta, ncol = 1L)))

    # Deviance computed from eta (not from clamped mu) for numerical
    # stability near saturated fits.
    dev_new <- -2 * sum(y * stats::plogis(eta, log.p = TRUE) +
                          (1 - y) * stats::plogis(-eta, log.p = TRUE))
    if (is.finite(dev) && abs(dev_new - dev) / (abs(dev_new) + 0.1) < tol) {
      dev <- dev_new
      converged <- TRUE
      break
    }
    dev <- dev_new
  }

  mu  <- stats::plogis(eta)
  eps <- 10 * .Machine$double.eps
  if (any(mu < eps | mu > 1 - eps)) {
    warning("fitted probabilities numerically 0 or 1 occurred; ",
            "possible separation, coefficients may be unbounded",
            call. = FALSE)
  }

  cn <- colnames(X)
  if (!is.null(cn)) names(beta) <- cn

  structure(
    list(
      coefficients = beta,
      fitted       = stats::plogis(eta),
      deviance     = dev,
      iter         = iter,
      converged    = converged
    ),
    class = "logit_fit"
  )
}

#' @export
print.logit_fit <- function(x, ...) {
  cat(sprintf(
    "Logistic regression (IRLS): %d coefficient(s), deviance %.4f, %d iteration(s)%s\n",
    length(x$coefficients), x$deviance, x$iter,
    if (x$converged) "" else " [NOT converged]"
  ))
  print(x$coefficients)
  invisible(x)
}
