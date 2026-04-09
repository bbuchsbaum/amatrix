#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(amatrix)
  }
  library(amatrix.models)
})

set.seed(101)
n <- 12L
p <- 5L
r <- 16L

X <- cbind(1, matrix(rnorm(n * (p - 1L)), nrow = n, ncol = p - 1L))
beta <- matrix(rnorm(p * r), nrow = p, ncol = r)
Y <- X %*% beta + matrix(rnorm(n * r, sd = 0.05), nrow = n, ncol = r)
Y_array <- array(Y[, seq_len(4L), drop = FALSE], dim = c(n, 2L, 2L))
Y_small <- Y[, seq_len(4L), drop = FALSE]

x_am <- adgeMatrix(X)

w <- seq(1, n) / n
many_fit <- many_lm(x_am, Y, method = "qr", include_residuals = TRUE, cache = TRUE)
many_fit_reuse <- many_lm(x_am, Y, method = "qr", include_residuals = FALSE, cache = TRUE)
many_weighted_fit <- many_lm(x_am, Y, weights = w, method = "qr", include_residuals = TRUE, cache = TRUE)
array_fit <- array_lm(x_am, Y_array, weights = w, method = "qr", include_residuals = TRUE)
wls <- wls_fit(x_am, Y_small, weights = w, method = "qr")
single_fit <- lm_fit(x_am, Y_small, method = "qr")
cov_mat <- covariance(x_am, weights = w)
cor_mat <- correlation(x_am)

cat("\nFlagship many-response QR fit\n")
cat(sprintf("responses: %d\n", many_fit$responses))
cat(sprintf("df.residual: %d\n", many_fit$df.residual))
cat(sprintf("qr helper path: %s\n", many_fit$qr_helper_path))
cat(sprintf("compact factor source: %s\n", many_fit$qr_compact_factor_source))
cat(sprintf("cache reused on second call: %s\n", many_fit_reuse$cache_reused))
cat("rss:\n")
print(many_fit$rss)

cat("\nWeighted many-response summary\n")
cat("weighted rss:\n")
print(many_weighted_fit$rss)

cat("\nSingle fit\n")
print(single_fit)

cat("\nWeighted fit summary\n")
print(wls)

cat("\nArray-response fit summary\n")
cat("response_dims:\n")
print(array_fit$response_dims)
cat("fitted dim:\n")
print(dim(fitted(array_fit)))
cat("residual dim:\n")
print(dim(residuals(array_fit)))

cat("\nCovariance summary\n")
print(dim(cov_mat))
cat("diag(covariance):\n")
print(diag(as.matrix(cov_mat)))

cat("\nCorrelation summary\n")
print(dim(cor_mat))
cat("diag(correlation):\n")
print(diag(as.matrix(cor_mat)))
