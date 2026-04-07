#!/usr/bin/env Rscript
# examples/krr.R
#
# Kernel Ridge Regression (KRR) using amatrix.
#
# Demonstrates GPU-accelerated tcrossprod — the operation CCA does not
# exercise heavily.  KRR is n-driven: forming the kernel matrix K = X X^T / p
# scales as O(n^2 p), so the crossover where GPU wins moves with n, not p.
#
# Algorithm (dual / kernel form):
#   1. K      = X_train X_train^T / p      [n×n, GPU tcrossprod]
#   2. alpha  = (K + lambda*I)^{-1} y      [n×n chol + solve]
#   3. K_test = X_test  X_train^T / p      [n_test×n_train, GPU tcrossprod cross-case]
#   4. y_hat  = K_test %*% alpha            [matmul]
#
# Equivalence with primal ridge:
#   beta = (X^T X / p + lambda I)^{-1} X^T y / p
#   X_test %*% beta  ==  K_test %*% alpha   (Woodbury identity)
#
# Usage:
#   Rscript -e 'Sys.setenv(AMATRIX_MLX_PROBE_GPU="1"); source("examples/krr.R")'

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(amatrix)
  }
  if (requireNamespace("amatrix.mlx", quietly = TRUE))
    amatrix.mlx::amatrix_mlx_register()
  if (requireNamespace("amatrix.arrayfire", quietly = TRUE))
    amatrix.arrayfire::amatrix_arrayfire_register()
})

# ── 1.  Implementation ────────────────────────────────────────────────────────

krr_amatrix <- function(X_train, y_train, lambda = 1.0, backend = "auto") {
  X_mat <- as.matrix(X_train)
  y_mat <- as.matrix(y_train)
  storage.mode(X_mat) <- "double"
  storage.mode(y_mat) <- "double"

  n <- nrow(X_mat); p <- ncol(X_mat)
  bk <- if (identical(backend, "auto")) "cpu" else backend

  Xc <- adgeMatrix(X_mat, preferred_backend = bk, precision = "fast")

  # ── GPU: form n×n kernel matrix via tcrossprod ────────────────────────────
  # tcrossprod(Xc) = X X^T; stays on device until as.matrix() below.
  K <- tcrossprod(Xc) / p

  # ── Cholesky: (K + lambda*I) alpha = y  ──────────────────────────────────
  # Download K, add lambda to diagonal (no GPU op for diagonal update),
  # re-wrap as adgeMatrix so chol() dispatches to GPU when n >= 256.
  K_h <- as.matrix(K)
  diag(K_h) <- diag(K_h) + lambda
  Kc    <- adgeMatrix(K_h, preferred_backend = bk, precision = "fast")
  R     <- as.matrix(chol(Kc))              # GPU chol when n >= 256
  alpha <- backsolve(R, forwardsolve(t(R), y_mat))

  list(
    alpha   = alpha,
    X_train = X_mat,
    n_train = n,
    p       = p,
    lambda  = lambda,
    backend = bk
  )
}

krr_predict <- function(model, X_test) {
  X_test_mat <- as.matrix(X_test)
  storage.mode(X_test_mat) <- "double"
  bk <- model$backend

  Xtst <- adgeMatrix(X_test_mat,    preferred_backend = bk, precision = "fast")
  Xtrn <- adgeMatrix(model$X_train, preferred_backend = bk, precision = "fast")

  # ── GPU: K_test = X_test X_train^T / p  (cross-case tcrossprod) ──────────
  # n_test × n_train; exercises the p != q path we recently fixed.
  K_test <- tcrossprod(Xtst, Xtrn) / model$p

  as.matrix(K_test) %*% model$alpha
}

# ── 2.  Correctness check ─────────────────────────────────────────────────────
#
# Verify KRR dual-form predictions match the primal ridge solution:
#   beta_primal = solve(X^T X / p + lambda I, X^T y / p)
#   y_hat_primal = X_test %*% beta_primal
#
# The Woodbury identity guarantees these are equal.

cat("── Correctness check (dual vs primal) ───────────────────────────────────\n")

set.seed(42)
n_tr <- 300L; n_ts <- 80L; p <- 25L; lambda <- 0.1

# Synthetic regression: y = X %*% beta_true + noise
beta_true <- rnorm(p)
X_tr   <- matrix(rnorm(n_tr * p), n_tr, p)
X_ts   <- matrix(rnorm(n_ts * p), n_ts, p)
y_tr   <- X_tr %*% beta_true + rnorm(n_tr, sd = 0.5)

# Dual form (GPU-accelerated crossproduct path)
model  <- krr_amatrix(X_tr, y_tr, lambda = lambda, backend = "cpu")
yhat_dual <- krr_predict(model, X_ts)

# Primal form (direct CPU reference)
p_scale <- ncol(X_tr)
beta_primal <- solve(
  t(X_tr) %*% X_tr / p_scale + lambda * diag(p_scale),
  t(X_tr) %*% y_tr / p_scale
)
yhat_primal <- X_ts %*% beta_primal

max_err <- max(abs(yhat_dual - yhat_primal))
cat(sprintf("  Max |delta| dual vs primal: %.2e  %s\n", max_err,
            if (max_err < 1e-3) "[PASS]" else "[FAIL]"))

# ── 3.  Backend timing comparison ─────────────────────────────────────────────

cat("\n── Backend timing ───────────────────────────────────────────────────────\n")
cat("  tcrossprod threshold: max(n, p) >= 2048 triggers GPU\n\n")

backends_to_test <- "cpu"
if (requireNamespace("amatrix.mlx", quietly = TRUE) &&
    isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE)))
  backends_to_test <- c(backends_to_test, "mlx")
if (requireNamespace("amatrix.arrayfire", quietly = TRUE) &&
    isTRUE(try(amatrix.arrayfire::amatrix_arrayfire_is_available(), silent = TRUE)))
  backends_to_test <- c(backends_to_test, "arrayfire")

time_krr <- function(n, p, lambda, backend, reps = 3L) {
  set.seed(1L)
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  elapsed <- numeric(reps)
  for (i in seq_len(reps)) {
    t0 <- proc.time()["elapsed"]
    m  <- krr_amatrix(X, y, lambda = lambda, backend = backend)
    krr_predict(m, X)     # include prediction (exercises cross-case tcrossprod)
    elapsed[i] <- proc.time()["elapsed"] - t0
  }
  median(elapsed) * 1e3
}

for (cfg in list(
  list(n =  512L, p = 100L, lambda = 0.5, label = "small  (n=512,  p=100)"),
  list(n = 2048L, p = 100L, lambda = 0.5, label = "medium (n=2048, p=100)  [GPU on]"),
  list(n = 4096L, p = 100L, lambda = 0.5, label = "large  (n=4096, p=100)  [GPU on]"),
  list(n = 2048L, p = 300L, lambda = 0.5, label = "wide   (n=2048, p=300)  [GPU on, more work/elem]")
)) {
  cat(sprintf("  %s\n", cfg$label))
  cpu_ms <- NULL
  for (bk in backends_to_test) {
    ms <- tryCatch(
      time_krr(cfg$n, cfg$p, cfg$lambda, bk),
      error = function(e) { cat("    [", bk, "error:", conditionMessage(e), "]\n"); NA_real_ }
    )
    if (is.na(ms)) next
    if (bk == "cpu") cpu_ms <- ms
    speedup <- if (!is.null(cpu_ms) && bk != "cpu")
      sprintf("  (%.1fx vs cpu)", cpu_ms / ms) else ""
    cat(sprintf("    %-12s %7.1f ms%s\n", paste0(bk, ":"), ms, speedup))
  }
  cat("\n")
}

# ── 4.  What the crossover looks like ─────────────────────────────────────────

cat("── Note on crossover ────────────────────────────────────────────────────\n")
cat("  KRR is n-driven: kernel matrix K = X X^T / p has n^2 entries.\n")
cat("  GPU wins when n >= 2048 (default tcrossprod threshold).\n")
cat("  Unlike CCA (p-driven), increasing p gives more work per element\n")
cat("  and improves GPU efficiency at the same n.\n")
cat("  The prediction step also uses cross-case tcrossprod(X_test, X_train)\n")
cat("  (n_test != n_train in general), which exercises the fixed non-square path.\n")
cat("  Lower the threshold:\n")
cat("    options(amatrix.mlx.tcrossprod_min_dim = 256L)\n")
cat("    options(amatrix.arrayfire.tcrossprod_min_dim = 256L)\n")
