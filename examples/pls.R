#!/usr/bin/env Rscript
# examples/pls.R
#
# Partial Least Squares Regression (PLS1 / NIPALS) using amatrix.
#
# Demonstrates GPU-accelerated matrix-vector products in an iterative
# algorithm.  Unlike CCA (dominated by large covariance crossproducts) and
# KRR (n-driven, O(n^2 p) kernel), PLS is p-driven: each NIPALS component
# costs three O(np) GPU matmuls and the bottleneck moves with max(n, p).
#
# Key implementation detail — implicit deflation:
#   Standard NIPALS deflates X in place each component, requiring a full
#   n×p re-upload to the GPU each iteration.  Instead we keep X_gpu fixed
#   (original matrix, uploaded once) and apply the accumulated deflation
#   as a small CPU correction:
#
#     X_k^T u  = X_0^T u  − P_{<k} (T_{<k}^T u)    [GPU + O(k·p) CPU]
#     X_k  w   = X_0  w   − T_{<k} (P_{<k}^T w)     [GPU + O(k·n) CPU]
#     X_k^T t  = X_0^T t  − P_{<k} (T_{<k}^T t)    [GPU + O(k·p) CPU]
#
#   For K components the corrections are O(K·(n+p)) total — negligible for
#   K << min(n, p).
#
# GPU dispatch threshold: max(n, p) >= 2048  (same as crossprod / CCA).
#
# Algorithm (NIPALS PLS1):
#   For k = 1 .. K:
#     1. w_k  = X_k^T r_{k-1} / ‖X_k^T r_{k-1}‖   [GPU crossprod, normalize]
#     2. t_k  = X_k w_k                              [GPU matmul]
#     3. p_k  = X_k^T t_k / (t_k^T t_k)             [GPU crossprod, normalize]
#     4. b_k  = t_k^T r_{k-1} / (t_k^T t_k)         [scalar]
#     5. Deflate residual: r_k = r_{k-1} − b_k t_k  [CPU, n-vector]
#   Prediction coefficients: β = W★ b, W★ = W (P^T W)^{-1}
#   Predict new X:           ŷ = X_test β
#
# Usage:
#   Rscript -e 'Sys.setenv(AMATRIX_MLX_PROBE_GPU="1"); source("examples/pls.R")'

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
  if (requireNamespace("amatrix.opencl", quietly = TRUE)) {
    Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = "1")
    options(amatrix.enable_opencl = TRUE)
    amatrix.opencl::amatrix_opencl_register()
  }
})

# ── 1.  Implementation ────────────────────────────────────────────────────────

pls_amatrix <- function(X_train, y_train, K = 10L, backend = "cpu") {
  X_mat <- as.matrix(X_train);  storage.mode(X_mat) <- "double"
  y_vec <- as.vector(y_train);  storage.mode(y_vec) <- "double"

  n <- nrow(X_mat); p <- ncol(X_mat)
  K <- min(K, p, n - 1L)  # can't have more components than rank

  # Allocate accumulators
  W <- matrix(0, p, K)   # weights     (unit-norm columns)
  P <- matrix(0, p, K)   # loadings
  T <- matrix(0, n, K)   # scores      (orthogonal columns)
  b <- numeric(K)         # regression coefficients for each score

  # Upload X once — implicit deflation avoids re-uploading each component.
  X_gpu <- adgeMatrix(X_mat, preferred_backend = backend, precision = "fast")
  r_cur <- y_vec  # running y residual

  for (k in seq_len(K)) {

    # ── (1) w_k = X_k^T r / ‖X_k^T r‖ ──────────────────────────────────────
    # GPU: X_0^T r_cur (p-vector)
    xt_r <- as.vector(crossprod(X_gpu, r_cur))
    # CPU correction for accumulated deflation: subtract P_{<k} (T_{<k}^T r)
    if (k > 1L) {
      Tk <- T[, seq_len(k - 1L), drop = FALSE]
      Pk <- P[, seq_len(k - 1L), drop = FALSE]
      xt_r <- xt_r - as.vector(Pk %*% (t(Tk) %*% r_cur))
    }
    w_k <- xt_r / sqrt(sum(xt_r^2))

    # ── (2) t_k = X_k w_k ────────────────────────────────────────────────────
    # GPU: X_0 w_k (n-vector)
    t_vec <- as.vector(X_gpu %*% w_k)
    # CPU correction: subtract T_{<k} (P_{<k}^T w_k)
    if (k > 1L) {
      t_vec <- t_vec - as.vector(Tk %*% (t(Pk) %*% w_k))
    }

    tt <- sum(t_vec^2)

    # ── (3) p_k = X_k^T t_k / (t_k^T t_k) ───────────────────────────────────
    # GPU: X_0^T t_vec (p-vector)
    xt_t <- as.vector(crossprod(X_gpu, t_vec))
    # CPU correction
    if (k > 1L) {
      xt_t <- xt_t - as.vector(Pk %*% (t(Tk) %*% t_vec))
    }
    p_k <- xt_t / tt

    # ── (4) scalar regression coefficient ────────────────────────────────────
    b[k] <- sum(t_vec * r_cur) / tt

    # ── (5) deflate y residual ────────────────────────────────────────────────
    r_cur <- r_cur - b[k] * t_vec

    W[, k] <- w_k
    P[, k] <- p_k
    T[, k] <- t_vec
  }

  # ── Prediction coefficient vector: beta = W* b, W* = W (P^T W)^{-1} ───────
  PtW    <- crossprod(P, W)      # K×K (small, CPU)
  W_star <- W %*% solve(PtW)    # p×K
  beta   <- as.vector(W_star %*% b)

  list(beta = beta, W = W, P = P, T = T, b = b, K = K,
       backend = backend, n_train = n, p = p)
}

pls_predict <- function(model, X_test) {
  X_mat <- as.matrix(X_test); storage.mode(X_mat) <- "double"
  as.vector(X_mat %*% model$beta)
}

# ── 2.  Correctness checks ────────────────────────────────────────────────────
#
# Check A: training predictions via scores  T b  match  X β  (algebraic identity).
# Check B: T^T T is diagonal (orthogonal scores — a hallmark of NIPALS correctness).
# Check C: K-component PLS minimises training MSE more than (K-1)-component.
# Check D: for K = p on a small overdetermined problem, PLS predictions equal OLS.

cat("── Correctness checks ───────────────────────────────────────────────────\n")

set.seed(42)
n_tr <- 200L; p_sm <- 15L; K_sm <- 5L

beta_true <- rnorm(p_sm)
X_tr_sm   <- matrix(rnorm(n_tr * p_sm), n_tr, p_sm)
y_tr_sm   <- X_tr_sm %*% beta_true + rnorm(n_tr, sd = 0.5)

m_sm <- pls_amatrix(X_tr_sm, y_tr_sm, K = K_sm, backend = "cpu")

# Check A — consistency of two prediction formulas
yhat_scores <- as.vector(m_sm$T %*% m_sm$b)
yhat_beta   <- pls_predict(m_sm, X_tr_sm)
err_A <- max(abs(yhat_scores - yhat_beta))
cat(sprintf("  A  T b vs X beta:          max|delta| = %.2e  %s\n",
            err_A, if (err_A < 1e-8) "[PASS]" else "[FAIL]"))

# Check B — orthogonal scores
TtT     <- crossprod(m_sm$T)
off_diag <- max(abs(TtT[upper.tri(TtT)]))
cat(sprintf("  B  Score orthogonality:    max off-diag = %.2e  %s\n",
            off_diag, if (off_diag < 1e-8) "[PASS]" else "[FAIL]"))

# Check C — monotone training MSE
mse_k <- vapply(seq_len(K_sm), function(k) {
  mk <- pls_amatrix(X_tr_sm, y_tr_sm, K = k, backend = "cpu")
  mean((pls_predict(mk, X_tr_sm) - y_tr_sm)^2)
}, numeric(1))
monotone <- all(diff(mse_k) <= 0)
cat(sprintf("  C  MSE decreases with K:   %s  (MSE K=1..%d: %s)\n",
            if (monotone) "[PASS]" else "[FAIL]", K_sm,
            paste(sprintf("%.3f", mse_k), collapse = ", ")))

# Check D — PLS with K = p matches OLS on a small overdetermined problem
K_full  <- p_sm
m_full  <- pls_amatrix(X_tr_sm, y_tr_sm, K = K_full, backend = "cpu")
yhat_pls_full <- pls_predict(m_full, X_tr_sm)
beta_ols <- as.vector(solve(crossprod(X_tr_sm), crossprod(X_tr_sm, y_tr_sm)))
yhat_ols  <- as.vector(X_tr_sm %*% beta_ols)
err_D <- max(abs(yhat_pls_full - yhat_ols))
cat(sprintf("  D  K=p PLS vs OLS:         max|delta| = %.2e  %s\n\n",
            err_D, if (err_D < 1e-5) "[PASS]" else "[FAIL]"))

# ── 3.  Backend timing ────────────────────────────────────────────────────────

cat("── Backend timing ───────────────────────────────────────────────────────\n")
cat("  GPU dispatch: max(n, p) >= 2048  (per-component crossprod + matmul)\n")
cat("  Algorithm: implicit deflation — X uploaded once, K GPU calls per comp.\n\n")

backends_to_test <- "cpu"
if (requireNamespace("amatrix.mlx", quietly = TRUE) &&
    isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE)))
  backends_to_test <- c(backends_to_test, "mlx")
if (requireNamespace("amatrix.arrayfire", quietly = TRUE) &&
    isTRUE(try(amatrix.arrayfire::amatrix_arrayfire_is_available(), silent = TRUE)))
  backends_to_test <- c(backends_to_test, "arrayfire")
if (requireNamespace("amatrix.opencl", quietly = TRUE) &&
    isTRUE(try(amatrix.opencl::amatrix_opencl_is_available(), silent = TRUE)))
  backends_to_test <- c(backends_to_test, "opencl")

time_pls <- function(n, p, K, backend, reps = 3L) {
  set.seed(1L)
  X <- matrix(rnorm(n * p), n, p)
  y <- X %*% rnorm(p) + rnorm(n, sd = 0.5)
  elapsed <- numeric(reps)
  for (i in seq_len(reps)) {
    t0 <- proc.time()["elapsed"]
    m  <- pls_amatrix(X, y, K = K, backend = backend)
    pls_predict(m, X)
    elapsed[i] <- proc.time()["elapsed"] - t0
  }
  median(elapsed) * 1e3
}

for (cfg in list(
  list(n = 1024L, p = 100L, K = 10L, label = "small  (n=1024, p=100,  K=10)"),
  list(n = 2048L, p = 200L, K = 10L, label = "medium (n=2048, p=200,  K=10)  [GPU on]"),
  list(n = 4096L, p = 200L, K = 10L, label = "large  (n=4096, p=200,  K=10)  [GPU on]"),
  list(n = 2048L, p = 500L, K = 10L, label = "wide   (n=2048, p=500,  K=10)  [GPU on, p-driven]"),
  list(n = 2048L, p = 200L, K = 30L, label = "deep   (n=2048, p=200,  K=30)  [GPU on, 3× iters]")
)) {
  cat(sprintf("  %s\n", cfg$label))
  cpu_ms <- NULL
  for (bk in backends_to_test) {
    ms <- tryCatch(
      time_pls(cfg$n, cfg$p, cfg$K, bk),
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

# ── 4.  Notes ─────────────────────────────────────────────────────────────────

cat("── Notes ────────────────────────────────────────────────────────────────\n")
cat("  PLS is p-driven: each component costs 3 × O(np) GPU matmuls.\n")
cat("  Unlike KRR (n-driven, O(n^2 p)), GPU wins when max(n,p) >= 2048.\n")
cat("  Unlike CCA (batch crossproducts), PLS is sequential: K components\n")
cat("  cannot be parallelised across components.\n")
cat("  Implicit deflation means X is uploaded to the GPU exactly once\n")
cat("  regardless of K; per-iteration CPU corrections are O(K·(n+p)).\n")
cat("  Lower the threshold for small-n benchmarking:\n")
cat("    options(amatrix.mlx.crossprod_min_dim = 256L)\n")
