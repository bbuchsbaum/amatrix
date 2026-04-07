#!/usr/bin/env Rscript
# tools/benchmark-lda.R
#
# Benchmarks Linear Discriminant Analysis across CPU / MLX / ArrayFire.
#
# Algorithm (Fisher / Sw^{-1}Sb eigenvalue form):
#   1. Class means   via segment_mean             [GPU grouped reduction]
#   2. Within-class scatter Sw  via crossprod(X_c)   [GPU p×p gemm]
#   3. Between-class scatter Sb via weighted outer   [CPU, small K×p]
#   4. Chol(Sw) → transform Sb → eigen              [CPU, p×p]
#   5. Project X onto discriminant axes              [GPU n×p gemm]
#
# Key primitives exercised:
#   segment_mean, crossprod, am_sweep, matmul, chol, solve
#
# Correctness reference: MASS::lda() (if available) or manual base-R path.
#
# Run:
#   Rscript -e 'Sys.setenv(AMATRIX_MLX_PROBE_GPU="1"); \
#     setwd("/Users/bbuchsbaum/code/amatrix"); \
#     source("tools/benchmark-lda.R")'

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

# ── Available backends ────────────────────────────────────────────────────────

active_backends <- "cpu"
mlx_ok <- requireNamespace("amatrix.mlx", quietly = TRUE) &&
  isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE))
af_ok  <- requireNamespace("amatrix.arrayfire", quietly = TRUE) &&
  isTRUE(try(amatrix.arrayfire::amatrix_arrayfire_is_available(), silent = TRUE))
if (mlx_ok) active_backends <- c(active_backends, "mlx")
if (af_ok)  active_backends <- c(active_backends, "arrayfire")

cat("Active backends:", paste(active_backends, collapse = ", "), "\n\n")

# ── amatrix LDA implementation ────────────────────────────────────────────────
#
# Returns list(loadings [p×d], scores [n×d], values [d], means [K×p])
# where d = min(K-1, p, n_components).
#
# Ridge parameter `reg` stabilises near-singular Sw (default 1e-6 relative).

lda_amatrix <- function(X_mat, labels, n_components = NULL, backend = "cpu",
                         reg = 1e-6) {
  n <- nrow(X_mat); p <- ncol(X_mat)
  K <- max(labels)
  d <- min(K - 1L, p)
  if (!is.null(n_components)) d <- min(d, n_components)

  X_gpu  <- adgeMatrix(X_mat, preferred_backend = backend, precision = "fast")
  counts <- tabulate(labels, nbins = K)          # n_k per class

  # Step 1: class means [K×p] via GPU segment mean ─────────────────────────
  means_gpu <- segment_mean(X_gpu, labels, K) # adgeMatrix K×p
  means_mat <- as.matrix(means_gpu)

  # Step 2: within-class scatter Sw via one-pass algebraic identity ────────
  #   Sw = X^T X − C^T diag(n_k) C
  #
  # Derivation: Sw = sum_k (X_k - mu_k)^T(X_k - mu_k)
  #   = X^T X - X^T L C - C^T L^T X + C^T L^T L C
  # where L = one_hot(labels).  L^T X = segment_sum = diag(n_k) C, so:
  #   = X^T X - C^T diag(n_k) C
  #
  # Advantage: no n×p centering matrix — just one GPU crossprod(X) and a
  # tiny K×p CPU correction.  Eliminates the center_X bottleneck entirely.
  Xsq      <- as.matrix(crossprod(X_gpu))          # p×p via GPU (main GEMM)
  wt_means <- sqrt(counts) * means_mat              # K×p, row-wise sqrt(n_k)
  Sw       <- Xsq - crossprod(wt_means)             # p×p, CPU (small K)

  # Step 3: between-class scatter Sb = sum_k n_k (mu_k-mu)(mu_k-mu)^T ─────
  grand_mean <- colSums(means_mat * counts) / n    # weighted grand mean
  mu_diff    <- means_mat - rep(grand_mean, each = K)   # K×p deviations
  wt_diff    <- sqrt(counts) * mu_diff             # K×p, weighted
  Sb         <- crossprod(wt_diff)                 # p×p, CPU (small K)

  # Step 4: generalised eigenvalue via Cholesky transform ───────────────────
  #   Regularise Sw: Sw_reg = Sw + reg * tr(Sw)/p * I
  trace_Sw   <- sum(diag(Sw))
  Sw_reg     <- Sw + diag(reg * trace_Sw / p, p)
  R          <- chol(Sw_reg)                       # R^T R = Sw_reg
  #   Sb_tilde = R^{-T} Sb R^{-1}  (symmetric)
  Rinv_Sb    <- solve(t(R), Sb)
  Sb_tilde   <- solve(t(R), t(Rinv_Sb))
  Sb_sym     <- 0.5 * (Sb_tilde + t(Sb_tilde))
  ev         <- eigen(Sb_sym, symmetric = TRUE)
  #   Back-transform: W = R^{-1} * V
  W          <- solve(R, ev$vectors[, seq_len(d), drop = FALSE])
  W          <- apply(W, 2L, function(w) w / sqrt(sum(w^2)))   # unit norm

  # Step 5: project X onto discriminant axes ────────────────────────────────
  # Size gate: GPU matmul for n×p × p×d is only faster when n*d >= gemm_min^2.
  gemm_min <- as.integer(getOption("amatrix.arrayfire.bdc_gemm_min", 256L))
  if (n * d >= gemm_min * gemm_min) {
    W_gpu  <- adgeMatrix(W, preferred_backend = backend, precision = "fast")
    scores <- as.matrix(X_gpu %*% W_gpu)
  } else {
    scores <- X_mat %*% W
  }

  list(loadings   = W,
       scores     = scores,
       values     = ev$values[seq_len(d)],
       means      = means_mat,
       grand_mean = grand_mean)
}

# Base-R reference (same algorithm, no GPU) ───────────────────────────────────

lda_base <- function(X_mat, labels, n_components = NULL, reg = 1e-6) {
  n <- nrow(X_mat); p <- ncol(X_mat)
  K <- max(labels)
  d <- min(K - 1L, p)
  if (!is.null(n_components)) d <- min(d, n_components)

  counts     <- tabulate(labels, nbins = K)
  means_mat  <- do.call(rbind, lapply(seq_len(K), function(k)
    colMeans(X_mat[labels == k, , drop = FALSE])))
  grand_mean <- colSums(means_mat * counts) / n
  # One-pass Sw: X^T X − C^T diag(n_k) C  (same formula as GPU path)
  wt_means <- sqrt(counts) * means_mat
  Sw       <- crossprod(X_mat) - crossprod(wt_means)
  mu_diff    <- means_mat - rep(grand_mean, each = K)
  Sb         <- crossprod(sqrt(counts) * mu_diff)
  trace_Sw   <- sum(diag(Sw))
  Sw_reg     <- Sw + diag(reg * trace_Sw / p, p)
  R          <- chol(Sw_reg)
  Rinv_Sb    <- solve(t(R), Sb)
  Sb_tilde   <- solve(t(R), t(Rinv_Sb))
  ev         <- eigen(0.5 * (Sb_tilde + t(Sb_tilde)), symmetric = TRUE)
  W          <- solve(R, ev$vectors[, seq_len(d), drop = FALSE])
  W          <- apply(W, 2L, function(w) w / sqrt(sum(w^2)))
  scores     <- X_mat %*% W

  list(loadings   = W,
       scores     = scores,
       values     = ev$values[seq_len(d)],
       means      = means_mat,
       grand_mean = grand_mean)
}

# ── Correctness checks ────────────────────────────────────────────────────────

check_lda <- function(n, p, K, backend, seed = 42L) {
  set.seed(seed)
  # Simulate K-class Gaussian data — balanced allocation so no class is empty
  class_means <- matrix(rnorm(K * p, sd = 2), K, p)
  base_labels <- rep(seq_len(K), ceiling(n / K))[seq_len(n)]
  labels      <- sample(base_labels)
  X_mat  <- class_means[labels, ] + matrix(rnorm(n * p), n, p)
  storage.mode(X_mat) <- "double"

  ref <- lda_base(X_mat, labels)
  am  <- lda_amatrix(X_mat, labels, backend = backend)

  # Eigenvalues (sign-invariant, should match up to numerical error)
  val_err <- max(abs(ref$values - am$values) / pmax(abs(ref$values), 1e-10))

  # Subspace angle: compare column spaces of W (sign/direction indeterminacy)
  # Use |cos theta| = |w_ref . w_am| >= threshold
  cos_angles <- abs(colSums(ref$loadings * am$loadings))
  subspace_ok <- all(cos_angles > 1 - 1e-3)

  # Class means from segment_mean
  mean_err <- max(abs(ref$means - am$means))

  # Reconstruction: scores should span same subspace (up to orthogonal rotation)
  # Check ||P_am - P_ref||_F where P = W W^T (projection matrices)
  P_ref  <- ref$loadings %*% t(ref$loadings)
  P_am   <- am$loadings  %*% t(am$loadings)
  proj_err <- norm(P_ref - P_am, "F") / norm(P_ref, "F")

  ok <- val_err < 1e-4 && mean_err < 1e-4 && proj_err < 1e-3

  cat(sprintf("  [%s] n=%d p=%d K=%d | eigenval_err=%.2e mean_err=%.2e proj_err=%.2e %s\n",
              backend, n, p, K, val_err, mean_err, proj_err,
              if (ok) "OK" else "FAIL"))
  ok
}

cat("── Correctness checks ──────────────────────────────────────────────────\n")
all_ok <- TRUE
for (bk in active_backends) {
  for (cfg in list(
    list(n = 200L,  p = 20L,  K = 4L),    # small standard
    list(n = 500L,  p = 50L,  K = 10L),   # medium multi-class
    list(n = 1000L, p = 5L,   K = 2L),    # binary, p << K-1
    list(n = 300L,  p = 100L, K = 3L),    # wide, p >> K
    list(n = 50L,   p = 8L,   K = 8L),    # n/K small (stress test)
    list(n = 2000L, p = 200L, K = 5L)     # larger
  )) {
    ok <- check_lda(cfg$n, cfg$p, cfg$K, bk)
    all_ok <- all_ok && ok
  }
}
cat(if (all_ok) "\nAll checks PASSED.\n" else "\nSome checks FAILED.\n")

# ── Performance benchmark ─────────────────────────────────────────────────────

if (!requireNamespace("bench", quietly = TRUE)) {
  cat("\nbench not available — skipping timing.\n")
  quit(save = "no")
}
library(bench)

bench_lda <- function(n, p, K, seed = 1L) {
  set.seed(seed)
  class_means <- matrix(rnorm(K * p, sd = 2), K, p)
  base_labels <- rep(seq_len(K), ceiling(n / K))[seq_len(n)]
  labels      <- sample(base_labels)
  X_mat  <- class_means[labels, ] + matrix(rnorm(n * p), n, p)
  storage.mode(X_mat) <- "double"

  exprs <- c(
    list(base = quote(lda_base(X_mat, labels))),
    setNames(lapply(active_backends, function(bk)
      substitute(lda_amatrix(X_mat, labels, backend = bk_), list(bk_ = bk))),
      active_backends)
  )

  res <- bench::mark(exprs = exprs, iterations = 5L,
                     check = FALSE, memory = FALSE, time_unit = "ms")
  res$size <- sprintf("%dx%d K=%d", n, p, K)
  res
}

cat("\n── Performance benchmarks ──────────────────────────────────────────────\n")
sizes <- list(
  list(n = 500L,  p = 40L,  K = 5L),
  list(n = 2000L, p = 100L, K = 8L),
  list(n = 8000L, p = 200L, K = 10L)
)

all_res <- lapply(sizes, function(s) {
  cat(sprintf("\n%dx%d  K=%d:\n", s$n, s$p, s$K))
  res <- bench_lda(s$n, s$p, s$K)
  print(res[, c("expression", "median", "itr/sec")])
  res
})

# ── Step profiling (largest size) ─────────────────────────────────────────────

cat("\n── Step profiling  (8000×200  K=10) ───────────────────────────────────\n")
{
  set.seed(99L); K <- 10L; n <- 8000L; p <- 200L
  class_means  <- matrix(rnorm(K * p, sd = 2), K, p)
  labels <- sample(rep(seq_len(K), ceiling(n / K))[seq_len(n)])
  X_mat  <- class_means[labels, ] + matrix(rnorm(n * p), n, p)
  counts     <- tabulate(labels, nbins = K)
  grand_mean <- colMeans(X_mat)
  bk <- active_backends[length(active_backends)]  # best available
  X_gpu <- adgeMatrix(X_mat, preferred_backend = bk, precision = "fast")

  steps <- bench::mark(
    segment_mean   = segment_mean(X_gpu, labels, K),
    # Old path: CPU centering (n×p) then GPU crossprod — was bottleneck
    Sw_twopass     = { mu <- as.matrix(segment_mean(X_gpu, labels, K));
                       Xc <- adgeMatrix(X_mat - mu[labels,], preferred_backend=bk,
                                        precision="fast");
                       crossprod(Xc) },
    # New path: one GPU crossprod(X) + tiny CPU correction — no n×p centering
    Sw_onepass     = { mu  <- as.matrix(segment_mean(X_gpu, labels, K));
                       Xsq <- as.matrix(crossprod(X_gpu));
                       Xsq - crossprod(sqrt(counts) * mu) },
    Sb_cpu         = { mu <- as.matrix(segment_mean(X_gpu, labels, K));
                       md <- mu - rep(colMeans(X_mat), each=K);
                       crossprod(sqrt(counts) * md) },
    project        = { W_dummy <- matrix(rnorm(p*9), p, 9);
                       X_gpu %*% adgeMatrix(W_dummy, preferred_backend=bk,
                                            precision="fast") },
    iterations = 5L, check = FALSE, memory = FALSE, time_unit = "ms"
  )
  print(steps[, c("expression", "median", "itr/sec")])
}

cat("\nDone.\n")
