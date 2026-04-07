#!/usr/bin/env Rscript
# tools/benchmark-spectral.R
#
# Benchmarks normalized spectral clustering (Ng-Jordan-Weiss) across backends.
#
# Algorithm:
#   1. Affinity  W  = kernel_matrix(X, "rbf", sigma)       [n×n, GPU]
#   2. Degree    d  = rowSums(W)                            [GPU rowSums]
#   3. Normalize L  = D^{-1/2} W D^{-1/2}                  [2× GPU sweep]
#   4. Eigenvecs V  = rsvd(L, k=K)$u                        [GPU rsvd]
#   5. Row-norm  Vn = V / ||V[i,:]||                        [GPU sweep]
#   6. K-means on Vn via pairwise_sqdist_argmin + segment_mean
#
# Key primitives exercised (in combination):
#   kernel_matrix, rowSums, am_sweep, rsvd, pairwise_sqdist_argmin,
#   segment_mean — the entire GPU pipeline chained end-to-end.
#
# Correctness: planted Gaussian blobs; check cluster recovery accuracy
# (best-permutation match) at each backend.
#
# Run:
#   Rscript -e 'setwd("/Users/bbuchsbaum/code/amatrix"); \
#     source("tools/benchmark-spectral.R")'

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

active_backends <- "cpu"
mlx_ok <- requireNamespace("amatrix.mlx", quietly = TRUE) &&
  isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE))
af_ok  <- requireNamespace("amatrix.arrayfire", quietly = TRUE) &&
  isTRUE(try(amatrix.arrayfire::amatrix_arrayfire_is_available(), silent = TRUE))
if (mlx_ok) active_backends <- c(active_backends, "mlx")
if (af_ok)  active_backends <- c(active_backends, "arrayfire")

cat("Active backends:", paste(active_backends, collapse = ", "), "\n\n")

# ── Data generation ───────────────────────────────────────────────────────────
# K Gaussian blobs with separation controlled by `spread`. Large spread → easy.

make_blobs <- function(n, p, K, spread = 4.0, seed = 1L) {
  set.seed(seed)
  centers <- matrix(rnorm(K * p, sd = spread), K, p)
  per_k   <- rep(seq_len(K), ceiling(n / K))[seq_len(n)]
  labels  <- sample(per_k)
  X       <- centers[labels, ] + matrix(rnorm(n * p), n, p)
  list(X = X, labels = labels, centers = centers)
}

# ── Spectral clustering implementation ───────────────────────────────────────
#
# Returns integer vector of 1-indexed cluster assignments (length n).

spectral_cluster <- function(X_mat, K, sigma = NULL, n_iter_km = 30L,
                              backend = "cpu", seed = 42L) {
  n <- nrow(X_mat); p <- ncol(X_mat)

  # Auto-select sigma as median pairwise distance / sqrt(2*log(n+1))
  if (is.null(sigma)) {
    sub  <- X_mat[seq_len(min(n, 200L)), , drop = FALSE]
    dmat <- as.matrix(dist(sub))
    sigma <- median(dmat[lower.tri(dmat)]) / sqrt(2 * log(n + 1))
    sigma <- max(sigma, 1e-3)
  }

  # Step 1: affinity matrix W [n×n] ──────────────────────────────────────────
  # kernel_matrix with preferred_backend returns a GPU-resident adgeMatrix and
  # zeros the diagonal on device — no CPU round-trip needed.
  W_gpu  <- kernel_matrix(X_mat, kernel = "rbf", sigma = sigma,
                          preferred_backend = backend, zero_diag = TRUE)
  if (!inherits(W_gpu, "adgeMatrix")) {
    diag(W_gpu) <- 0
    W_gpu <- adgeMatrix(W_gpu, preferred_backend = backend, precision = "fast")
  }

  # Step 2 & 3: normalized Laplacian L = D^{-1/2} W D^{-1/2} ───────────────
  d      <- rowSums(W_gpu)                      # n-vector, GPU rowSums
  d_inv  <- 1 / sqrt(pmax(d, 1e-12))            # D^{-1/2}, avoid /0
  L_gpu  <- am_sweep(W_gpu, 1L, d_inv, "*")     # D^{-1/2} W
  L_gpu  <- am_sweep(L_gpu, 2L, d_inv, "*")     # D^{-1/2} W D^{-1/2}

  # Step 4: top-K eigenvectors via rsvd ──────────────────────────────────────
  # L_sym is symmetric PSD → leading left singular vectors = eigenvectors.
  sv     <- rsvd(L_gpu, k = K)
  V      <- sv$u                                # n×K, plain matrix after rsvd

  # Step 5: row-normalise V ──────────────────────────────────────────────────
  row_norms <- sqrt(rowSums(V^2))
  V_norm    <- V / pmax(row_norms, 1e-12)       # unit-norm rows

  # Step 6: K-means on V_norm with k-means++ initialisation ─────────────────
  # k-means++ picks centroids proportional to d² to nearest chosen centroid,
  # greatly improving convergence for sparse/small-n cases.
  set.seed(seed)
  chosen    <- sample.int(n, 1L)
  for (j in seq_len(K - 1L)) {
    cents_so_far <- V_norm[chosen, , drop = FALSE]   # j×K_emb
    # Min squared distance from each point to any centroid so far
    D_cross <- V_norm %*% t(cents_so_far)            # n × j
    c_sq    <- rowSums(cents_so_far^2)               # j-vector
    v_sq    <- rowSums(V_norm^2)                     # n-vector
    D2_mat  <- outer(v_sq, c_sq, "+") - 2 * D_cross  # n × j
    min_d2  <- apply(D2_mat, 1L, min)
    min_d2[chosen] <- 0                              # exclude already chosen
    probs   <- min_d2 / max(sum(min_d2), 1e-15)
    chosen  <- c(chosen, sample.int(n, 1L, prob = probs))
  }
  centroids <- V_norm[chosen, , drop = FALSE]
  labels_km <- integer(n)

  for (iter in seq_len(n_iter_km)) {
    # Use pairwise_sqdist_argmin for the assignment step
    new_labels <- pairwise_sqdist_argmin(V_norm, t(centroids))
    if (identical(new_labels, labels_km)) break
    labels_km <- new_labels
    # Update centroids via segment_mean
    cm <- segment_mean(
      adgeMatrix(V_norm, preferred_backend = backend, precision = "fast"),
      labels_km, K
    )
    centroids <- as.matrix(cm)
    # Handle empty clusters: reinit from random points
    empty <- which(is.na(centroids[, 1L]))
    if (length(empty)) {
      centroids[empty, ] <- V_norm[sample.int(n, length(empty)), , drop = FALSE]
    }
  }
  labels_km
}

# ── Accuracy via best-permutation matching ────────────────────────────────────
# For each permutation of predicted labels, find the one maximising overlap
# with true labels. Exact for K ≤ 8; uses greedy for larger K.

best_accuracy <- function(pred, true, K) {
  # Build K×K confusion matrix
  C <- matrix(0L, K, K)
  for (i in seq_len(length(pred))) C[pred[i], true[i]] <- C[pred[i], true[i]] + 1L
  # Greedy maximum-weight matching (linear assignment via max per row)
  used <- logical(K)
  total <- 0L
  for (k in seq_len(K)) {
    row  <- C[k, ]
    row[used] <- -1L
    best <- which.max(row)
    used[best] <- TRUE
    total <- total + C[k, best]
  }
  total / length(pred)
}

# ── Correctness checks ────────────────────────────────────────────────────────

cat("── Correctness checks ──────────────────────────────────────────────────\n")
all_ok <- TRUE

configs <- list(
  list(n = 300L,  p = 10L, K = 3L, spread = 5.0),   # small, easy
  list(n = 600L,  p = 15L, K = 4L, spread = 4.0),   # medium
  list(n = 300L,  p = 2L,  K = 3L, spread = 3.0),   # very low-dim (p=2)
  list(n = 400L,  p = 30L, K = 5L, spread = 4.5),   # higher-dim
  list(n = 400L,  p = 5L,  K = 6L, spread = 6.0),   # many clusters (k-means++ init needed)
  list(n = 500L,  p = 20L, K = 3L, spread = 3.5)    # moderate overlap
)

for (cfg in configs) {
  dat <- make_blobs(cfg$n, cfg$p, cfg$K, cfg$spread)
  for (bk in active_backends) {
    pred  <- spectral_cluster(dat$X, cfg$K, backend = bk)
    acc   <- best_accuracy(pred, dat$labels, cfg$K)
    ok    <- acc >= 0.85
    all_ok <- all_ok && ok
    cat(sprintf("  [%s] n=%d p=%d K=%d spread=%.1f | acc=%.3f %s\n",
                bk, cfg$n, cfg$p, cfg$K, cfg$spread,
                acc, if (ok) "OK" else "FAIL"))
  }
}
cat(if (all_ok) "\nAll checks PASSED.\n" else "\nSome checks FAILED.\n")

# ── Performance benchmark ─────────────────────────────────────────────────────

if (!requireNamespace("bench", quietly = TRUE)) {
  cat("\nbench not available — skipping timing.\n")
  quit(save = "no")
}
library(bench)

bench_spectral <- function(n, p, K, seed = 1L) {
  dat <- make_blobs(n, p, K, spread = 4.0, seed = seed)
  X   <- dat$X; storage.mode(X) <- "double"

  exprs <- c(
    list(base = quote(spectral_cluster(X, K, backend = "cpu"))),
    setNames(lapply(active_backends, function(bk)
      substitute(spectral_cluster(X, K, backend = bk_), list(bk_ = bk))),
      active_backends)
  )
  res <- bench::mark(exprs = exprs, iterations = 3L,
                     check = FALSE, memory = FALSE, time_unit = "ms")
  res$size <- sprintf("%dx%d K=%d", n, p, K)
  res
}

cat("\n── Performance benchmarks ──────────────────────────────────────────────\n")
for (sz in list(
  list(n = 500L,  p = 10L, K = 4L),
  list(n = 1500L, p = 20L, K = 5L),
  list(n = 3000L, p = 30L, K = 6L)
)) {
  cat(sprintf("\n%dx%d  K=%d:\n", sz$n, sz$p, sz$K))
  res <- bench_spectral(sz$n, sz$p, sz$K)
  print(res[, c("expression", "median", "itr/sec")])
}

# ── Step profiling (n=1500, p=20, K=5) ────────────────────────────────────────

cat("\n── Step profiling  (1500×20  K=5) ─────────────────────────────────────\n")
{
  dat <- make_blobs(1500L, 20L, 5L, spread = 4.0, seed = 99L)
  X   <- dat$X; storage.mode(X) <- "double"; K <- 5L
  bk  <- active_backends[length(active_backends)]

  # Pre-compute intermediate results for isolated step timing
  W_plain <- kernel_matrix(X, kernel = "rbf")
  diag(W_plain) <- 0
  W_gpu   <- adgeMatrix(W_plain, preferred_backend = bk, precision = "fast")
  d       <- rowSums(W_gpu)
  d_inv   <- 1 / sqrt(pmax(d, 1e-12))
  L_gpu   <- am_sweep(am_sweep(W_gpu, 1L, d_inv, "*"), 2L, d_inv, "*")
  sv      <- rsvd(L_gpu, k = K)
  V       <- sv$u
  V_n     <- V / pmax(sqrt(rowSums(V^2)), 1e-12)
  V_gpu   <- adgeMatrix(V_n, preferred_backend = bk, precision = "fast")
  cents   <- V_n[sample.int(1500L, K), ]

  steps <- bench::mark(
    # Old path: compute kernel → download → diag<-0 → re-upload
    kernel_dl_upload = { w <- kernel_matrix(X, kernel="rbf");
                         diag(w) <- 0;
                         adgeMatrix(w, preferred_backend=bk, precision="fast") },
    # New path: compute + zero_diag on GPU → one download (for host copy) → resident
    kernel_resident  = kernel_matrix(X, kernel="rbf",
                                     preferred_backend=bk, zero_diag=TRUE),
    # Laplacian sweep using already-warm resident W
    laplacian_sweep  = { dv <- 1/sqrt(pmax(rowSums(W_gpu), 1e-12));
                         am_sweep(am_sweep(W_gpu, 1L, dv, "*"), 2L, dv, "*") },
    rsvd_L          = rsvd(L_gpu, k = K),
    row_normalize   = { rn <- sqrt(rowSums(V^2)); V / pmax(rn, 1e-12) },
    km_assign       = pairwise_sqdist_argmin(V_n, t(cents)),
    km_update       = segment_mean(V_gpu, sample(K, 1500L, replace=TRUE), K),
    iterations = 5L, check = FALSE, memory = FALSE, time_unit = "ms"
  )
  print(steps[, c("expression", "median", "itr/sec")])
}

cat("\nDone.\n")
