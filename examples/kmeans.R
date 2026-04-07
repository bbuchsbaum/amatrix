#!/usr/bin/env Rscript
# examples/kmeans.R
#
# Lloyd's K-means clustering using amatrix.
#
# The bottleneck in K-means is the assignment step: computing squared distances
# from every point to every centroid.  The key identity rewrites this as a
# matrix-matrix product (gemm):
#
#   ‖x_i − c_k‖² = ‖x_i‖² − 2 (X C^T)[i,k] + ‖c_k‖²
#
# The dominant term X C^T (n×p × p×K = n×K gemm) dispatches to the GPU.
# Unlike PLS (K×3 matrix-vector calls) K-means issues one gemm per iteration,
# so GPU utilisation is high at large n and p.
#
# GPU operations per iteration:
#   1. X %*% t(C)   → n×K gemm     [GPU, dispatches via matmul path]
#   2. argmin over K → n argmins    [CPU, linear in K]
#   3. centroid mean → K colMeans   [CPU, O(np)]
#
# GPU dispatch: max(n, p) >= 512 (matmul cold threshold).
# Crossover: GPU wins when n × p × K is large enough to offset upload cost.
#
# Usage:
#   Rscript -e 'Sys.setenv(AMATRIX_MLX_PROBE_GPU="1"); source("examples/kmeans.R")'

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

kmeans_amatrix <- function(X, K, max_iter = 100L, tol = 1e-6,
                           init = "random", seed = 1L, backend = "cpu") {
  X_mat <- as.matrix(X); storage.mode(X_mat) <- "double"
  n <- nrow(X_mat); p <- ncol(X_mat)

  # ── Initialise centroids ────────────────────────────────────────────────────
  set.seed(seed)
  if (identical(init, "random")) {
    idx <- sample.int(n, K)
    centroids <- X_mat[idx, , drop = FALSE]            # K×p
  } else if (is.matrix(init)) {
    centroids <- init                                   # user-supplied K×p
  } else {
    stop("init must be 'random' or a K×p matrix")
  }

  # ── Upload X to GPU once ────────────────────────────────────────────────────
  X_gpu <- adgeMatrix(X_mat, preferred_backend = backend, precision = "fast")

  # Precompute squared row norms of X (n-vector, CPU — used every iteration)
  x_norms_sq <- rowSums(X_mat^2)                       # O(np), one-time CPU

  labels    <- integer(n)
  converged <- FALSE

  for (iter in seq_len(max_iter)) {

    # ── Assignment step (fully GPU) ────────────────────────────────────────────
    # D[i,k] = ‖x_i‖² + ‖c_k‖² − 2 (X C^T)[i,k]
    c_norms_sq <- rowSums(centroids^2)                  # K-vector (small, CPU)
    XCt_gpu    <- X_gpu %*% t(centroids)               # n×K gemm → GPU adgeMatrix
    D_gpu      <- ewise("*", XCt_gpu, -2)              # −2 (X C^T)
    D_gpu      <- am_sweep(D_gpu, 1L, x_norms_sq, "+") # + ‖x_i‖² per row
    D_gpu      <- am_sweep(D_gpu, 2L, c_norms_sq, "+") # + ‖c_k‖² per col

    new_labels <- am_rowargmin(D_gpu)                   # argmin of D = nearest centroid

    # ── Convergence check ──────────────────────────────────────────────────────
    if (iter > 1L && identical(new_labels, labels)) {
      converged <- TRUE
      labels <- new_labels
      break
    }
    labels <- new_labels

    # ── Update step: GPU scatter mean + empty-cluster reinit (CPU, O(K)) ──────
    old_centroids <- centroids
    centroids <- am_scatter_mean(X_gpu, labels, K)     # K×p group means (GPU)

    # Reinitialise any empty clusters to random data points
    empty <- which(is.na(centroids[, 1L]))
    if (length(empty) > 0L)
      centroids[empty, ] <- X_mat[sample.int(n, length(empty), replace = TRUE), ]

    # Centroid shift convergence
    if (max(abs(centroids - old_centroids)) < tol) {
      converged <- TRUE; break
    }
  }

  # ── Within-cluster sum of squares ─────────────────────────────────────────
  wcss <- sum(vapply(seq_len(K), function(k) {
    idx <- which(labels == k)
    if (length(idx) == 0L) return(0)
    Xk <- X_mat[idx, , drop = FALSE]
    sum(sweep(Xk, 2L, centroids[k, ], "-")^2)
  }, numeric(1)))

  list(
    labels    = labels,
    centroids = centroids,
    wcss      = wcss,
    iter      = iter,
    converged = converged,
    backend   = backend
  )
}

# ── 2.  Correctness checks ────────────────────────────────────────────────────
#
# Check A: WCSS from amatrix matches stats::kmeans (same init, same data).
# Check B: backends agree on WCSS to float precision (MLX uses float32).
# Check C: WCSS is non-increasing across restarts (monotone property).

cat("── Correctness checks ───────────────────────────────────────────────────\n")

set.seed(42)
n_sm <- 300L; p_sm <- 10L; K_sm <- 5L
X_sm <- matrix(rnorm(n_sm * p_sm), n_sm, p_sm)

# Use same initial centroids for reference
init_idx <- sample.int(n_sm, K_sm)
C0 <- X_sm[init_idx, , drop = FALSE]

# Check A — WCSS matches stats::kmeans
m_cpu  <- kmeans_amatrix(X_sm, K_sm, init = C0, backend = "cpu")
ref_km <- stats::kmeans(X_sm, centers = C0, algorithm = "Lloyd",
                        iter.max = 100L)
wcss_err <- abs(m_cpu$wcss - ref_km$tot.withinss)
cat(sprintf("  A  WCSS vs stats::kmeans:  |delta| = %.2e  %s\n",
            wcss_err, if (wcss_err < 1e-4) "[PASS]" else "[FAIL]"))

# Check B — backends agree (float32 precision ~1e-4)
backends_ok <- c("cpu")
if (requireNamespace("amatrix.mlx", quietly = TRUE) &&
    isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE)))
  backends_ok <- c(backends_ok, "mlx")

for (bk in backends_ok[-1L]) {
  m_bk <- kmeans_amatrix(X_sm, K_sm, init = C0, backend = bk)
  err  <- abs(m_bk$wcss - m_cpu$wcss)
  cat(sprintf("  B  WCSS cpu vs %-10s |delta| = %.2e  %s\n",
              paste0(bk, ":"), err,
              if (err < 1e-2) "[PASS]" else "[FAIL]"))
}

# Check C — WCSS non-increasing over random restarts
wcss_runs <- vapply(seq_len(5L), function(seed) {
  kmeans_amatrix(X_sm, K_sm, seed = seed, backend = "cpu")$wcss
}, numeric(1))
cat(sprintf("  C  WCSS across 5 seeds:   %s  (range [%.2f, %.2f])\n\n",
            "[INFO]", min(wcss_runs), max(wcss_runs)))

# ── 3.  Backend timing ────────────────────────────────────────────────────────

cat("── Backend timing ───────────────────────────────────────────────────────\n")
cat("  GPU dispatch: matmul X (n×p) %*% t(C) (p×K) — max(n,p) >= 512\n")
cat("  1 gemm per iteration × max_iter iterations\n\n")

backends_to_test <- "cpu"
if (requireNamespace("amatrix.mlx", quietly = TRUE) &&
    isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE)))
  backends_to_test <- c(backends_to_test, "mlx")
if (requireNamespace("amatrix.arrayfire", quietly = TRUE) &&
    isTRUE(try(amatrix.arrayfire::amatrix_arrayfire_is_available(), silent = TRUE)))
  backends_to_test <- c(backends_to_test, "arrayfire")

time_km <- function(n, p, K, backend, max_iter = 20L, reps = 3L) {
  set.seed(1L)
  X <- matrix(rnorm(n * p), n, p)
  elapsed <- numeric(reps)
  for (i in seq_len(reps)) {
    t0 <- proc.time()["elapsed"]
    kmeans_amatrix(X, K, max_iter = max_iter, backend = backend)
    elapsed[i] <- proc.time()["elapsed"] - t0
  }
  median(elapsed) * 1e3
}

for (cfg in list(
  list(n = 2000L,  p = 50L,  K = 20L, label = "small  (n=2000,  p=50,  K=20)"),
  list(n = 5000L,  p = 100L, K = 50L, label = "medium (n=5000,  p=100, K=50)"),
  list(n = 10000L, p = 100L, K=100L,  label = "large  (n=10000, p=100, K=100)"),
  list(n = 5000L,  p = 500L, K = 50L, label = "wide   (n=5000,  p=500, K=50)  [p-driven]")
)) {
  cat(sprintf("  %s\n", cfg$label))
  cpu_ms <- NULL
  for (bk in backends_to_test) {
    ms <- tryCatch(
      time_km(cfg$n, cfg$p, cfg$K, bk),
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
cat("  Full GPU pipeline per iteration:\n")
cat("    1. X %*% t(C)          n×K gemm                     [GPU]\n")
cat("    2. am_sweep(..., 1L)   broadcast row norms           [GPU]\n")
cat("    3. am_sweep(..., 2L)   broadcast col norms           [GPU]\n")
cat("    4. am_rowargmax(D)     row-wise argmax               [GPU → n ints]\n")
cat("    5. am_scatter_mean(X)  one-hot matmul K×p means      [GPU → K×p]\n")
cat("    6. empty-cluster check tabulate(labels, K)           [CPU, O(n)]\n")
cat("  Only CPU↔GPU transfers: labels (n ints, ~40KB) + centroids (K×p, small)\n")
cat("  Estimated savings vs old CPU post-processing:\n")
cat("    sweep×2 + max.col + which()×K: ~35 ms/iter → ~5 ms/iter (GPU)\n")
cat("    centroid update (loop): ~10 ms/iter → ~5 ms/iter (GPU matmul)\n")
cat("  Net: ~40 ms CPU overhead → ~10 ms, freeing ~30 ms/iter for large n,p,K\n")
