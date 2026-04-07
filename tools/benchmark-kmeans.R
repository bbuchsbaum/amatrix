#!/usr/bin/env Rscript
# tools/benchmark-kmeans.R
#
# Benchmarks Lloyd's K-means across CPU / MLX / ArrayFire.
#
# Assignment step: D[i,k] = ‖x_i‖² + ‖c_k‖² − 2(XC^T)[i,k]
# The dominant term X C^T (n×p × p×K = n×K) is a true gemm per iteration.
# CPU post-processing (sweep/argmin/update) is O(nK + np) and dilutes speedup.
# GPU wins when p is large enough that the gemm dominates the CPU overhead.
#
# Run:
#   Rscript -e 'Sys.setenv(AMATRIX_MLX_PROBE_GPU="1"); \
#     setwd("/Users/bbuchsbaum/code/amatrix"); \
#     source("tools/benchmark-kmeans.R")'

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

# ── Implementation ────────────────────────────────────────────────────────────

kmeans_amatrix <- function(X_mat, K, max_iter = 20L, backend = "cpu", seed = 1L) {
  n <- nrow(X_mat); p <- ncol(X_mat)
  set.seed(seed)
  centroids  <- X_mat[sample.int(n, K), , drop = FALSE]
  X_gpu      <- adgeMatrix(X_mat, preferred_backend = backend, precision = "fast")
  x_norms_sq <- rowSums(X_mat^2)
  labels     <- integer(n)

  for (iter in seq_len(max_iter)) {
    # Full GPU pipeline: gemm + broadcast sweep + argmin + scatter mean
    c_norms_sq <- rowSums(centroids^2)
    XCt_gpu    <- X_gpu %*% t(centroids)               # n×K gemm
    D_gpu      <- ewise("*", XCt_gpu, -2)              # −2 XCt
    D_gpu      <- am_sweep(D_gpu, 1L, x_norms_sq, "+") # + row norms
    D_gpu      <- am_sweep(D_gpu, 2L, c_norms_sq, "+") # + col norms
    new_labels <- am_rowargmin(D_gpu)                   # argmin of D

    if (identical(new_labels, labels)) break
    labels <- new_labels

    centroids <- am_scatter_mean(X_gpu, labels, K)      # GPU scatter mean
    empty <- which(is.na(centroids[, 1L]))
    if (length(empty) > 0L)
      centroids[empty, ] <- X_mat[sample.int(n, length(empty), replace = TRUE), ]
  }
  invisible(labels)
}

kmeans_base <- function(X_mat, K, max_iter = 20L, seed = 1L) {
  n <- nrow(X_mat); p <- ncol(X_mat)
  set.seed(seed)
  centroids  <- X_mat[sample.int(n, K), , drop = FALSE]
  x_norms_sq <- rowSums(X_mat^2)
  labels     <- integer(n)

  for (iter in seq_len(max_iter)) {
    XCt        <- X_mat %*% t(centroids)                     # plain BLAS gemm
    D          <- sweep((-2) * XCt, 1L, x_norms_sq,          "+")
    D          <- sweep(D,          2L, rowSums(centroids^2), "+")
    new_labels <- max.col(-D, ties.method = "first")

    if (identical(new_labels, labels)) break
    labels <- new_labels

    for (k in seq_len(K)) {
      idx <- which(labels == k)
      centroids[k, ] <- if (length(idx) == 0L) X_mat[sample.int(n, 1L), ]
                        else .colMeans(X_mat[idx, , drop = FALSE], length(idx), p)
    }
  }
  invisible(labels)
}

time_km <- function(n, p, K, backend, max_iter = 20L, reps = 3L) {
  set.seed(42L)
  # Gaussian mixture — more realistic than uniform noise
  K_true <- min(10L, K)
  means  <- matrix(rnorm(K_true * p, sd = 3), K_true, p)
  cl     <- sample.int(K_true, n, replace = TRUE)
  X      <- means[cl, ] + matrix(rnorm(n * p), n, p)
  elapsed <- numeric(reps)
  for (i in seq_len(reps)) {
    t0 <- proc.time()["elapsed"]
    if (identical(backend, "base")) kmeans_base(X, K, max_iter)
    else                            kmeans_amatrix(X, K, max_iter, backend)
    elapsed[i] <- proc.time()["elapsed"] - t0
  }
  median(elapsed) * 1e3
}

print_row <- function(backends, n, p, K, label) {
  cat(sprintf("  %s\n", label))
  cpu_ms <- NULL
  for (bk in c(backends, "base")) {
    ms <- tryCatch(time_km(n, p, K, bk),
                   error = function(e) { cat("    [", bk, "error:", conditionMessage(e), "]\n"); NA_real_ })
    if (is.na(ms)) next
    if (identical(bk, "cpu")) cpu_ms <- ms
    vs <- if (!is.null(cpu_ms) && !identical(bk, "cpu") && !identical(bk, "base"))
            sprintf("  %+.1fx vs cpu", cpu_ms / ms) else ""
    cat(sprintf("    %-14s %7.1f ms%s\n", bk, ms, vs))
  }
  cat("\n")
}

# ── Sweep 1: vary n, fixed p=200, K=50 ───────────────────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat(" Sweep 1: vary n, p=200, K=50, 20 iterations\n")
cat("  Gemm: n×200 × 200×50 — scales as O(n × p × K) per iter\n")
cat("  CPU overhead (sweep+argmin+update) also O(n): dilutes speedup\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

for (n in c(2000L, 5000L, 10000L, 20000L)) {
  gflop <- as.numeric(n) * 200 * 50 * 20 / 1e9
  print_row(active_backends, n, 200L, 50L,
            sprintf("n = %6d  (%.1f Gops over 20 iters)", n, gflop))
}

# ── Sweep 2: vary p, fixed n=10000, K=50 ─────────────────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat(" Sweep 2: vary p (features), n=10000, K=50, 20 iterations\n")
cat("  Larger p → gemm dominates CPU overhead → GPU crossover\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

for (p in c(50L, 100L, 200L, 500L, 1000L)) {
  ratio <- 10000L * p * 50L / (10000L * 50L)  # gemm / argmin work ratio
  print_row(active_backends, 10000L, p, 50L,
            sprintf("p = %5d  (gemm:argmin work ratio ≈ %d:1)", p, ratio))
}

# ── Sweep 3: vary K, fixed n=10000, p=200 ────────────────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat(" Sweep 3: vary K (clusters), n=10000, p=200, 20 iterations\n")
cat("  Both gemm and argmin scale with K — crossover shifts?\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

for (K in c(10L, 25L, 50L, 100L, 200L)) {
  print_row(active_backends, 10000L, 200L, K,
            sprintf("K = %4d  [n×K argmin matrix = %.0fM elements]", K, 10000*K/1e6))
}

cat("══ Summary ═════════════════════════════════════════════════════\n")
cat("  K-means: ONE gemm per iteration (vs PLS's K×3 gemv).\n")
cat("  GPU speedup grows with p (gemm compute-bound vs CPU overhead).\n")
cat("  CPU post-processing (sweep + argmin + centroid update) is O(nK+np).\n")
cat("  Crossover: GPU wins when p × K / (K + p) >> 1.\n")
cat("  vs KRR: KRR sees 7× (pure O(n²p) gemm, no CPU post-processing).\n")
cat("  vs CCA: CCA ~10× (large batch crossproducts, GPU-resident path).\n")
cat("  vs PLS: PLS ~1× (gemv-bound, many small GPU calls).\n")
