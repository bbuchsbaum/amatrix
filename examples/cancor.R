#!/usr/bin/env Rscript
# examples/cancor.R
#
# Canonical Correlation Analysis using amatrix.
#
# Demonstrates how GPU-accelerated crossproducts drop into a standard CCA
# algorithm with no change to the calling code.  The three crossproducts
# (Sxx, Syy, Sxy) dispatch to the active backend; everything else (Cholesky,
# triangular solves, SVD of the small k×k matrix) runs on CPU.
#
# GPU acceleration engages automatically when max(n, p) or max(n, q) >= 2048
# (the default crossprod threshold for both MLX and ArrayFire).  At smaller
# sizes the dispatch falls back to CPU with no penalty.
#
# Usage:
#   Rscript -e 'Sys.setenv(AMATRIX_MLX_PROBE_GPU="1"); source("examples/cancor.R")'
#
# Or interactively after pkgload::load_all("."):
#   Sys.setenv(AMATRIX_MLX_PROBE_GPU = "1")
#   source("examples/cancor.R")

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
#
# cancor_amatrix(X, Y, k, center, backend)
#
# Computes the first k canonical correlation pairs between X (n×p) and Y (n×q).
#
# Algorithm:
#   1. Center X, Y (optional).
#   2. Compute sample covariance matrices Sxx (p×p), Syy (q×q), Sxy (p×q)
#      via crossprod — dispatched to the active backend.
#   3. Cholesky: Rx = chol(Sxx), Ry = chol(Syy)   [CPU, small]
#   4. Form M = Rx^{-T} Sxy Ry^{-1}               [CPU, small]
#   5. SVD(M, nu=k, nv=k)                          [CPU, small]
#   6. Recover canonical loadings A = Rx^{-1} U, B = Ry^{-1} V
#
# Returns a list matching base::cancor():
#   cor     — k canonical correlations (descending)
#   xcoef   — p×k X canonical loadings
#   ycoef   — q×k Y canonical loadings
#   xcenter — column means of X (zeros if center=FALSE)
#   ycenter — column means of Y (zeros if center=FALSE)

cancor_amatrix <- function(X, Y, k = NULL, center = TRUE, backend = "auto") {
  X_mat <- as.matrix(X)
  Y_mat <- as.matrix(Y)
  storage.mode(X_mat) <- "double"
  storage.mode(Y_mat) <- "double"

  n <- nrow(X_mat)
  p <- ncol(X_mat)
  q <- ncol(Y_mat)
  if (nrow(Y_mat) != n) stop("X and Y must have the same number of rows")
  k <- min(if (is.null(k)) min(p, q) else k, p, q)

  # Center (CPU — colMeans is fast)
  xcenter <- if (center) colMeans(X_mat) else rep(0.0, p)
  ycenter <- if (center) colMeans(Y_mat) else rep(0.0, q)
  if (center) {
    X_mat <- X_mat - rep(xcenter, each = n)
    Y_mat <- Y_mat - rep(ycenter, each = n)
  }

  # Wrap as adgeMatrix so the backend can be selected
  bk <- if (identical(backend, "auto")) "cpu" else backend
  Xc <- adgeMatrix(X_mat, preferred_backend = bk, precision = "fast")
  Yc <- adgeMatrix(Y_mat, preferred_backend = bk, precision = "fast")

  # ── GPU-accelerated crossproducts + division + Cholesky ───────────────────
  # crossprod() and / both dispatch to the active backend.
  # chol() dispatches through amatrix (GPU when p >= chol threshold, else CPU).
  # All three Gram matrices stay on device until as.matrix() is needed below.
  denom <- n - 1L
  Sxx <- crossprod(Xc) / denom        # p×p  adgeMatrix — stays on device
  Syy <- crossprod(Yc) / denom        # q×q  adgeMatrix — stays on device
  Sxy <- crossprod(Xc, Yc) / denom   # p×q  adgeMatrix — stays on device
  Rx  <- as.matrix(chol(Sxx))         # p×p  upper triangular; downloads here
  Ry  <- as.matrix(chol(Syy))         # q×q  upper triangular; downloads here

  # ── CPU: triangular solves + SVD (no GPU path in either backend) ──────────
  Sxy_h <- as.matrix(Sxy)             # p×q  download for CPU solve
  Z  <- solve(t(Rx), Sxy_h)           # Rx^{-T} Sxy               p×q
  M  <- t(solve(t(Ry), t(Z)))         # Z Ry^{-1}                 p×q
  s  <- base::svd(M, nu = k, nv = k)

  list(
    cor     = s$d[seq_len(k)],
    xcoef   = solve(Rx, s$u),          # p×k
    ycoef   = solve(Ry, s$v),          # q×k
    xcenter = xcenter,
    ycenter = ycenter
  )
}

# ── 2.  Correctness check against base::cancor() ─────────────────────────────

cat("── Correctness check ────────────────────────────────────────────────\n")
set.seed(42)
n <- 200L; p <- 12L; q <- 8L; k <- 4L

# Correlated synthetic data
Sigma <- matrix(0.4, p + q, p + q); diag(Sigma) <- 1
L     <- t(chol(Sigma))
raw   <- matrix(rnorm((p + q) * n), n) %*% t(L)
X_ex  <- raw[, seq_len(p)]
Y_ex  <- raw[, p + seq_len(q)]

ref  <- cancor(X_ex, Y_ex)
am   <- cancor_amatrix(X_ex, Y_ex, k = k, backend = "cpu")

cat(sprintf("  First %d canonical correlations:\n", k))
cat(sprintf("    base::cancor : %s\n",
            paste(round(ref$cor[seq_len(k)], 6), collapse = "  ")))
cat(sprintf("    cancor_amatrix: %s\n",
            paste(round(am$cor, 6), collapse = "  ")))
max_err <- max(abs(am$cor - ref$cor[seq_len(k)]))
cat(sprintf("  Max |delta| in correlations: %.2e  %s\n", max_err,
            if (max_err < 1e-5) "[PASS]" else "[FAIL]"))

# ── 3.  Backend timing comparison ────────────────────────────────────────────

cat("\n── Backend timing (small n=500, p=80, q=60 — CPU expected to win) ──\n")

time_cancor <- function(n, p, q, k, backend, reps = 3L) {
  set.seed(1L)
  X_t <- matrix(rnorm(n * p), n, p)
  Y_t <- matrix(rnorm(n * q), n, q)
  elapsed <- numeric(reps)
  for (i in seq_len(reps)) {
    t0 <- proc.time()["elapsed"]
    cancor_amatrix(X_t, Y_t, k = k, backend = backend)
    elapsed[i] <- proc.time()["elapsed"] - t0
  }
  median(elapsed) * 1e3  # ms
}

backends_to_test <- "cpu"
if (requireNamespace("amatrix.mlx", quietly = TRUE) &&
    isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE)))
  backends_to_test <- c(backends_to_test, "mlx")
if (requireNamespace("amatrix.arrayfire", quietly = TRUE) &&
    isTRUE(try(amatrix.arrayfire::amatrix_arrayfire_is_available(), silent = TRUE)))
  backends_to_test <- c(backends_to_test, "arrayfire")

for (cfg in list(
  list(n = 500L,  p = 80L,  q = 60L,  k = 8L,  label = "small  (n=500,  p=80,  q=60)"),
  list(n = 3000L, p = 200L, q = 150L, k = 20L, label = "medium (n=3000, p=200, q=150)"),
  list(n = 8000L, p = 300L, q = 200L, k = 20L, label = "large  (n=8000, p=300, q=200)")
)) {
  cat(sprintf("\n  %s\n", cfg$label))
  cpu_ms <- NULL
  for (bk in backends_to_test) {
    ms <- tryCatch(
      time_cancor(cfg$n, cfg$p, cfg$q, cfg$k, bk),
      error = function(e) NA_real_
    )
    if (bk == "cpu") cpu_ms <- ms
    if (is.na(ms)) {
      cat(sprintf("    %-12s  [skipped — backend error]\n", paste0(bk, ":")))
      next
    }
    speedup <- if (!is.null(cpu_ms) && bk != "cpu") sprintf("  (%.1fx vs cpu)", cpu_ms / ms) else ""
    cat(sprintf("    %-12s %6.1f ms%s\n", paste0(bk, ":"), ms, speedup))
  }
}

# ── 4.  What the crossover looks like ────────────────────────────────────────

cat("\n── Note on crossover ────────────────────────────────────────────────\n")
cat("  GPU acceleration (MLX / ArrayFire) engages for the crossproduct step\n")
cat("  when max(n, p) >= 2048 (default threshold).\n")
cat("  For n < 2048: all three crossproducts run on CPU — no GPU benefit.\n")
cat("  For n >= 2048 with p,q in the hundreds: GPU wins on the Sxx/Syy/Sxy\n")
cat("  step; Cholesky + SVD remain on CPU and are negligible by comparison.\n")
cat("  Lower the threshold for testing:\n")
cat("    options(amatrix.mlx.crossprod_min_dim = 64L)\n")
cat("    options(amatrix.arrayfire.crossprod_min_dim = 64L)\n")
