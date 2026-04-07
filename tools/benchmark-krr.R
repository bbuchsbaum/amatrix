#!/usr/bin/env Rscript
# tools/benchmark-krr.R
#
# Benchmarks Kernel Ridge Regression across CPU / MLX / ArrayFire,
# sweeping n (observations) and p (features).
#
# KRR kernel matrix K = X X^T / p uses tcrossprod — the GPU path CCA
# does not stress.  Crossover threshold: max(n, p) >= 2048.
#
# Run (MLX):
#   Rscript -e 'Sys.setenv(AMATRIX_MLX_PROBE_GPU="1"); \
#     setwd("/Users/bbuchsbaum/code/amatrix"); \
#     source("tools/benchmark-krr.R")'

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(amatrix)
  }
  if (!requireNamespace("bench", quietly = TRUE))
    stop("install.packages('bench') required for this benchmark")
  library(bench)
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

# ── Helpers ───────────────────────────────────────────────────────────────────

make_krr_data <- function(n, p, seed = 1L) {
  set.seed(seed)
  beta <- rnorm(p)
  X    <- matrix(rnorm(n * p), n, p)
  y    <- X %*% beta + rnorm(n, sd = 0.5)
  list(X = X, y = y)
}

# GPU-timed kernel formation only (tcrossprod + scalar div + download)
bench_tcrossprod <- function(X_mat, backend, n) {
  bk  <- if (identical(backend, "base")) "cpu" else backend
  Xc  <- adgeMatrix(X_mat, preferred_backend = bk, precision = "fast")
  p   <- ncol(X_mat)
  function() {
    K <- as.matrix(tcrossprod(Xc)) / p
    invisible(K)
  }
}

# Full KRR: form kernel + Cholesky solve + predict (K_test = X X^T cross-case)
krr_full <- function(X_mat, y_mat, lambda, backend) {
  n <- nrow(X_mat); p <- ncol(X_mat)
  bk <- if (identical(backend, "base")) "cpu" else backend
  Xc <- adgeMatrix(X_mat, preferred_backend = bk, precision = "fast")

  # GPU tcrossprod then GPU chol (both dispatch through amatrix)
  K   <- tcrossprod(Xc) / p
  K_h <- as.matrix(K)
  diag(K_h) <- diag(K_h) + lambda
  Kc    <- adgeMatrix(K_h, preferred_backend = bk, precision = "fast")
  R     <- as.matrix(chol(Kc))   # GPU chol when n >= 256
  alpha <- backsolve(R, forwardsolve(t(R), y_mat))

  # Prediction on training set: cross-case tcrossprod(X_test, X_train)
  # Here X_test == X_train so n_test == n_train, but uses the cross-case bridge
  Xtrn <- adgeMatrix(X_mat, preferred_backend = bk, precision = "fast")
  K_test <- tcrossprod(Xc, Xtrn) / p
  as.matrix(K_test) %*% alpha
}

# Base R reference (no amatrix)
krr_base <- function(X_mat, y_mat, lambda) {
  p  <- ncol(X_mat)
  K  <- tcrossprod(X_mat) / p
  diag(K) <- diag(K) + lambda
  R     <- chol(K)
  alpha <- backsolve(R, forwardsolve(t(R), y_mat))
  K %*% alpha
}

# ── Sweep 1: vary n, fixed p=100 ─────────────────────────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat(" Sweep 1: vary n (observations),  p=100, lambda=0.5\n")
cat("  tcrossprod threshold: max(n,p) >= 2048 triggers GPU\n")
cat("  Work scales as O(n^2 p) — very different from CCA's O(n p^2)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

p1 <- 100L; lam1 <- 0.5
ns <- c(512L, 1024L, 2048L, 4096L)

sweep1_rows <- list()
for (n in ns) {
  d <- make_krr_data(n, p1)
  exprs_total <- list()
  for (bk in active_backends) {
    fn <- local({
      bk_ <- bk; X_ <- d$X; y_ <- d$y; lam_ <- lam1
      function() krr_full(X_, y_, lam_, bk_)
    })
    exprs_total[[bk]] <- fn
  }
  exprs_total[["base_krr"]] <- local({
    X_ <- d$X; y_ <- d$y; lam_ <- lam1
    function() krr_base(X_, y_, lam_)
  })

  res <- tryCatch(
    bench::mark(!!!exprs_total, iterations = 5L, check = FALSE,
                memory = FALSE, time_unit = "ms"),
    error = function(e) {
      message("  [bench error for n=", n, ": ", conditionMessage(e), "]")
      NULL
    }
  )
  if (is.null(res)) next
  res$n <- n
  res$expression <- as.character(res$expression)
  sweep1_rows <- c(sweep1_rows, list(res[, c("n", "expression", "median", "min", "max")]))
}

s1 <- do.call(rbind, sweep1_rows)
for (n in ns) {
  sub <- s1[s1$n == n, ]
  cpu_ms <- as.numeric(sub$median[sub$expression == "cpu"])
  cat(sprintf("  n = %5d  (p=%d)\n", n, p1))
  for (i in seq_len(nrow(sub))) {
    bk    <- sub$expression[i]
    ms    <- as.numeric(sub$median[i])
    vs    <- if (bk != "cpu" && !is.na(cpu_ms))
               sprintf("  %+.1fx vs cpu", cpu_ms / ms)
             else ""
    above <- if (n >= 2048L && bk %in% c("mlx", "arrayfire")) " [GPU on]" else ""
    cat(sprintf("    %-14s %7.1f ms%s%s\n", bk, ms, vs, above))
  }
  cat("\n")
}

# ── Sweep 2: vary p (features), fixed n=4096 ─────────────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat(" Sweep 2: tcrossprod-only timing — vary p, fixed n=4096\n")
cat("  More p = more work per kernel element = better GPU utilisation\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

n2 <- 4096L; lam2 <- 0.5
ps <- c(50L, 100L, 200L, 400L, 800L)

sweep2_rows <- list()
for (p in ps) {
  d <- make_krr_data(n2, p)
  exprs_cp <- list()
  for (bk in active_backends) {
    exprs_cp[[bk]] <- bench_tcrossprod(d$X, bk, n2)
  }
  exprs_cp[["base_tcrossprod"]] <- local({
    X_ <- d$X; p_ <- p
    function() invisible(tcrossprod(X_) / p_)
  })

  res <- bench::mark(!!!exprs_cp, iterations = 5L, check = FALSE,
                     memory = FALSE, time_unit = "ms")
  res$p <- p
  res$expression <- as.character(res$expression)
  sweep2_rows <- c(sweep2_rows, list(res[, c("p", "expression", "median")]))
}

s2 <- do.call(rbind, sweep2_rows)
cat(sprintf("  tcrossprod-only timing (n=%d)\n\n", n2))
for (p in ps) {
  sub <- s2[s2$p == p, ]
  cpu_ms <- as.numeric(sub$median[sub$expression == "cpu"])
  cat(sprintf("  p = %4d  [kernel cost O(n^2 p) = %.0fM ops]\n",
              p, n2^2 * p / 1e6))
  for (i in seq_len(nrow(sub))) {
    bk <- sub$expression[i]
    ms <- as.numeric(sub$median[i])
    vs <- if (bk != "cpu" && !is.na(cpu_ms))
            sprintf("  %+.1fx vs cpu", cpu_ms / ms)
          else ""
    cat(sprintf("    %-18s %7.1f ms%s\n", bk, ms, vs))
  }
  cat("\n")
}

# ── Sweep 3: crossover (GPU forced on for all sizes) ──────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat(" Sweep 3: crossover — GPU forced on (tcrossprod_min_dim = 0)\n")
cat("  Shows GPU overhead at small n vs break-even point\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

op_old_mlx <- getOption("amatrix.mlx.tcrossprod_min_dim")
op_old_af  <- getOption("amatrix.arrayfire.tcrossprod_min_dim")
options(amatrix.mlx.tcrossprod_min_dim = 0L,
        amatrix.arrayfire.tcrossprod_min_dim = 0L)
on.exit({
  options(amatrix.mlx.tcrossprod_min_dim = op_old_mlx,
          amatrix.arrayfire.tcrossprod_min_dim = op_old_af)
}, add = TRUE)

p3 <- 100L; lam3 <- 0.5
ns3 <- c(128L, 256L, 512L, 1024L, 2048L, 4096L)

cat(sprintf("  Full KRR timing  (p=%d, lambda=%.1f)\n\n", p3, lam3))
for (n in ns3) {
  d <- make_krr_data(n, p3)
  exprs <- list()
  for (bk in active_backends) {
    exprs[[bk]] <- local({
      bk_ <- bk; X_ <- d$X; y_ <- d$y; lam_ <- lam3
      function() krr_full(X_, y_, lam_, bk_)
    })
  }
  exprs[["base_krr"]] <- local({
    X_ <- d$X; y_ <- d$y; lam_ <- lam3
    function() krr_base(X_, y_, lam_)
  })

  res <- bench::mark(!!!exprs, iterations = 5L, check = FALSE,
                     memory = FALSE, time_unit = "ms")
  res$expression <- as.character(res$expression)
  cpu_ms <- as.numeric(res$median[res$expression == "cpu"])
  cat(sprintf("  n = %5d\n", n))
  for (i in seq_len(nrow(res))) {
    bk <- res$expression[i]
    ms <- as.numeric(res$median[i])
    vs <- if (bk != "cpu" && !is.na(cpu_ms))
            sprintf("  %+.1fx vs cpu", cpu_ms / ms)
          else ""
    cat(sprintf("    %-14s %7.1f ms%s\n", bk, ms, vs))
  }
  cat("\n")
}

cat("══ Summary ═════════════════════════════════════════════════════\n")
cat("  KRR is n-driven: GPU wins when n >= 2048 (kernel O(n^2 p)).\n")
cat("  Higher p improves GPU efficiency at the same n.\n")
cat("  Cholesky solve on n×n (CPU for n < 256, GPU above) is secondary.\n")
cat("  Prediction uses cross-case tcrossprod(X_test, X_train) —\n")
cat("  the same non-square path fixed for CCA cross-covariance.\n")
cat("  To force GPU at any n: options(amatrix.mlx.tcrossprod_min_dim=0)\n")
