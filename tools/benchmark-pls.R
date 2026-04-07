#!/usr/bin/env Rscript
# tools/benchmark-pls.R
#
# Benchmarks NIPALS PLS1 across CPU / MLX / ArrayFire.
#
# PLS is sequential (K components, each depending on the previous deflation)
# and p-driven: each component costs 3 × O(np) GPU matmuls (gemv).
# Unlike CCA/KRR, PLS issues many small GPU kernel launches, so GPU wins
# only at large p (>> 2000) where gemv parallelism compensates.
#
# Implicit-deflation variant: X uploaded to GPU once; accumulated T,P
# corrections applied on CPU each iteration (O(K·(n+p)) total overhead).
#
# Sweep 1: vary n, fixed p=500, K=10  — n-scaling at moderate p
# Sweep 2: vary p, fixed n=4096, K=10 — p-scaling (GPU crossover finder)
# Sweep 3: vary K, fixed n=4096, p=500 — component-count scaling
#
# Run:
#   Rscript -e 'Sys.setenv(AMATRIX_MLX_PROBE_GPU="1"); \
#     setwd("/Users/bbuchsbaum/code/amatrix"); \
#     source("tools/benchmark-pls.R")'

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

# ── NIPALS PLS1 implementation (implicit deflation) ───────────────────────────

pls_amatrix <- function(X_mat, y_vec, K = 10L, backend = "cpu") {
  n <- nrow(X_mat); p <- ncol(X_mat)
  K <- min(K, p, n - 1L)
  W <- matrix(0, p, K); P <- matrix(0, p, K)
  T <- matrix(0, n, K); b <- numeric(K)

  X_gpu <- adgeMatrix(X_mat, preferred_backend = backend, precision = "fast")
  r_cur <- y_vec

  for (k in seq_len(K)) {
    Tk <- if (k > 1L) T[, seq_len(k - 1L), drop = FALSE] else NULL
    Pk <- if (k > 1L) P[, seq_len(k - 1L), drop = FALSE] else NULL

    xt_r <- as.vector(crossprod(X_gpu, r_cur))
    if (k > 1L) xt_r <- xt_r - as.vector(Pk %*% (t(Tk) %*% r_cur))
    w_k <- xt_r / sqrt(sum(xt_r^2))

    t_vec <- as.vector(X_gpu %*% w_k)
    if (k > 1L) t_vec <- t_vec - as.vector(Tk %*% (t(Pk) %*% w_k))
    tt <- sum(t_vec^2)

    xt_t <- as.vector(crossprod(X_gpu, t_vec))
    if (k > 1L) xt_t <- xt_t - as.vector(Pk %*% (t(Tk) %*% t_vec))
    p_k <- xt_t / tt

    b[k] <- sum(t_vec * r_cur) / tt
    r_cur <- r_cur - b[k] * t_vec
    W[, k] <- w_k; P[, k] <- p_k; T[, k] <- t_vec
  }

  PtW    <- crossprod(P, W)
  W_star <- W %*% solve(PtW)
  as.vector(W_star %*% b)
}

# Base R reference (no amatrix, explicit deflation)
pls_base <- function(X_mat, y_vec, K = 10L) {
  n <- nrow(X_mat); p <- ncol(X_mat); K <- min(K, p, n - 1L)
  W <- matrix(0, p, K); P <- matrix(0, p, K); b <- numeric(K)
  X_cur <- X_mat; r_cur <- y_vec

  for (k in seq_len(K)) {
    w_k  <- crossprod(X_cur, r_cur)
    w_k  <- w_k / sqrt(sum(w_k^2))
    t_vec <- X_cur %*% w_k
    tt    <- sum(t_vec^2)
    p_k   <- crossprod(X_cur, t_vec) / tt
    b[k]  <- sum(t_vec * r_cur) / tt
    r_cur <- r_cur - b[k] * t_vec
    X_cur <- X_cur - outer(t_vec, p_k)
    W[, k] <- w_k; P[, k] <- p_k
  }
  as.vector(W %*% solve(crossprod(P, W), b))
}

make_data <- function(n, p, seed = 1L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  y <- X %*% rnorm(p) + rnorm(n, sd = 0.5)
  list(X = X, y = y)
}

# ── Sweep 1: vary n, fixed p=500, K=10 ───────────────────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat(" Sweep 1: vary n,  p=500, K=10\n")
cat("  Each component: 3 × O(np) GPU gemv + O(K·(n+p)) CPU correction\n")
cat("  GPU crossover: max(n,p) >= 2048  (crossprod cold threshold)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

p1 <- 500L; K1 <- 10L
ns <- c(1024L, 2048L, 4096L)

for (n in ns) {
  d <- make_data(n, p1)
  exprs <- list()
  for (bk in active_backends) {
    exprs[[bk]] <- local({
      bk_ <- bk; X_ <- d$X; y_ <- d$y; K_ <- K1
      function() pls_amatrix(X_, y_, K_, bk_)
    })
  }
  exprs[["base_pls"]] <- local({
    X_ <- d$X; y_ <- d$y; K_ <- K1
    function() pls_base(X_, y_, K_)
  })

  res <- tryCatch(
    bench::mark(!!!exprs, iterations = 5L, check = FALSE,
                memory = FALSE, time_unit = "ms"),
    error = function(e) {
      message("  [bench error n=", n, ": ", conditionMessage(e), "]"); NULL
    }
  )
  if (is.null(res)) next
  res$n <- n; res$expression <- as.character(res$expression)
  cpu_ms <- as.numeric(res$median[res$expression == "cpu"])

  cat(sprintf("  n = %5d  (p=%d, K=%d)\n", n, p1, K1))
  for (i in seq_len(nrow(res))) {
    bk <- res$expression[i]; ms <- as.numeric(res$median[i])
    vs <- if (bk != "cpu" && !is.na(cpu_ms))
            sprintf("  %+.1fx vs cpu", cpu_ms / ms) else ""
    cat(sprintf("    %-14s %7.1f ms%s\n", bk, ms, vs))
  }
  cat("\n")
}

# ── Sweep 2: vary p, fixed n=4096, K=10 ──────────────────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat(" Sweep 2: vary p (features), n=4096, K=10\n")
cat("  GPU wins when p is large enough (gemv becomes compute-bound)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

n2 <- 4096L; K2 <- 10L
ps <- c(100L, 500L, 1000L, 2000L)

for (p in ps) {
  d <- make_data(n2, p)
  exprs <- list()
  for (bk in active_backends) {
    exprs[[bk]] <- local({
      bk_ <- bk; X_ <- d$X; y_ <- d$y; K_ <- K2
      function() pls_amatrix(X_, y_, K_, bk_)
    })
  }
  exprs[["base_pls"]] <- local({
    X_ <- d$X; y_ <- d$y; K_ <- K2
    function() pls_base(X_, y_, K_)
  })

  res <- tryCatch(
    bench::mark(!!!exprs, iterations = 5L, check = FALSE,
                memory = FALSE, time_unit = "ms"),
    error = function(e) {
      message("  [bench error p=", p, ": ", conditionMessage(e), "]"); NULL
    }
  )
  if (is.null(res)) next
  res$p <- p; res$expression <- as.character(res$expression)
  cpu_ms <- as.numeric(res$median[res$expression == "cpu"])

  cat(sprintf("  p = %5d  (n=%d, K=%d)  [%dM gemv ops/iter]\n",
              p, n2, K2, as.integer(n2 * p / 1e6)))
  for (i in seq_len(nrow(res))) {
    bk <- res$expression[i]; ms <- as.numeric(res$median[i])
    vs <- if (bk != "cpu" && !is.na(cpu_ms))
            sprintf("  %+.1fx vs cpu", cpu_ms / ms) else ""
    cat(sprintf("    %-14s %7.1f ms%s\n", bk, ms, vs))
  }
  cat("\n")
}

# ── Sweep 3: vary K (components), fixed n=4096, p=500 ────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat(" Sweep 3: vary K (components), n=4096, p=500\n")
cat("  Time should scale linearly with K\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

n3 <- 4096L; p3 <- 500L
Ks <- c(5L, 10L, 20L, 40L)
d3 <- make_data(n3, p3)

for (K in Ks) {
  exprs <- list()
  for (bk in active_backends) {
    exprs[[bk]] <- local({
      bk_ <- bk; X_ <- d3$X; y_ <- d3$y; K_ <- K
      function() pls_amatrix(X_, y_, K_, bk_)
    })
  }
  exprs[["base_pls"]] <- local({
    X_ <- d3$X; y_ <- d3$y; K_ <- K
    function() pls_base(X_, y_, K_)
  })

  res <- tryCatch(
    bench::mark(!!!exprs, iterations = 5L, check = FALSE,
                memory = FALSE, time_unit = "ms"),
    error = function(e) {
      message("  [bench error K=", K, ": ", conditionMessage(e), "]"); NULL
    }
  )
  if (is.null(res)) next
  res$K <- K; res$expression <- as.character(res$expression)
  cpu_ms <- as.numeric(res$median[res$expression == "cpu"])

  cat(sprintf("  K = %3d  (n=%d, p=%d)\n", K, n3, p3))
  for (i in seq_len(nrow(res))) {
    bk <- res$expression[i]; ms <- as.numeric(res$median[i])
    vs <- if (bk != "cpu" && !is.na(cpu_ms))
            sprintf("  %+.1fx vs cpu", cpu_ms / ms) else ""
    cat(sprintf("    %-14s %7.1f ms%s\n", bk, ms, vs))
  }
  cat("\n")
}

cat("══ Summary ═════════════════════════════════════════════════════\n")
cat("  PLS issues K×3 sequential gemv calls — GPU latency dominates.\n")
cat("  GPU wins only at large p (>> 2000) where gemv is compute-bound.\n")
cat("  Compare to KRR (n-driven, O(n^2 p) gemm) which sees 7× at n=4096.\n")
cat("  To benchmark at the GPU crossover:\n")
cat("    options(amatrix.mlx.crossprod_min_dim = 512L)\n")
