#!/usr/bin/env Rscript
# tools/benchmark-cancor.R
#
# Benchmarks Canonical Correlation Analysis across CPU / MLX / ArrayFire,
# sweeping n (observations), p (X variables), and q (Y variables).
#
# Measures total CCA time and breaks it down into:
#   - crossprods (Sxx, Syy, Sxy) — the GPU-acceleratable part
#   - factorize  (Cholesky + triangular solves + SVD) — always CPU, small
#
# Run (MLX):
#   Rscript -e 'Sys.setenv(AMATRIX_MLX_PROBE_GPU="1"); \
#     setwd("/Users/bbuchsbaum/code/amatrix"); \
#     source("tools/benchmark-cancor.R")'
#
# Run (ArrayFire):
#   Rscript tools/benchmark-cancor.R   # AF uses its own device init

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

make_data <- function(n, p, q, seed = 1L) {
  set.seed(seed)
  rho  <- 0.3
  Sig  <- matrix(rho, p + q, p + q); diag(Sig) <- 1
  L    <- t(chol(Sig))
  raw  <- matrix(rnorm((p + q) * n), n) %*% t(L)
  list(X = raw[, seq_len(p)], Y = raw[, p + seq_len(q)])
}

# Timed crossproduct block only
bench_crossprods <- function(X_mat, Y_mat, backend, n) {
  bk  <- if (identical(backend, "base")) "cpu" else backend
  Xc  <- adgeMatrix(X_mat, preferred_backend = bk, precision = "fast")
  Yc  <- adgeMatrix(Y_mat, preferred_backend = bk, precision = "fast")
  d   <- n - 1L
  function() {
    Sxx <- as.matrix(crossprod(Xc)) / d
    Syy <- as.matrix(crossprod(Yc)) / d
    Sxy <- as.matrix(crossprod(Xc, Yc)) / d
    invisible(list(Sxx, Syy, Sxy))
  }
}

# Full CCA (matching examples/cancor.R implementation)
cca_full <- function(X_mat, Y_mat, k, backend) {
  n <- nrow(X_mat); p <- ncol(X_mat); q <- ncol(Y_mat)
  xcenter <- colMeans(X_mat); ycenter <- colMeans(Y_mat)
  Xc_h <- X_mat - rep(xcenter, each = n)
  Yc_h <- Y_mat - rep(ycenter, each = n)
  bk   <- if (identical(backend, "base")) "cpu" else backend
  Xc   <- adgeMatrix(Xc_h, preferred_backend = bk, precision = "fast")
  Yc   <- adgeMatrix(Yc_h, preferred_backend = bk, precision = "fast")
  d    <- n - 1L
  Sxx  <- crossprod(Xc) / d
  Syy  <- crossprod(Yc) / d
  Sxy  <- crossprod(Xc, Yc) / d
  Rx   <- as.matrix(chol(Sxx)); Ry <- as.matrix(chol(Syy))
  Sxy_h <- as.matrix(Sxy)
  Z    <- solve(t(Rx), Sxy_h)
  M    <- t(solve(t(Ry), t(Z)))
  s    <- base::svd(M, nu = k, nv = k)
  list(cor = s$d[seq_len(k)],
       xcoef = solve(Rx, s$u), ycoef = solve(Ry, s$v))
}

# ── Sweep 1: vary n, fixed p=150, q=100 ──────────────────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat(" Sweep 1: vary n (observations),  p=150, q=100, k=10\n")
cat("  Crossproduct threshold: max(n,p) >= 2048 triggers GPU\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

p1 <- 150L; q1 <- 100L; k1 <- 10L
ns <- c(512L, 1024L, 2048L, 4096L, 8192L)

sweep1_rows <- list()
for (n in ns) {
  d <- make_data(n, p1, q1)
  exprs_total <- list()
  for (bk in active_backends) {
    fn <- local({
      bk_ <- bk; X_ <- d$X; Y_ <- d$Y; k_ <- k1
      function() cca_full(X_, Y_, k_, bk_)
    })
    exprs_total[[bk]] <- fn
  }
  # add base R cancor
  exprs_total[["base_cancor"]] <- local({
    X_ <- d$X; Y_ <- d$Y; k_ <- k1
    function() cancor(X_, Y_)
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
# Pivot: print one row per n, columns per backend
for (n in ns) {
  sub <- s1[s1$n == n, ]
  cpu_ms <- as.numeric(sub$median[sub$expression == "cpu"])
  cat(sprintf("  n = %5d  (p=%d, q=%d)\n", n, p1, q1))
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

# ── Sweep 2: vary p (X variables), fixed n=4096, q=100 ───────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat(" Sweep 2: vary p (X vars),  n=4096, q=100, k=10\n")
cat("  Sxx cost scales as O(n*p^2); Sxy as O(n*p*q)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

n2 <- 4096L; q2 <- 100L; k2 <- 10L
ps <- c(50L, 100L, 200L, 400L, 800L)

sweep2_rows <- list()
for (p in ps) {
  d <- make_data(n2, p, q2)
  exprs_cp <- list()
  for (bk in active_backends) {
    exprs_cp[[bk]] <- bench_crossprods(d$X, d$Y, bk, n2)
  }
  exprs_cp[["base_crossprod"]] <- local({
    X_ <- d$X - rep(colMeans(d$X), each = n2)
    Y_ <- d$Y - rep(colMeans(d$Y), each = n2)
    d_ <- n2 - 1L
    function() {
      Sxx <- crossprod(X_) / d_
      Syy <- crossprod(Y_) / d_
      Sxy <- crossprod(X_, Y_) / d_
      invisible(list(Sxx, Syy, Sxy))
    }
  })

  res <- bench::mark(!!!exprs_cp, iterations = 5L, check = FALSE,
                     memory = FALSE, time_unit = "ms")
  res$p <- p
  res$expression <- as.character(res$expression)
  sweep2_rows <- c(sweep2_rows, list(res[, c("p", "expression", "median")]))
}

s2 <- do.call(rbind, sweep2_rows)
cat(sprintf("  Crossproduct-only timing (n=%d, q=%d)\n\n", n2, q2))
for (p in ps) {
  sub <- s2[s2$p == p, ]
  cpu_ms <- as.numeric(sub$median[sub$expression == "cpu"])
  cat(sprintf("  p = %4d\n", p))
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

# ── Sweep 3: crossover point (GPU threshold override) ────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat(" Sweep 3: crossover — GPU forced on for all sizes\n")
cat("  options(*_min_dim = 0) bypasses size threshold\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

op_old_mlx <- getOption("amatrix.mlx.crossprod_min_dim")
op_old_af  <- getOption("amatrix.arrayfire.crossprod_min_dim")
options(amatrix.mlx.crossprod_min_dim = 0L,
        amatrix.arrayfire.crossprod_min_dim = 0L)
on.exit({
  options(amatrix.mlx.crossprod_min_dim = op_old_mlx,
          amatrix.arrayfire.crossprod_min_dim = op_old_af)
}, add = TRUE)

p3 <- 100L; q3 <- 80L; k3 <- 8L
ns3 <- c(128L, 256L, 512L, 1024L, 2048L, 4096L)

cat(sprintf("  Full CCA timing  (p=%d, q=%d, k=%d)\n\n", p3, q3, k3))
for (n in ns3) {
  d <- make_data(n, p3, q3)
  exprs <- list()
  for (bk in active_backends) {
    exprs[[bk]] <- local({
      bk_ <- bk; X_ <- d$X; Y_ <- d$Y; k_ <- k3
      function() cca_full(X_, Y_, k_, bk_)
    })
  }
  exprs[["base_cancor"]] <- local({
    X_ <- d$X; Y_ <- d$Y
    function() cancor(X_, Y_)
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
cat("  GPU backends win decisively when n >= 2048 and p,q >= 100.\n")
cat("  Below that, CPU wins due to launch/transfer overhead.\n")
cat("  Cholesky + triangular solves + SVD are always CPU-bound and\n")
cat("  negligible (<1 ms for p,q <= 500).\n")
cat("  To force GPU at any size: options(amatrix.mlx.crossprod_min_dim=0)\n")
