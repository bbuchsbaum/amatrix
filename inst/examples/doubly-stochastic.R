#!/usr/bin/env Rscript
# inst/examples/doubly-stochastic.R
#
# Sinkhorn-Knopp algorithm: make a non-negative matrix doubly stochastic
# (all row sums = 1, all col sums = 1).
#
# Algorithm:
#   Given non-negative A [n×n]:
#   Repeat until convergence:
#     1. Row-normalise: A <- A / rowSums(A)    [GPU am_sweep]
#     2. Col-normalise: A <- A / colSums(A)    [GPU am_sweep]
#
# Key primitives exercised (iterated):
#   rowSums, colSums, am_sweep — chained GPU reductions + broadcast scaling.
#
# Correctness: row sums ≈ 1, col sums ≈ 1, non-negativity preserved.
#
# Run:
#   Rscript inst/examples/doubly-stochastic.R

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(amatrix)
  }
  if (requireNamespace("amatrix.arrayfire", quietly = TRUE))
    amatrix.arrayfire::amatrix_arrayfire_register()
  if (requireNamespace("amatrix.mlx", quietly = TRUE))
    amatrix.mlx::amatrix_mlx_register()
})

active_backends <- "cpu"
mlx_ok <- requireNamespace("amatrix.mlx", quietly = TRUE) &&
  isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE))
af_ok  <- requireNamespace("amatrix.arrayfire", quietly = TRUE) &&
  isTRUE(try(amatrix.arrayfire::amatrix_arrayfire_is_available(), silent = TRUE))
if (mlx_ok) active_backends <- c(active_backends, "mlx")
if (af_ok)  active_backends <- c(active_backends, "arrayfire")

cat("Active backends:", paste(active_backends, collapse = ", "), "\n\n")

# ── Sinkhorn-Knopp (amatrix GPU path) ────────────────────────────────────────

sinkhorn_amatrix <- function(A_mat, max_iter = 200L, tol = 1e-8,
                              backend = "cpu") {
  n <- nrow(A_mat)
  A <- adgeMatrix(A_mat, preferred_backend = backend, precision = "fast")

  for (iter in seq_len(max_iter)) {
    # Row normalise
    rs <- rowSums(A)
    A  <- am_sweep(A, 1L, pmax(rs, 1e-15), "/")
    # Col normalise
    cs <- colSums(A)
    A  <- am_sweep(A, 2L, pmax(cs, 1e-15), "/")
    # Convergence: check every 5 iters to reduce GPU→CPU sync overhead
    if (iter %% 5L == 0L || iter == max_iter) {
      rs2 <- rowSums(A)
      cs2 <- colSums(A)
      err <- max(max(abs(rs2 - 1)), max(abs(cs2 - 1)))
      if (err < tol) break
    }
  }
  list(matrix = as.matrix(A), iterations = iter,
       row_err = max(abs(rowSums(as.matrix(A)) - 1)),
       col_err = max(abs(colSums(as.matrix(A)) - 1)))
}

# ── Sinkhorn-Knopp (base R reference) ────────────────────────────────────────

sinkhorn_base <- function(A_mat, max_iter = 200L, tol = 1e-8) {
  A <- A_mat
  for (iter in seq_len(max_iter)) {
    rs <- rowSums(A)
    A  <- A / pmax(rs, 1e-15)
    cs <- colSums(A)
    A  <- t(t(A) / pmax(cs, 1e-15))
    if (iter %% 5L == 0L || iter == max_iter) {
      err <- max(max(abs(rowSums(A) - 1)), max(abs(colSums(A) - 1)))
      if (err < tol) break
    }
  }
  list(matrix = A, iterations = iter,
       row_err = max(abs(rowSums(A) - 1)),
       col_err = max(abs(colSums(A) - 1)))
}

# ── Data generation ──────────────────────────────────────────────────────────

make_nonneg <- function(n, seed = 42L) {
  set.seed(seed)
  # Positive-entry matrix (exp of Gaussian ensures strict positivity)
  A <- exp(matrix(rnorm(n * n), n, n))
  storage.mode(A) <- "double"
  A
}

# ── Correctness checks ──────────────────────────────────────────────────────

cat("── Correctness checks ──────────────────────────────────────────────────\n")
all_ok <- TRUE

configs <- list(
  list(n = 10L,   tol = 1e-10),   # tiny — full convergence
  list(n = 50L,   tol = 1e-8),    # small
  list(n = 200L,  tol = 1e-6),    # medium
  list(n = 500L,  tol = 1e-6),    # larger
  list(n = 100L,  tol = 1e-8)     # moderate, tight tol
)

for (cfg in configs) {
  A <- make_nonneg(cfg$n)
  ref <- sinkhorn_base(A, tol = cfg$tol)

  for (bk in active_backends) {
    am <- sinkhorn_amatrix(A, tol = cfg$tol, backend = bk)

    # Check doubly stochastic property
    rs_ok <- am$row_err < cfg$tol * 10
    cs_ok <- am$col_err < cfg$tol * 10
    # Check non-negativity
    nn_ok <- all(am$matrix >= -1e-12)
    # Check agreement with base R (after same iterations)
    diff  <- max(abs(am$matrix - ref$matrix))
    # GPU float32 accumulation → 1e-4 tolerance for cross-backend agreement
    agree <- diff < 1e-3

    ok <- rs_ok && cs_ok && nn_ok && agree
    all_ok <- all_ok && ok

    cat(sprintf("  [%s] n=%3d tol=%.0e | iters=%d row_err=%.1e col_err=%.1e diff=%.1e %s\n",
                bk, cfg$n, cfg$tol, am$iterations,
                am$row_err, am$col_err, diff,
                if (ok) "OK" else "FAIL"))
  }
}
cat(if (all_ok) "\nAll checks PASSED.\n" else "\nSome checks FAILED.\n")

# ── Performance benchmark ────────────────────────────────────────────────────

if (!requireNamespace("bench", quietly = TRUE)) {
  cat("\nbench not available — skipping timing.\n")
  quit(save = "no")
}
library(bench)

# ── Sinkhorn via resident_handle (zero S4 overhead per step) ─────────────────

sinkhorn_handle <- function(A_mat, max_iter = 200L, tol = 1e-8, backend = "cpu") {
  h <- resident_handle(adgeMatrix(A_mat, preferred_backend = backend, precision = "fast"))
  for (iter in seq_len(max_iter)) {
    rs <- rh_rowSums(h)
    am_sweep_inplace(h, 1L, pmax(rs, 1e-15), "/")
    cs <- rh_colSums(h)
    am_sweep_inplace(h, 2L, pmax(cs, 1e-15), "/")
    if (iter %% 5L == 0L || iter == max_iter) {
      rs2 <- rh_rowSums(h)
      cs2 <- rh_colSums(h)
      err <- max(max(abs(rs2 - 1)), max(abs(cs2 - 1)))
      if (err < tol) break
    }
  }
  mat <- as.matrix(h)
  list(matrix = mat, iterations = iter,
       row_err = max(abs(rowSums(mat) - 1)),
       col_err = max(abs(colSums(mat) - 1)))
}

bench_sinkhorn <- function(n, n_iter_fixed = 50L, seed = 1L) {
  A <- make_nonneg(n, seed = seed)

  # Handle path only works on GPU backends with residency support
  gpu_backends <- setdiff(active_backends, "cpu")
  exprs <- c(
    list(base = substitute(sinkhorn_base(A, max_iter = NI, tol = 0),
                           list(NI = n_iter_fixed))),
    setNames(lapply(gpu_backends, function(bk)
      substitute(sinkhorn_handle(A, max_iter = NI, tol = 0, backend = bk_),
                 list(NI = n_iter_fixed, bk_ = bk))),
      paste0(gpu_backends, "_handle"))
  )
  res <- bench::mark(exprs = exprs, iterations = 3L,
                     check = FALSE, memory = FALSE, time_unit = "ms")
  res$size <- sprintf("%dx%d", n, n)
  res
}

cat("\n── Performance benchmarks (fixed 50 Sinkhorn iterations) ────────────\n")
for (sz in list(
  list(n = 200L),
  list(n = 500L),
  list(n = 1000L),
  list(n = 2000L)
)) {
  cat(sprintf("\n%dx%d:\n", sz$n, sz$n))
  res <- bench_sinkhorn(sz$n)
  print(res[, c("expression", "median", "itr/sec")])
}

# ── Step profiling ───────────────────────────────────────────────────────────

cat("\n── Step profiling  (1000×1000, single iteration) ──────────────────────\n")
{
  A <- make_nonneg(1000L, seed = 99L)
  bk <- active_backends[length(active_backends)]
  A_gpu <- adgeMatrix(A, preferred_backend = bk, precision = "fast")
  rs <- rowSums(A_gpu)
  cs <- colSums(A_gpu)

  steps <- bench::mark(
    rowSums       = rowSums(A_gpu),
    row_sweep     = am_sweep(A_gpu, 1L, pmax(rs, 1e-15), "/"),
    colSums       = colSums(A_gpu),
    col_sweep     = am_sweep(A_gpu, 2L, pmax(cs, 1e-15), "/"),
    full_iter     = { Ag <- am_sweep(A_gpu, 1L, pmax(rowSums(A_gpu), 1e-15), "/");
                      am_sweep(Ag, 2L, pmax(colSums(Ag), 1e-15), "/") },
    iterations = 10L, check = FALSE, memory = FALSE, time_unit = "ms"
  )
  print(steps[, c("expression", "median", "itr/sec")])
}

cat("\nDone.\n")
