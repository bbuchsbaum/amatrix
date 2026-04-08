#!/usr/bin/env Rscript
# tools/benchmark-pagerank.R
#
# Benchmarks PageRank / power iteration on dense and sparse graphs.
#
# Why this is useful:
#   - Unlike CCA/KRR/LDA, PageRank is repeated matrix-vector multiply, not GEMM.
#   - Arithmetic intensity is low, so bandwidth, dispatch, and launch overhead
#     matter more than peak flop rate.
#   - The sparse variant exercises adgCMatrix -> dense-vector dispatch and the
#     ArrayFire sparse SpMM path, which the current algorithm scripts barely touch.
#
# Algorithm:
#   x_{t+1} = alpha * P x_t + (1 - alpha) * v
# where P is column-stochastic and v is the uniform teleportation vector.
#
# Key primitives exercised:
#   dense:  matmul, scalar ewise, vector add
#   sparse: sparse matmul, dense vector updates, repeated host/device sync for
#           convergence checks
#
# Run:
#   Rscript -e 'Sys.setenv(AMATRIX_MLX_PROBE_GPU="1"); \
#     setwd("/Users/bbuchsbaum/code/amatrix"); \
#     source("tools/benchmark-pagerank.R")'
#
# Optional env vars:
#   AMATRIX_PAGERANK_ENABLE_GPU=1
#     Include GPU backends in the performance sweeps. Off by default because
#     long multi-sweep GPU sessions can expose backend/runtime instability.
#   AMATRIX_PAGERANK_SWEEP=dense|sparse_n|sparse_degree|all
#     Run only a subset of the timing sweeps. Default: all.

run_gpu_perf <- identical(Sys.getenv("AMATRIX_PAGERANK_ENABLE_GPU", "0"), "1")
sweep_mode <- Sys.getenv("AMATRIX_PAGERANK_SWEEP", "all")

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(amatrix)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("install.packages(\"Matrix\") required for this benchmark")
  }
  library(Matrix)
  if (run_gpu_perf && requireNamespace("amatrix.mlx", quietly = TRUE)) {
    amatrix.mlx::amatrix_mlx_register()
  }
  if (run_gpu_perf && requireNamespace("amatrix.arrayfire", quietly = TRUE)) {
    amatrix.arrayfire::amatrix_arrayfire_register()
  }
})

# ── Available backends ────────────────────────────────────────────────────────

active_backends <- "cpu"
mlx_ok <- run_gpu_perf &&
  requireNamespace("amatrix.mlx", quietly = TRUE) &&
  isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE))
af_ok <- run_gpu_perf &&
  requireNamespace("amatrix.arrayfire", quietly = TRUE) &&
  isTRUE(try(amatrix.arrayfire::amatrix_arrayfire_is_available(), silent = TRUE))
if (mlx_ok) active_backends <- c(active_backends, "mlx")
if (af_ok)  active_backends <- c(active_backends, "arrayfire")

sparse_backends <- unique(c("cpu", intersect(active_backends, "arrayfire")))
dense_check_backends <- if (run_gpu_perf) active_backends else "cpu"
sparse_check_backends <- if (run_gpu_perf) sparse_backends else "cpu"
dense_perf_backends <- if (run_gpu_perf) active_backends else "cpu"
sparse_perf_backends <- if (run_gpu_perf) sparse_backends else "cpu"

cat("Active dense backends: ", paste(active_backends, collapse = ", "), "\n", sep = "")
cat("Active sparse backends:", paste(sparse_backends, collapse = ", "), "\n\n")
if (!run_gpu_perf && length(active_backends) > 1L) {
  cat("Correctness checks default to CPU only.\n")
  cat("Performance sweeps default to CPU/base only.\n")
  cat("Set AMATRIX_PAGERANK_ENABLE_GPU=1 to include GPU timings.\n\n")
}

# ── Graph generation ──────────────────────────────────────────────────────────

make_dense_graph <- function(n, seed = 1L) {
  set.seed(seed)
  P <- matrix(rexp(n * n), n, n)
  P <- sweep(P, 2L, colSums(P), "/")
  storage.mode(P) <- "double"
  P
}

make_sparse_graph <- function(n, degree = 16L, seed = 1L) {
  set.seed(seed)
  degree <- min(as.integer(degree), n)
  rows <- as.vector(replicate(n, sample.int(n, degree, replace = FALSE)))
  cols <- rep(seq_len(n), each = degree)
  vals <- rexp(n * degree)
  A <- Matrix::sparseMatrix(
    i = rows, j = cols, x = vals,
    dims = c(n, n), giveCsparse = TRUE
  )
  P <- A %*% Matrix::Diagonal(x = 1 / Matrix::colSums(A))
  as(P, "dgCMatrix")
}

# ── PageRank implementations ──────────────────────────────────────────────────

pagerank_base <- function(P, alpha = 0.85, max_iter = 100L,
                          tol = 1e-10, check_every = 5L) {
  n <- nrow(P)
  x <- matrix(1 / n, n, 1L)
  teleport <- matrix((1 - alpha) / n, n, 1L)
  delta <- Inf

  for (iter in seq_len(max_iter)) {
    x_new <- alpha * (P %*% x) + teleport
    if (iter %% check_every == 0L || iter == max_iter) {
      delta <- sum(abs(x_new - x))
      x <- x_new
      if (tol > 0 && delta < tol) break
    } else {
      x <- x_new
    }
  }

  list(score = as.numeric(x), iterations = iter, delta = delta)
}

pagerank_amatrix <- function(P, backend = "cpu", alpha = 0.85,
                             max_iter = 100L, tol = 1e-8, check_every = 5L) {
  n <- nrow(P)
  P_am <- if (inherits(P, "dgCMatrix")) {
    adgCMatrix(P, preferred_backend = backend, precision = "fast")
  } else {
    adgeMatrix(P, preferred_backend = backend, precision = "fast")
  }
  x <- adgeMatrix(matrix(1 / n, n, 1L),
                  preferred_backend = backend, precision = "fast")
  teleport <- adgeMatrix(matrix((1 - alpha) / n, n, 1L),
                         preferred_backend = backend, precision = "fast")
  delta <- Inf

  for (iter in seq_len(max_iter)) {
    x_new <- alpha * (P_am %*% x) + teleport
    if (iter %% check_every == 0L || iter == max_iter) {
      delta <- sum(abs(as.matrix(x_new) - as.matrix(x)))
      x <- x_new
      if (tol > 0 && delta < tol) break
    } else {
      x <- x_new
    }
  }

  list(score = as.numeric(as.matrix(x)), iterations = iter, delta = delta)
}

plan_for_pagerank <- function(P, backend) {
  x0 <- matrix(1 / nrow(P), nrow(P), 1L)
  P_arg <- if (inherits(P, "dgCMatrix")) {
    adgCMatrix(P, preferred_backend = backend, precision = "fast")
  } else {
    adgeMatrix(P, preferred_backend = backend, precision = "fast")
  }
  amatrix_backend_plan(P_arg, "matmul", y = x0)$chosen
}

tol_for_plan <- function(chosen_backend) {
  if (identical(chosen_backend, "cpu")) 1e-10 else 1e-4
}

# ── Correctness checks ────────────────────────────────────────────────────────

check_dense_case <- function(n, backend, seed = 1L) {
  P <- make_dense_graph(n, seed = seed)
  chosen <- plan_for_pagerank(P, backend)
  tol <- tol_for_plan(chosen)
  ref <- pagerank_base(P, tol = 1e-12)
  am <- pagerank_amatrix(P, backend = backend, tol = tol)
  max_err <- max(abs(ref$score - am$score))
  l1_err <- sum(abs(ref$score - am$score))
  mass_err <- abs(sum(am$score) - 1)
  ok <- max_err < tol * 10 && mass_err < tol * 10
  cat(sprintf(
    "  [dense %s->%s] n=%4d | iters=%3d max_err=%.2e l1_err=%.2e mass_err=%.2e %s\n",
    backend, chosen, n, am$iterations, max_err, l1_err, mass_err,
    if (ok) "OK" else "FAIL"
  ))
  ok
}

check_sparse_case <- function(n, degree, backend, seed = 1L) {
  P <- make_sparse_graph(n, degree = degree, seed = seed)
  chosen <- plan_for_pagerank(P, backend)
  tol <- tol_for_plan(chosen)
  ref <- pagerank_base(P, tol = 1e-12)
  am <- pagerank_amatrix(P, backend = backend, tol = tol)
  max_err <- max(abs(ref$score - am$score))
  l1_err <- sum(abs(ref$score - am$score))
  mass_err <- abs(sum(am$score) - 1)
  ok <- max_err < tol * 10 && mass_err < tol * 10
  cat(sprintf(
    "  [sparse %s->%s] n=%4d deg=%2d nnz=%6d | iters=%3d max_err=%.2e l1_err=%.2e mass_err=%.2e %s\n",
    backend, chosen, n, degree, length(P@x), am$iterations, max_err, l1_err,
    mass_err, if (ok) "OK" else "FAIL"
  ))
  ok
}

cat("── Correctness checks ──────────────────────────────────────────────────\n")
all_ok <- TRUE

for (bk in dense_check_backends) {
  for (n in c(64L, 128L, 256L)) {
    all_ok <- all_ok && check_dense_case(n, bk, seed = 100L + n)
  }
}

for (bk in sparse_check_backends) {
  for (cfg in list(
    list(n = 128L, degree = 4L),
    list(n = 256L, degree = 8L),
    list(n = 512L, degree = 16L)
  )) {
    all_ok <- all_ok && check_sparse_case(cfg$n, cfg$degree, bk,
                                          seed = cfg$n + cfg$degree)
  }
}

cleanup_amatrix <- function() {
  invisible(gc(FALSE))
  invisible(try(amatrix_gc(cache = TRUE), silent = TRUE))
  invisible(gc(FALSE))
}

cleanup_amatrix()

cat(if (all_ok) "\nAll checks PASSED.\n" else "\nSome checks FAILED.\n")

# ── Performance benchmark helpers ────────────────────────────────────────────

time_expr_ms <- function(fn, reps = 3L) {
  elapsed <- numeric(reps)
  for (i in seq_len(reps)) {
    cleanup_amatrix()
    t0 <- proc.time()[["elapsed"]]
    invisible(fn())
    elapsed[i] <- (proc.time()[["elapsed"]] - t0) * 1e3
    cleanup_amatrix()
  }
  stats::median(elapsed)
}

time_named_exprs <- function(exprs, reps = 3L) {
  rows <- lapply(names(exprs), function(name) {
    ms <- tryCatch(
      time_expr_ms(exprs[[name]], reps = reps),
      error = function(e) {
        message("  [", name, " error: ", conditionMessage(e), "]")
        NA_real_
      }
    )
    data.frame(expression = name, median_ms = ms, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

format_plan_suffix <- function(name, plans) {
  if (is.null(plans) || is.na(plans[name]) || identical(name, "cpu")) return("")
  sprintf("  [plan:%s]", plans[name])
}

print_bench_block <- function(res, order, plans = NULL) {
  cpu_ms <- as.numeric(res$median_ms[res$expression == "cpu"])
  for (name in order) {
    idx <- which(res$expression == name)
    if (!length(idx)) next
    ms <- as.numeric(res$median_ms[idx[1L]])
    if (is.na(ms)) {
      cat(sprintf("    %-14s %8s%s\n",
                  name, "ERROR", format_plan_suffix(name, plans)))
      next
    }
    vs <- if (!identical(name, "cpu") && !is.na(cpu_ms)) {
      sprintf("  %+.1fx vs cpu", cpu_ms / ms)
    } else {
      ""
    }
    cat(sprintf("    %-14s %8.1f ms%s%s\n",
                name, ms, vs, format_plan_suffix(name, plans)))
  }
  cat("\n")
}

bench_dense_pagerank <- function(n, max_iter = 50L, seed = 1L) {
  P <- make_dense_graph(n, seed = seed)
  exprs <- list()
  for (bk in dense_perf_backends) {
    exprs[[bk]] <- local({
      P_ <- P
      bk_ <- bk
      function() pagerank_amatrix(P_, backend = bk_, max_iter = max_iter, tol = 0)
    })
  }
  exprs[["base_dense"]] <- local({
    P_ <- P
    function() pagerank_base(P_, max_iter = max_iter, tol = 0)
  })
  res <- time_named_exprs(exprs, reps = 3L)
  plans <- setNames(
    vapply(dense_perf_backends, function(bk) plan_for_pagerank(P, bk), character(1)),
    dense_perf_backends
  )
  list(result = res, plans = plans)
}

bench_sparse_pagerank <- function(n, degree = 16L, max_iter = 50L, seed = 1L) {
  P <- make_sparse_graph(n, degree = degree, seed = seed)
  exprs <- list()
  for (bk in sparse_perf_backends) {
    exprs[[bk]] <- local({
      P_ <- P
      bk_ <- bk
      function() pagerank_amatrix(P_, backend = bk_, max_iter = max_iter, tol = 0)
    })
  }
  exprs[["base_sparse"]] <- local({
    P_ <- P
    function() pagerank_base(P_, max_iter = max_iter, tol = 0)
  })
  res <- time_named_exprs(exprs, reps = 3L)
  plans <- setNames(
    vapply(sparse_perf_backends, function(bk) plan_for_pagerank(P, bk), character(1)),
    sparse_perf_backends
  )
  list(result = res, plans = plans, nnz = length(P@x))
}

# ── Dense sweep ───────────────────────────────────────────────────────────────

if (sweep_mode %in% c("all", "dense")) {
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat(" Dense Sweep: vary n, fixed 50 PageRank iterations\n")
  cat("  Repeated dense matvec: low arithmetic intensity vs GEMM-heavy scripts\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  for (n in c(512L, 1024L, 2048L, 4096L)) {
    cat(sprintf("  n = %5d  [matrix bytes ~= %.1f MiB]\n", n, 8 * n * n / 1024^2))
    dense <- bench_dense_pagerank(n, max_iter = 50L, seed = n)
    print_bench_block(dense$result, c(dense_perf_backends, "base_dense"), dense$plans)
    cleanup_amatrix()
  }
}

# ── Sparse sweep 1 ────────────────────────────────────────────────────────────

spmm_min_nnz <- getOption("amatrix.arrayfire.spmm_min_nnz", 10000L)

if (sweep_mode %in% c("all", "sparse_n")) {
  cat("═══════════════════════════════════════════════════════════════\n")
  cat(" Sparse Sweep 1: vary n, fixed degree=16, fixed 50 iterations\n")
  cat(sprintf("  ArrayFire sparse SpMM threshold: nnz >= %d\n", spmm_min_nnz))
  cat("═══════════════════════════════════════════════════════════════\n\n")

  for (n in c(2048L, 4096L, 8192L, 16384L)) {
    sparse <- bench_sparse_pagerank(n, degree = 16L, max_iter = 50L, seed = n)
    gpu_note <- if (sparse$nnz >= spmm_min_nnz) " [sparse GPU eligible]" else ""
    cat(sprintf("  n = %5d  degree = 16  nnz = %7d%s\n",
                n, sparse$nnz, gpu_note))
    print_bench_block(sparse$result, c(sparse_perf_backends, "base_sparse"), sparse$plans)
    cleanup_amatrix()
  }
}

# ── Sparse sweep 2 ────────────────────────────────────────────────────────────

if (sweep_mode %in% c("all", "sparse_degree")) {
  cat("═══════════════════════════════════════════════════════════════\n")
  cat(" Sparse Sweep 2: vary degree, fixed n=8192, fixed 50 iterations\n")
  cat("  More degree -> more work per iteration and less overhead-dominated runtime\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")

  for (degree in c(4L, 8L, 16L, 32L, 64L)) {
    sparse <- bench_sparse_pagerank(8192L, degree = degree,
                                    max_iter = 50L, seed = 500L + degree)
    gpu_note <- if (sparse$nnz >= spmm_min_nnz) " [sparse GPU eligible]" else ""
    cat(sprintf("  n = %5d  degree = %2d  nnz = %7d%s\n",
                8192L, degree, sparse$nnz, gpu_note))
    print_bench_block(sparse$result, c(sparse_perf_backends, "base_sparse"), sparse$plans)
    cleanup_amatrix()
  }
}

cat("══ Summary ═════════════════════════════════════════════════════\n")
cat("  Dense PageRank is matvec-bound, so GPU crossover should be later than GEMM-heavy methods.\n")
cat("  Sparse PageRank stresses adgCMatrix dispatch, repeated vector updates, and sparse GPU routing.\n")
cat("  If sparse ArrayFire is slower than CPU above the nnz threshold, the likely issue is launch/upload overhead or resident-cache misses.\n")
cat("  If dense GPU plans fall back to cpu, matvec size is still below the backend's useful threshold.\n")
