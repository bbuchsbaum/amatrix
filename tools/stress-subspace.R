#!/usr/bin/env Rscript
# tools/stress-subspace.R
#
# Subspace-iteration stress test.  Exercises two large GEMMs per Krylov step,
# QR reorthogonalization, a small SVD, and reconstruction accuracy.
#
# This script reveals:
#   - Float32 precision loss (rel_err > ~1.5 after q iterations)
#   - Backend dispatch regressions (crash mid-loop or silent CPU fallback)
#   - Speedup headroom on each matrix size
#   - The as.matrix() boundary between device and host
#
# Usage:
#   Rscript tools/stress-subspace.R                      # all backends, all sizes
#   Rscript tools/stress-subspace.R --backend=mlx        # single backend
#   Rscript tools/stress-subspace.R --size=large         # single size
#   Rscript tools/stress-subspace.R --q=5                # deeper power iteration

# ── Boilerplate: repo root + library path ────────────────────────────────────

bench_lib <- Sys.getenv("AMATRIX_BENCH_LIB", "")
if (nzchar(bench_lib)) .libPaths(c(normalizePath(bench_lib), .libPaths()))

script_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
repo_root <- if (length(script_file)) {
  normalizePath(dirname(dirname(sub("^--file=", "", script_file[[1L]]))), mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  stop("Run from the amatrix repo root or pass --file= to Rscript")
}
suppressPackageStartupMessages(
  pkgload::load_all(repo_root, quiet = TRUE)
)

# ── Argument parsing ─────────────────────────────────────────────────────────

args  <- commandArgs(trailingOnly = TRUE)
arg   <- function(key, default) {
  m <- regmatches(args, regexpr(paste0("(?<=--", key, "=)\\S+"), args, perl = TRUE))
  if (length(m)) m[[1L]] else default
}
flag  <- function(key) any(args == paste0("--", key))

SIZE_FILTER    <- arg("size",    "all")
BACKEND_FILTER <- arg("backend", "all")
Q_ITER         <- as.integer(arg("q", "3"))
N_TRIALS       <- as.integer(arg("trials", "5"))
K_RANK         <- as.integer(arg("k", "20"))
REL_ERR_WARN   <- 1.5   # relative reconstruction error: warn above this
REL_ERR_FAIL   <- 3.0   # fail above this

# ── Sizes ────────────────────────────────────────────────────────────────────

ALL_SIZES <- list(
  small  = list(n = 256L,  p =  64L, label = "256×64"),
  medium = list(n = 1024L, p = 128L, label = "1024×128"),
  large  = list(n = 4096L, p = 256L, label = "4096×256")
)
SIZES <- if (SIZE_FILTER == "all") ALL_SIZES else ALL_SIZES[SIZE_FILTER]
if (!length(SIZES)) stop("Unknown --size=", SIZE_FILTER,
                         ".  Choose: small, medium, large")

# ── Backend detection ────────────────────────────────────────────────────────

detect_backends <- function() {
  out <- "cpu"
  mlx_ok <- tryCatch({
    Sys.setenv(AMATRIX_MLX_PROBE_GPU = "1")
    requireNamespace("amatrix.mlx", quietly = TRUE) &&
      amatrix.mlx::amatrix_mlx_is_available()
  }, error = function(e) FALSE)
  if (mlx_ok) out <- c(out, "mlx")
  af_ok <- tryCatch({
    requireNamespace("amatrix.arrayfire", quietly = TRUE) &&
      amatrix.arrayfire::amatrix_arrayfire_is_available()
  }, error = function(e) FALSE)
  if (af_ok) out <- c(out, "arrayfire")
  out
}

ALL_BACKENDS <- detect_backends()
BACKENDS <- if (BACKEND_FILTER == "all") ALL_BACKENDS else {
  if (!BACKEND_FILTER %in% ALL_BACKENDS)
    stop("Backend '", BACKEND_FILTER, "' not available.  Found: ",
         paste(ALL_BACKENDS, collapse = ", "))
  BACKEND_FILTER
}

# ── Core algorithm ───────────────────────────────────────────────────────────
#
# Subspace iteration for rank-k approximation of X.
# Returns: list(U, d, V, rel_err) where rel_err = ||X - UdV'||_F / ||X - X_opt||_F
# rel_err close to 1.0 means GPU result matches optimal rank-k.

.make_low_rank <- function(n, p, k, seed) {
  set.seed(seed)
  U0 <- qr.Q(qr(matrix(rnorm(n * k), n, k)))
  V0 <- qr.Q(qr(matrix(rnorm(p * k), p, k)))
  X  <- U0 %*% diag(seq(k, 1)) %*% t(V0) +
        matrix(rnorm(n * p, sd = 0.05), n, p)
  storage.mode(X) <- "double"
  X
}

subspace_iter <- function(X_am, k, q, X_host) {
  p <- ncol(X_am)
  n <- nrow(X_am)
  oversampling <- min(10L, p - k)
  k_over <- k + oversampling

  # Random starting block
  Omega <- matrix(rnorm(p * k_over), p, k_over)

  # Power iteration: each step applies X then X'
  Y <- X_am %*% Omega                          # forward:  n × k_over
  for (i in seq_len(q)) {
    Z <- crossprod(X_am, Y)                    # backward: p × k_over
    Y <- X_am %*% qr.Q(qr(as.matrix(Z)))      # forward again, reortho'd
  }

  Q  <- qr.Q(qr(as.matrix(Y)))                # n × k_over orthonormal basis
  B  <- t(Q) %*% X_host                       # k_over × p  (small, on CPU)
  sv <- svd(B, nu = k, nv = k)
  U  <- Q %*% sv$u
  d  <- sv$d[seq_len(k)]
  V  <- sv$v

  # Reconstruction error vs optimal rank-k SVD
  X_rec  <- U %*% diag(d) %*% t(V)
  sv_opt <- base::svd(X_host, nu = k, nv = k)
  X_opt  <- sv_opt$u %*% diag(sv_opt$d[seq_len(k)]) %*% t(sv_opt$v)
  rel_err <- norm(X_host - X_rec, "F") / norm(X_host - X_opt, "F")

  list(d = d, rel_err = rel_err)
}

# ── Timing helper ─────────────────────────────────────────────────────────────

time_ms <- function(fn, n = N_TRIALS) {
  times <- numeric(n)
  for (i in seq_len(n)) {
    t0 <- Sys.time()
    fn()
    times[[i]] <- as.numeric(difftime(Sys.time(), t0, units = "secs")) * 1e3
  }
  sort(times)[[max(1L, as.integer(n / 2L))]]  # median
}

# ── Output helpers ────────────────────────────────────────────────────────────

.hr <- function(char = "─", width = 76) cat(strrep(char, width), "\n")
.status <- function(val, warn, fail) {
  if (val >= fail) "\033[31mFAIL\033[0m"
  else if (val >= warn) "\033[33mWARN\033[0m"
  else "\033[32mOK\033[0m"
}
.speedup_str <- function(cpu_ms, ms) {
  if (is.na(cpu_ms) || identical(cpu_ms, ms)) return("  1.0×")
  if (cpu_ms < 0.5 || ms < 0.5) return("  <1ms")
  sprintf("%5.1f×", cpu_ms / ms)
}

# ── Main loop ─────────────────────────────────────────────────────────────────

cat(sprintf(
  "\nSubspace iteration  k=%d  q=%d  trials=%d\n",
  K_RANK, Q_ITER, N_TRIALS
))
cat("Backends:", paste(BACKENDS, collapse = ", "), "\n\n")

all_pass <- TRUE
results  <- list()

for (sz_name in names(SIZES)) {
  sz  <- SIZES[[sz_name]]
  n   <- sz$n; p <- sz$p
  .hr()
  cat(sprintf("  %s  (%s)  k=%d  q=%d\n", sz_name, sz$label, K_RANK, Q_ITER))
  .hr()
  cat(sprintf("  %-12s  %8s  %8s  %7s  %s\n",
              "backend", "med_ms", "speedup", "rel_err", "status"))

  X_host <- .make_low_rank(n, p, K_RANK, seed = 7L)
  cpu_med <- NA_real_

  for (backend in BACKENDS) {
    X_am <- adgeMatrix(X_host, preferred_backend = backend, precision = "fast")

    fn <- local({
      X_am_ <- X_am; k_ <- K_RANK; q_ <- Q_ITER; Xh_ <- X_host
      function() subspace_iter(X_am_, k = k_, q = q_, X_host = Xh_)
    })
    result <- tryCatch({
      fn()                    # warmup (not timed)
      med <- time_ms(fn)
      res <- fn()             # one final run to capture rel_err
      list(med_ms = med, rel_err = res$rel_err, error = NULL)
    }, error = function(e) list(med_ms = NA, rel_err = NA, error = conditionMessage(e)))

    if (!is.null(result$error)) {
      cat(sprintf("  %-12s  %8s  %8s  %7s  \033[31mERROR: %s\033[0m\n",
                  backend, "—", "—", "—", result$error))
      all_pass <- FALSE
      next
    }

    if (identical(backend, "cpu")) cpu_med <- result$med_ms
    status_str <- .status(result$rel_err, REL_ERR_WARN, REL_ERR_FAIL)
    if (result$rel_err >= REL_ERR_FAIL) all_pass <- FALSE

    cat(sprintf("  %-12s  %8.1f  %8s  %7.3f  %s\n",
                backend,
                result$med_ms,
                .speedup_str(cpu_med, result$med_ms),
                result$rel_err,
                status_str))

    results[[paste(sz_name, backend, sep = "/")]] <- c(
      list(size = sz_name, backend = backend),
      result
    )
  }
  cat("\n")
}

.hr("═")
cat(if (all_pass) "\033[32m✓ ALL PASS\033[0m\n" else "\033[31m✗ FAILURES DETECTED\033[0m\n")
.hr("═")
cat("\n")
cat("rel_err: reconstruction error relative to optimal rank-k SVD.\n")
cat("  1.0 = perfect.  < 1.5 = OK for float32.  > 3.0 = precision problem.\n")
cat("speedup: median time vs cpu backend at this size.\n\n")

if (!all_pass) quit(status = 1L)
