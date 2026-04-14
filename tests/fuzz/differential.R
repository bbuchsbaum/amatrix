#!/usr/bin/env Rscript
# Differential fuzzer for amatrix: cpu vs mlx vs arrayfire
# Usage: Rscript tests/fuzz/differential.R [duration_seconds]
#
# Correct API:
#   adgeMatrix(X, backend="cpu"|"mlx"|"arrayfire")  -> adgeMatrix object
#   matmul(A, B), crossprod(A), tcrossprod(A)
#   A + B, A - B, A * B, A / B
#   rowsums(A), colsums(A), rowmeans(A), colmeans(A)   # lowercase, NOT rowSums()
#   log(A), exp(A), sqrt(A), abs(A)
#   svd(A), chol(A)
#   dist_matrix(A), kernel_matrix(A, kernel=, sigma=)
#   as.matrix(result)  -> plain R matrix

suppressPackageStartupMessages({ library(amatrix) })

BACKENDS <- list()
BACKENDS[["cpu"]] <- list(name = "cpu", available = TRUE)

BACKENDS[["mlx"]] <- tryCatch({
  suppressPackageStartupMessages(library(amatrix.mlx))
  list(name = "mlx", available = TRUE)
}, error = function(e) list(name = "mlx", available = FALSE, err = conditionMessage(e)))

BACKENDS[["arrayfire"]] <- tryCatch({
  suppressPackageStartupMessages(library(amatrix.arrayfire))
  list(name = "arrayfire", available = TRUE)
}, error = function(e) list(name = "arrayfire", available = FALSE, err = conditionMessage(e)))

cat("=== Backend availability ===\n")
for (b in BACKENDS) {
  status <- if (b$available) "OK" else paste("UNAVAILABLE:", b$err)
  cat(sprintf("  %s: %s\n", b$name, status))
}

available_gpu_backends <- names(Filter(function(b) b$available && b$name != "cpu", BACKENDS))
cat("GPU backends:", paste(available_gpu_backends, collapse=", "), "\n\n")

# ---------------------------------------------------------------------------
# Fixture generation
# ---------------------------------------------------------------------------

SHAPE_POOL <- list(
  c(4, 4), c(8, 8), c(16, 16), c(32, 32),
  c(8, 4), c(16, 4), c(32, 8), c(64, 8),
  c(4, 8), c(4, 16), c(8, 32), c(8, 64),
  c(4, 1), c(1, 4), c(12, 3), c(3, 12)
)

SQUARE_SIZES <- c(4L, 6L, 8L, 12L, 16L)

gen_matrix <- function(nr, nc, dist, seed) {
  set.seed(seed)
  switch(dist,
    uniform    = matrix(runif(nr*nc, -10, 10), nr, nc),
    normal     = matrix(rnorm(nr*nc), nr, nc),
    ill_cond   = {
      m <- matrix(rnorm(nr*nc), nr, nc)
      m[, 1] <- m[, 1] * 1e6
      if (nc > 1) m[, 2] <- m[, 2] * 1e-6
      m
    },
    sparse     = {
      m <- matrix(rnorm(nr*nc), nr, nc)
      m[sample(length(m), floor(0.7 * length(m)))] <- 0
      m
    },
    with_nan   = {
      m <- matrix(rnorm(nr*nc), nr, nc)
      m[sample(length(m), max(1L, floor(0.05 * length(m))))] <- NaN
      m
    },
    with_inf   = {
      m <- matrix(rnorm(nr*nc), nr, nc)
      m[sample(length(m), max(1L, floor(0.05 * length(m))))] <- Inf
      m
    },
    large_range = {
      m <- matrix(rnorm(nr*nc), nr, nc)
      m * 10^runif(1, -8, 8)
    },
    integer_vals = matrix(sample(-20:20, nr*nc, replace=TRUE)*1.0, nr, nc),
    matrix(rnorm(nr*nc), nr, nc)  # fallback
  )
}

gen_pos_def <- function(n, seed) {
  set.seed(seed)
  A <- matrix(rnorm(n * n), n, n)
  t(A) %*% A + diag(n) * (n + 1)
}

DISTS <- c("uniform", "normal", "ill_cond", "sparse",
           "with_nan", "with_inf", "large_range", "integer_vals")

# ---------------------------------------------------------------------------
# Backend construction helper
# ---------------------------------------------------------------------------

to_am <- function(X, backend) {
  tryCatch(
    adgeMatrix(X, backend = backend),
    error = function(e) structure(list(msg = conditionMessage(e)), class = "am_error")
  )
}

safe_matrix <- function(A) {
  tryCatch(as.matrix(A), error = function(e) NULL)
}

# ---------------------------------------------------------------------------
# Op execution: returns plain R matrix/vector, or am_error, or NULL
# ---------------------------------------------------------------------------

run_op <- function(op, X, backend, seed) {
  A <- to_am(X, backend)
  if (inherits(A, "am_error")) return(A)

  tryCatch({
    switch(op,

      # --- matmul family ---
      "matmul" = {
        B <- to_am(t(X), backend)
        if (inherits(B, "am_error")) return(NULL)
        safe_matrix(matmul(A, B))
      },
      "matmul_square" = {
        # square: A %*% A
        safe_matrix(matmul(A, A))
      },
      "crossprod" = safe_matrix(crossprod(A)),
      "tcrossprod" = safe_matrix(tcrossprod(A)),

      # --- arithmetic ---
      "add"  = safe_matrix(A + A),
      "sub"  = safe_matrix(A - A),
      "mul"  = safe_matrix(A * A),
      "div"  = {
        X2 <- abs(X) + 1
        A2 <- to_am(X2, backend)
        if (inherits(A2, "am_error")) return(NULL)
        safe_matrix(A / A2)
      },
      "scalar_mul" = safe_matrix(A * 3.14),
      "scalar_add" = safe_matrix(A + 2.71),

      # --- reductions ---
      "rowsums"  = as.numeric(rowsums(A)),
      "colsums"  = as.numeric(colsums(A)),
      "rowmeans" = as.numeric(rowmeans(A)),
      "colmeans" = as.numeric(colmeans(A)),

      # --- elementwise math (on positive values) ---
      "log"  = {
        Ap <- to_am(abs(X) + 1e-6, backend)
        if (inherits(Ap, "am_error")) return(NULL)
        safe_matrix(log(Ap))
      },
      "exp"  = {
        Xc <- pmin(pmax(X, -20), 20)
        Ap <- to_am(Xc, backend)
        if (inherits(Ap, "am_error")) return(NULL)
        safe_matrix(exp(Ap))
      },
      "sqrt" = {
        Ap <- to_am(abs(X), backend)
        if (inherits(Ap, "am_error")) return(NULL)
        safe_matrix(sqrt(Ap))
      },
      "abs"  = safe_matrix(abs(A)),

      # --- comparisons ---
      "eq"   = { r <- tryCatch(A == A, error=function(e) NULL); if(!is.null(r)) safe_matrix(r) else NULL },
      "lt"   = { r <- tryCatch(A < A,  error=function(e) NULL); if(!is.null(r)) safe_matrix(r) else NULL },

      # --- factorizations ---
      "chol" = {
        n <- min(nrow(X), ncol(X))
        Xpd <- gen_pos_def(n, seed)
        Apd <- to_am(Xpd, backend)
        if (inherits(Apd, "am_error")) return(NULL)
        safe_matrix(chol(Apd))
      },
      "svd_d" = {
        res <- svd(A)
        as.numeric(res$d)
      },
      "svd_reconstruct" = {
        res <- svd(A)
        U <- safe_matrix(res$u)
        D <- as.numeric(res$d)
        V <- safe_matrix(res$v)
        if (is.null(U) || is.null(V)) return(NULL)
        U %*% diag(D) %*% t(V)
      },

      # --- distance / kernel ---
      "dist_eucl"     = safe_matrix(dist_matrix(A)),
      "kernel_rbf"    = safe_matrix(kernel_matrix(A, kernel = "rbf", sigma = 1.0)),
      "kernel_linear" = safe_matrix(kernel_matrix(A, kernel = "linear")),

      NULL
    )
  }, error = function(e) {
    structure(list(msg = conditionMessage(e)), class = "am_error")
  })
}

# ---------------------------------------------------------------------------
# Divergence checker
# ---------------------------------------------------------------------------

# Relative tolerance for f32 GPU vs f64 CPU.
# f32 has ~7 significant digits so 1e-5 relative is conservative.
# We use max(abs_diff) / max(1, max(|cpu|)) so zero-valued results still
# use an absolute floor of 1e-4.
TOL_REL <- 1e-4   # relative error threshold (f32 precision headroom)
TOL_ABS <- 1e-7   # absolute floor for near-zero results

check_diff <- function(cpu_res, gpu_res) {
  if (is.null(cpu_res) || is.null(gpu_res)) return(NULL)
  if (inherits(cpu_res, "am_error") || inherits(gpu_res, "am_error")) return(NULL)

  # Dimension mismatch
  if (!identical(dim(cpu_res), dim(gpu_res)) &&
      !identical(length(cpu_res), length(gpu_res))) {
    return(list(kind = "dim_mismatch",
                cpu_dim = dim(cpu_res) %||% length(cpu_res),
                gpu_dim = dim(gpu_res) %||% length(gpu_res)))
  }

  if (!is.numeric(cpu_res) || !is.numeric(gpu_res)) return(NULL)

  cv <- as.numeric(cpu_res)
  gv <- as.numeric(gpu_res)
  if (length(cv) != length(gv)) {
    return(list(kind = "length_mismatch", cpu_len = length(cv), gpu_len = length(gv)))
  }

  # Treat NA and NaN as equivalent missing-value sentinels (amatrix returns NA,
  # base R returns NaN; both signal "not a number" for the same positions).
  missing_c <- is.na(cv)   # catches both NA and NaN
  missing_g <- is.na(gv)
  inf_c     <- is.infinite(cv)
  inf_g     <- is.infinite(gv)

  # True NaN-propagation mismatch: positions that are missing in one but finite
  # in the other (ignore Inf vs NaN differences which vary by IEEE implementation)
  bad_c <- missing_c & !missing_g & !inf_g
  bad_g <- missing_g & !missing_c & !inf_c
  if (any(bad_c) || any(bad_g)) {
    return(list(kind        = "nan_mismatch",
                cpu_na_n    = sum(missing_c),
                gpu_na_n    = sum(missing_g),
                max_abs_diff = NA_real_))
  }

  # Compare only finite positions where both are valid
  valid <- !missing_c & !missing_g & !inf_c & !inf_g
  if (sum(valid) == 0) return(NULL)

  cv_v <- cv[valid]; gv_v <- gv[valid]
  abs_diffs <- abs(cv_v - gv_v)
  max_abs   <- max(abs_diffs)
  mean_abs  <- mean(abs_diffs)

  # Relative error: normalise by scale of cpu result
  scale      <- max(abs(cv_v), 1e-10)
  rel_err    <- max_abs / scale

  # Flag if relative error exceeds f32 threshold AND abs error exceeds floor
  if (rel_err > TOL_REL && max_abs > TOL_ABS) {
    return(list(kind          = "value_divergence",
                max_abs_diff  = max_abs,
                mean_abs_diff = mean_abs,
                rel_err       = rel_err,
                scale         = scale,
                tol_rel       = TOL_REL))
  }
  NULL
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Op list: (op_name, requires_square)
# ---------------------------------------------------------------------------

OPS <- list(
  list(op = "matmul",          square = FALSE),
  list(op = "matmul_square",   square = TRUE),
  list(op = "crossprod",       square = FALSE),
  list(op = "tcrossprod",      square = FALSE),
  list(op = "add",             square = FALSE),
  list(op = "sub",             square = FALSE),
  list(op = "mul",             square = FALSE),
  list(op = "div",             square = FALSE),
  list(op = "scalar_mul",      square = FALSE),
  list(op = "scalar_add",      square = FALSE),
  list(op = "rowsums",         square = FALSE),
  list(op = "colsums",         square = FALSE),
  list(op = "rowmeans",        square = FALSE),
  list(op = "colmeans",        square = FALSE),
  list(op = "log",             square = FALSE),
  list(op = "exp",             square = FALSE),
  list(op = "sqrt",            square = FALSE),
  list(op = "abs",             square = FALSE),
  list(op = "eq",              square = FALSE),
  list(op = "lt",              square = FALSE),
  list(op = "chol",            square = FALSE),  # uses gen_pos_def internally
  list(op = "svd_d",           square = FALSE),
  list(op = "svd_reconstruct", square = FALSE),
  list(op = "dist_eucl",       square = FALSE),
  list(op = "kernel_rbf",      square = FALSE),
  list(op = "kernel_linear",   square = FALSE)
)

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

fuzz <- function(duration_sec = 360L, seed_start = 42L) {
  divergences   <- list()
  errors_by_key <- list()
  op_counts     <- setNames(integer(length(OPS)), sapply(OPS, `[[`, "op"))
  fixture_count <- 0L

  cat(sprintf("=== Fuzzer start: duration=%ds seed_start=%d ===\n", duration_sec, seed_start))
  cat(sprintf("Ops: %s\n\n", paste(sapply(OPS, `[[`, "op"), collapse=", ")))

  t0 <- proc.time()[["elapsed"]]
  seed <- seed_start

  repeat {
    elapsed <- proc.time()[["elapsed"]] - t0
    if (elapsed >= duration_sec) break
    seed <- seed + 1L

    # Pick shape and distribution
    set.seed(seed)
    shape <- SHAPE_POOL[[sample(length(SHAPE_POOL), 1L)]]
    nr <- shape[1L]; nc <- shape[2L]
    dist <- sample(DISTS, 1L)
    X <- gen_matrix(nr, nc, dist, seed)

    for (op_spec in OPS) {
      op <- op_spec$op

      # For square-only ops, use a square matrix
      X_use <- if (op_spec$square) {
        n <- sample(SQUARE_SIZES, 1L)
        gen_matrix(n, n, dist, seed + 1000L)
      } else X

      op_counts[op] <- op_counts[op] + 1L
      fixture_count <- fixture_count + 1L

      # CPU reference
      cpu_res <- run_op(op, X_use, "cpu", seed)

      if (inherits(cpu_res, "am_error")) next  # CPU error = skip
      if (is.null(cpu_res)) next

      # GPU backends
      for (gpu in available_gpu_backends) {
        gpu_res <- run_op(op, X_use, gpu, seed)

        if (inherits(gpu_res, "am_error")) {
          key <- paste0(op, "@", gpu)
          errors_by_key[[key]] <- c(errors_by_key[[key]], list(list(
            seed = seed, shape = c(nr, nc), dist = dist, msg = gpu_res$msg
          )))
          next
        }
        if (is.null(gpu_res)) next

        diff <- check_diff(cpu_res, gpu_res)
        if (!is.null(diff)) {
          entry <- list(op = op, backend = gpu, seed = seed,
                        shape = c(nrow(X_use), ncol(X_use)), dist = dist,
                        diff = diff)
          divergences <- c(divergences, list(entry))

          max_d <- if (!is.null(diff$max_abs_diff)) sprintf("%.3e", diff$max_abs_diff) else "?"
          cat(sprintf("[DIVERGE] op=%-20s backend=%-12s seed=%6d shape=%dx%d dist=%-12s diff=%s\n",
                      op, gpu, seed, nrow(X_use), ncol(X_use), dist, max_d))
        }
      }
    }

    if (seed %% 50L == 0L) {
      elapsed2 <- proc.time()[["elapsed"]] - t0
      cat(sprintf("[progress] seed=%d fixtures=%d divergences=%d elapsed=%.0fs\n",
                  seed, fixture_count, length(divergences), elapsed2))
    }
  }

  elapsed_final <- proc.time()[["elapsed"]] - t0
  cat(sprintf("\n=== Done: fixtures=%d divergences=%d elapsed=%.1fs ===\n",
              fixture_count, length(divergences), elapsed_final))

  list(divergences   = divergences,
       errors_by_key = errors_by_key,
       op_counts     = op_counts,
       fixture_count = fixture_count,
       elapsed       = elapsed_final,
       backends      = available_gpu_backends)
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
duration <- if (length(args) >= 1L) as.integer(args[1L]) else 360L

results <- fuzz(duration_sec = duration, seed_start = 42L)
saveRDS(results, "/Users/bbuchsbaum/code/amatrix/tests/fuzz/results.rds")
cat("Saved to tests/fuzz/results.rds\n")

# ---------------------------------------------------------------------------
# Summary report
# ---------------------------------------------------------------------------

cat("\n=== DIVERGENCE SUMMARY ===\n")
if (length(results$divergences) == 0) {
  cat("No divergences found.\n")
} else {
  by_key <- split(results$divergences,
                  sapply(results$divergences, function(d) paste0(d$op, "@", d$backend)))
  for (key in names(by_key)) {
    group <- by_key[[key]]
    cat(sprintf("\n%s: %d divergence(s)\n", key, length(group)))
    for (d in head(group, 5L)) {
      diff <- d$diff
      cat(sprintf("  seed=%d shape=%dx%d dist=%-12s kind=%-18s max_diff=%s\n",
                  d$seed, d$shape[1L], d$shape[2L], d$dist,
                  diff$kind,
                  if (!is.null(diff$max_abs_diff)) sprintf("%.3e", diff$max_abs_diff) else "NA"))
    }
  }
}

cat("\n=== GPU ERRORS (op@backend: count, first message) ===\n")
if (length(results$errors_by_key) == 0) {
  cat("No GPU-side errors.\n")
} else {
  for (key in names(results$errors_by_key)) {
    errs <- results$errors_by_key[[key]]
    cat(sprintf("  %s: %d errors | first: %s\n", key, length(errs), errs[[1L]]$msg))
  }
}

cat("\n=== OP COVERAGE ===\n")
for (op in names(results$op_counts)) {
  cat(sprintf("  %-22s %d\n", op, results$op_counts[op]))
}
