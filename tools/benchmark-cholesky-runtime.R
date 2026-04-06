#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(amatrix)
  library(amatrix.mlx)
})

if (!amatrix_mlx_native_available()) {
  stop("benchmark-cholesky-runtime.R requires the native MLX backend", call. = FALSE)
}

old_available <- getOption("amatrix.mlx.available")
options(amatrix.mlx.available = TRUE)
on.exit(options(amatrix.mlx.available = old_available), add = TRUE)

benchmark_elapsed <- function(fn, reps = 3L, iterations = 5L, warmup = NULL) {
  if (is.function(warmup)) {
    warmup()
  }

  timings <- numeric(reps)
  for (idx in seq_len(reps)) {
    gc()
    timings[[idx]] <- system.time({
      for (iter in seq_len(iterations)) {
        fn()
      }
    })[["elapsed"]] / iterations
  }
  median(timings)
}

frob_norm <- function(x) {
  sqrt(sum(x * x))
}

make_ridge_spd_case <- function(n_obs = 4096L, p = 768L, rhs_cols = 64L, lambda = 0.75, seed = 20260406L) {
  set.seed(seed)
  X <- matrix(rnorm(n_obs * p), nrow = n_obs, ncol = p)
  B <- matrix(rnorm(p * rhs_cols), nrow = p, ncol = rhs_cols)
  list(
    workload = "ridge_spd",
    case = sprintf("%dx%d", p, p),
    rhs_cols = rhs_cols,
    A = crossprod(X) + diag(lambda, p),
    B = B
  )
}

make_kernel_spd_case <- function(n = 640L, p = 12L, rhs_cols = 32L, sigma = 1.1, jitter = 0.2, seed = 20260407L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  B <- matrix(rnorm(n * rhs_cols), nrow = n, ncol = rhs_cols)
  list(
    workload = "kernel_spd",
    case = sprintf("%dx%d", n, n),
    rhs_cols = rhs_cols,
    A = am_kernel(X, kernel = "rbf", sigma = sigma) + diag(jitter, n),
    B = B
  )
}

make_amatrix <- function(A, backend = c("cpu", "mlx")) {
  backend <- match.arg(backend)
  if (identical(backend, "mlx")) {
    return(adgeMatrix(A, preferred_backend = "mlx", precision = "fast"))
  }
  adgeMatrix(A, precision = "strict")
}

benchmark_case <- function(case) {
  iterations <- if (nrow(case$A) >= 768L) 3L else 5L
  fac_cpu <- am_chol_factor(make_amatrix(case$A, "cpu"))
  fac_mlx <- am_chol_factor(make_amatrix(case$A, "mlx"))
  ref_sol <- solve(case$A, case$B)
  sol_cpu <- am_chol_solve(fac_cpu, case$B)
  sol_mlx <- am_chol_solve(fac_mlx, case$B)

  factor_rows <- lapply(c("cpu", "mlx"), function(backend_name) {
    fac <- if (identical(backend_name, "cpu")) fac_cpu else fac_mlx
    elapsed <- benchmark_elapsed(
      function() {
        out <- am_chol_factor(make_amatrix(case$A, backend_name))
        invisible(out)
      },
      iterations = iterations,
      warmup = if (identical(backend_name, "mlx")) {
        function() invisible(am_chol_factor(make_amatrix(case$A, backend_name)))
      } else {
        NULL
      }
    )

    data.frame(
      workload = case$workload,
      case = case$case,
      rhs_cols = case$rhs_cols,
      phase = "factor",
      runtime = backend_name,
      elapsed = elapsed,
      rel_error = frob_norm(crossprod(fac@factor) - case$A) / frob_norm(case$A),
      stringsAsFactors = FALSE
    )
  })

  solve_rows <- lapply(c("cpu", "mlx"), function(backend_name) {
    fac <- if (identical(backend_name, "cpu")) fac_cpu else fac_mlx
    sol <- if (identical(backend_name, "cpu")) sol_cpu else sol_mlx
    elapsed <- benchmark_elapsed(
      function() {
        out <- am_chol_solve(fac, case$B)
        invisible(out)
      },
      iterations = max(5L, iterations),
      warmup = if (identical(backend_name, "mlx")) {
        function() invisible(am_chol_solve(fac, case$B))
      } else {
        NULL
      }
    )

    data.frame(
      workload = case$workload,
      case = case$case,
      rhs_cols = case$rhs_cols,
      phase = "batched_solve",
      runtime = backend_name,
      elapsed = elapsed,
      rel_error = frob_norm(sol - ref_sol) / frob_norm(ref_sol),
      stringsAsFactors = FALSE
    )
  })

  total_rows <- lapply(c("cpu", "mlx"), function(backend_name) {
    elapsed <- benchmark_elapsed(
      function() {
        fac <- am_chol_factor(make_amatrix(case$A, backend_name))
        out <- am_chol_solve(fac, case$B)
        invisible(out)
      },
      iterations = iterations,
      warmup = if (identical(backend_name, "mlx")) {
        function() {
          fac <- am_chol_factor(make_amatrix(case$A, backend_name))
          invisible(am_chol_solve(fac, case$B))
        }
      } else {
        NULL
      }
    )

    data.frame(
      workload = case$workload,
      case = case$case,
      rhs_cols = case$rhs_cols,
      phase = "factor_plus_batched_solve",
      runtime = backend_name,
      elapsed = elapsed,
      rel_error = if (identical(backend_name, "cpu")) {
        frob_norm(sol_cpu - ref_sol) / frob_norm(ref_sol)
      } else {
        frob_norm(sol_mlx - ref_sol) / frob_norm(ref_sol)
      },
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, c(factor_rows, solve_rows, total_rows))
}

cases <- list(
  make_ridge_spd_case(),
  make_kernel_spd_case()
)

results <- do.call(rbind, lapply(cases, benchmark_case))
row.names(results) <- NULL
results$elapsed <- sprintf("%.6f", results$elapsed)
results$rel_error <- ifelse(
  is.na(results$rel_error),
  NA_character_,
  sprintf("%.3e", results$rel_error)
)

cat("Notes:\n")
cat("- factor benchmarks are cold builds on fresh adgeMatrix objects, so am_chol_factor() cache reuse is not counted.\n")
cat("- batched_solve benchmarks reuse a precomputed factor and isolate the many-RHS triangular-solve path.\n")
cat("- ridge_spd uses crossprod(X) + lambda*I; kernel_spd uses an RBF kernel matrix plus diagonal jitter.\n")
cat("- rel_error is the Cholesky reconstruction residual for factor rows and the relative solution error versus CPU solve for solve rows.\n\n")
print(results)
