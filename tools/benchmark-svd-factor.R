#!/usr/bin/env Rscript

r_string <- function(x) {
  paste0("\"", gsub("([\"\\\\])", "\\\\\\1", x), "\"")
}

bench_lib <- Sys.getenv("AMATRIX_BENCH_LIB", "")
if (nzchar(bench_lib)) {
  .libPaths(c(normalizePath(bench_lib), .libPaths()))
}

script_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", script_args, value = TRUE)
direct_file_entry <- length(script_arg) > 0L
script_path <- if (direct_file_entry) {
  normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path("tools", "benchmark-svd-factor.R"), mustWork = FALSE)
}
repo_root <- normalizePath(
  if (direct_file_entry) dirname(dirname(script_path)) else ".",
  mustWork = TRUE
)
safe_mlx_command <- sprintf(
  "Rscript -e 'Sys.setenv(AMATRIX_MLX_PROBE_GPU = \"1\", AMATRIX_SVD_FACTOR_ALLOW_NESTED_MLX = \"true\", AMATRIX_SVD_FACTOR_MLX_MAIN = \"true\"); setwd(%s); source(\"tools/benchmark-svd-factor.R\")'",
  r_string(repo_root)
)

allow_nested_mlx <- identical(
  tolower(Sys.getenv("AMATRIX_SVD_FACTOR_ALLOW_NESTED_MLX", "false")),
  "true"
)
mlx_main_process <- identical(
  tolower(Sys.getenv("AMATRIX_SVD_FACTOR_MLX_MAIN", "false")),
  "true"
)

if (allow_nested_mlx && direct_file_entry && !mlx_main_process) {
  stop(
    paste0(
      "Direct file-entry MLX launch is unstable on this machine.\nUse:\n  ",
      safe_mlx_command
    ),
    call. = FALSE
  )
}

benchmark_mlx_in_process <- allow_nested_mlx && (mlx_main_process || !direct_file_entry)

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(amatrix)
  }
})

load_optional_backend <- function(pkg) {
  source_dir <- file.path(repo_root, "backends", pkg)
  if (requireNamespace("pkgload", quietly = TRUE) && dir.exists(source_dir)) {
    pkgload::load_all(source_dir, quiet = TRUE)
    return(TRUE)
  }

  if (requireNamespace(pkg, quietly = TRUE)) {
    invisible(loadNamespace(pkg))
    return(TRUE)
  }

  FALSE
}

mlx_source_dir <- file.path(repo_root, "backends", "amatrix.mlx")
have_mlx <- if (benchmark_mlx_in_process) {
  dir.exists(mlx_source_dir) || requireNamespace("amatrix.mlx", quietly = TRUE)
} else {
  FALSE
}

if (!requireNamespace("irlba", quietly = TRUE)) {
  stop("Package 'irlba' is required for this benchmark", call. = FALSE)
}

options(
  amatrix.optional_backends = FALSE,
  amatrix.mlx.available = have_mlx
)

if (have_mlx) {
  load_optional_backend("amatrix.mlx")
}

clear_factor_cache <- function() {
  cache <- amatrix:::.amatrix_state$model_cache
  keys <- ls(envir = cache, all.names = TRUE)
  if (length(keys) > 0L) {
    rm(list = keys, envir = cache)
  }
  invisible(NULL)
}

benchmark_elapsed <- function(fn, reps = 3L, warmup = NULL) {
  if (is.function(warmup)) {
    warmup()
  }
  timings <- numeric(reps)
  for (idx in seq_len(reps)) {
    clear_factor_cache()
    gc()
    timings[[idx]] <- system.time(invisible(fn()))[["elapsed"]]
  }
  median(timings)
}

relative_sv_error <- function(actual, expected) {
  max(abs(actual - expected) / pmax(abs(expected), 1e-12))
}

measure_mlx_warmup <- function() {
  if (!have_mlx) {
    return(NA_real_)
  }

  x <- adgeMatrix(
    matrix(rnorm(32L * 16L), nrow = 32L, ncol = 16L),
    preferred_backend = "mlx",
    precision = "fast"
  )
  clear_factor_cache()
  gc()
  system.time(
    invisible(amatrix:::svd_factor(
      x,
      k = 5L,
      method = "rsvd",
      n_oversamples = 4L,
      n_iter = 1L
    ))
  )[["elapsed"]]
}

benchmark_factor_path <- function(host,
                                  k,
                                  n_oversamples,
                                  n_iter,
                                  preferred_backend,
                                  method,
                                  phase,
                                  reference_d,
                                  warmup = NULL,
                                  reps = 3L) {
  precision <- if (identical(preferred_backend, "cpu")) "strict" else "fast"
  result <- NULL

  elapsed <- benchmark_elapsed(
    fn = function() {
      x <- adgeMatrix(
        host,
        preferred_backend = preferred_backend,
        precision = precision
      )
      result <<- amatrix:::svd_factor(
        x,
        k = k,
        method = method,
        n_oversamples = n_oversamples,
        n_iter = n_iter
      )
    },
    reps = reps,
    warmup = warmup
  )

  data.frame(
    implementation = sprintf("svd_factor(%s)", method),
    preferred_backend = preferred_backend,
    phase = phase,
    selected_method = result@method,
    selected_engine = result@engine,
    selected_backend = result@backend,
    elapsed = elapsed,
    rel_sv_err = relative_sv_error(result@d, reference_d),
    stringsAsFactors = FALSE
  )
}

benchmark_case <- function(n, p, k = 20L, n_oversamples = 10L, n_iter = 2L) {
  set.seed(20260405L + n + p + k)
  host <- matrix(rnorm(n * p), nrow = n, ncol = p)

  exact_elapsed <- system.time(
    exact_ref <- base::svd(host, nu = k, nv = k)
  )[["elapsed"]]
  reference_d <- exact_ref$d[seq_len(k)]

  svdr_result <- NULL
  svdr_elapsed <- benchmark_elapsed(
    fn = function() {
      svdr_result <<- irlba::svdr(
        host,
        k = k,
        extra = n_oversamples,
        it = n_iter
      )
    },
    reps = 3L
  )

  rows <- list(
    data.frame(
      implementation = "base::svd",
      preferred_backend = "cpu",
      phase = "reference",
      selected_method = "exact",
      selected_engine = "exact_svd",
      selected_backend = "cpu",
      elapsed = exact_elapsed,
      rel_sv_err = 0,
      stringsAsFactors = FALSE
    ),
    data.frame(
      implementation = "irlba::svdr",
      preferred_backend = "cpu",
      phase = "reference",
      selected_method = "rsvd",
      selected_engine = "irlba_svdr",
      selected_backend = "cpu",
      elapsed = svdr_elapsed,
      rel_sv_err = relative_sv_error(svdr_result$d[seq_len(k)], reference_d),
      stringsAsFactors = FALSE
    ),
    benchmark_factor_path(
      host = host,
      k = k,
      n_oversamples = n_oversamples,
      n_iter = n_iter,
      preferred_backend = "cpu",
      method = "exact",
      phase = "reference",
      reference_d = reference_d
    )
  )

  if (have_mlx) {
    rows[[length(rows) + 1L]] <- benchmark_factor_path(
      host = host,
      k = k,
      n_oversamples = n_oversamples,
      n_iter = n_iter,
      preferred_backend = "mlx",
      method = "rsvd",
      phase = "steady_state",
      reference_d = reference_d
    )
    rows[[length(rows) + 1L]] <- benchmark_factor_path(
      host = host,
      k = k,
      n_oversamples = n_oversamples,
      n_iter = n_iter,
      preferred_backend = "mlx",
      method = "auto",
      phase = "steady_state",
      reference_d = reference_d
    )
  }

  out <- do.call(rbind, rows)
  out$case <- sprintf("%dx%d", n, p)
  out$k <- as.integer(k)
  out$n_oversamples <- as.integer(n_oversamples)
  out$n_iter <- as.integer(n_iter)
  out[, c(
    "case",
    "k",
    "n_oversamples",
    "n_iter",
    "implementation",
    "preferred_backend",
    "phase",
    "selected_method",
    "selected_engine",
    "selected_backend",
    "elapsed",
    "rel_sv_err"
  )]
}

cases <- list(
  c(300L, 240L),
  c(500L, 400L),
  c(1000L, 800L),
  c(2000L, 1600L)
)

mlx_warmup_elapsed <- measure_mlx_warmup()
rows <- do.call(
  rbind,
  lapply(cases, function(case) benchmark_case(case[[1]], case[[2]]))
)

row.names(rows) <- NULL
rows$elapsed <- sprintf("%.3f", rows$elapsed)
rows$rel_sv_err <- sprintf("%.4f", rows$rel_sv_err)

cat("Notes:\n")
cat("- This harness benchmarks rank-k factorization, not downstream projection/reconstruction.\n")
cat("- `selected_method` is the public factorization mode; `selected_engine` is the concrete implementation path.\n")
cat("- `steady_state` means the MLX backend was warmed once before timing to avoid first-call compile cost.\n")
cat(sprintf("- Current auto policy uses `min(dim) >= %d` and `k / min(dim) <= %.2f` before selecting rsvd.\n",
            amatrix:::.amatrix_svd_factor_rsvd_min_dim(),
            amatrix:::.amatrix_svd_factor_rsvd_max_rank_ratio()))
if (have_mlx) {
  cat(sprintf("- One-time MLX warm-up on this run: %.3f s\n", mlx_warmup_elapsed))
} else if (
  !allow_nested_mlx &&
    (dir.exists(mlx_source_dir) ||
       nzchar(system.file(package = "amatrix.mlx")))
) {
  cat("- MLX rows are skipped here by default to keep direct file-entry benchmarking stable on this machine.\n")
  cat("- Re-run with the direct command below to benchmark MLX through the stable `Rscript -e 'source(...)'` launch path.\n")
  cat(sprintf("  %s\n", safe_mlx_command))
}
if (nzchar(bench_lib)) {
  cat(sprintf("- Prepended library path from AMATRIX_BENCH_LIB=%s\n", bench_lib))
}
cat("\n")
print(rows)
