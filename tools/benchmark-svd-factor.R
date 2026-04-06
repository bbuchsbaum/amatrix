#!/usr/bin/env Rscript

bench_lib <- Sys.getenv("AMATRIX_BENCH_LIB", "")
if (nzchar(bench_lib)) {
  .libPaths(c(normalizePath(bench_lib), .libPaths()))
}
repo_root <- normalizePath(".")

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

allow_nested_mlx <- identical(
  tolower(Sys.getenv("AMATRIX_SVD_FACTOR_ALLOW_NESTED_MLX", "false")),
  "true"
)
have_mlx <- if (allow_nested_mlx) load_optional_backend("amatrix.mlx") else FALSE

if (!requireNamespace("irlba", quietly = TRUE)) {
  stop("Package 'irlba' is required for this benchmark", call. = FALSE)
}

options(
  amatrix.mlx.available = have_mlx
)

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

r_string <- function(x) {
  paste0("\"", gsub("([\"\\\\])", "\\\\\\1", x), "\"")
}

run_rscript_expr <- function(code, args = character()) {
  out <- suppressWarnings(system2(
    command = file.path(R.home("bin"), "Rscript"),
    # Launch MLX workers through `Rscript -e`; quote the expression payload so
    # the shell does not split it.
    args = c("-e", shQuote(code), shQuote(args)),
    stdout = TRUE,
    stderr = TRUE
  ))
  list(
    output = out,
    status = attr(out, "status", exact = TRUE)
  )
}

write_worker_payload <- function(payload) {
  path <- tempfile("amatrix-svd-bench-", fileext = ".rds")
  saveRDS(payload, path)
  path
}

parse_worker_result <- function(output) {
  line <- grep("^RESULT\t", output, value = TRUE)
  if (length(line) == 0L) {
    return(NULL)
  }
  fields <- strsplit(tail(line, 1L), "\t", fixed = TRUE)[[1L]]
  fields[-1L]
}

mlx_worker_preamble <- function() {
  setup <- character()
  if (nzchar(bench_lib)) {
    setup <- c(
      setup,
      sprintf(".libPaths(c(%s, .libPaths()));", r_string(normalizePath(bench_lib)))
    )
  }

  setup <- c(
    setup,
    sprintf("pkgload::load_all(%s, quiet = TRUE);", r_string(repo_root)),
    sprintf(
      "if (requireNamespace(\"pkgload\", quietly = TRUE) && dir.exists(%s)) { pkgload::load_all(%s, quiet = TRUE); } else { invisible(loadNamespace(\"amatrix.mlx\")); };",
      r_string(file.path(repo_root, "backends", "amatrix.mlx")),
      r_string(file.path(repo_root, "backends", "amatrix.mlx"))
    ),
    "options(amatrix.mlx.available = TRUE);",
    "payload <- readRDS(commandArgs(trailingOnly = TRUE)[1L]);"
  )

  paste(setup, collapse = " ")
}

measure_mlx_warmup <- function() {
  if (!have_mlx) {
    return(NA_real_)
  }

  payload_path <- write_worker_payload(list(
    x = matrix(rnorm(32L * 16L), nrow = 32L, ncol = 16L)
  ))
  on.exit(unlink(payload_path), add = TRUE)

  code <- paste(
    mlx_worker_preamble(),
    "elapsed <- system.time(invisible(amatrix.mlx:::amatrix_mlx_rsvd(payload$x, k = 5L, n_oversamples = 4L, n_iter = 1L)))[[\"elapsed\"]];",
    "cat(sprintf(\"RESULT\\t%.6f\\n\", unname(elapsed)));"
  )

  res <- run_rscript_expr(code, args = payload_path)
  parsed <- parse_worker_result(res$output)
  if (!is.null(res$status) || is.null(parsed)) {
    warning("MLX warm-up worker failed; steady-state timings may be unavailable.", call. = FALSE)
    return(NA_real_)
  }
  as.numeric(parsed[[1L]])
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
      x <- amatrix::adgeMatrix(
        host,
        preferred_backend = preferred_backend,
        precision = precision
      )
      result <<- amatrix::svd_factor(
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

benchmark_mlx_factor_path <- function(host,
                                      k,
                                      n_oversamples,
                                      n_iter,
                                      method,
                                      reference_d,
                                      reps = 3L) {
  if (!have_mlx) {
    return(NULL)
  }

  payload_path <- write_worker_payload(list(
    host = host,
    k = as.integer(k),
    n_oversamples = as.integer(n_oversamples),
    n_iter = as.integer(n_iter),
    method = method,
    reference_d = reference_d,
    reps = as.integer(reps)
  ))
  on.exit(unlink(payload_path), add = TRUE)

  code <- paste(
    mlx_worker_preamble(),
    "set.seed(20260405L);",
    "warm <- matrix(rnorm(32L * 16L), nrow = 32L, ncol = 16L);",
    "invisible(amatrix.mlx:::amatrix_mlx_rsvd(warm, k = 5L, n_oversamples = 4L, n_iter = 1L));",
    "clear_cache <- function() {",
    "  cache <- amatrix:::.amatrix_state$model_cache;",
    "  keys <- ls(envir = cache, all.names = TRUE);",
    "  if (length(keys) > 0L) rm(list = keys, envir = cache);",
    "  invisible(NULL)",
    "};",
    "timings <- numeric(payload$reps);",
    "err <- NA_real_;",
    "for (idx in seq_len(payload$reps)) {",
    "  clear_cache(); gc();",
    "  timings[[idx]] <- system.time({",
    "    x <- adgeMatrix(payload$host, preferred_backend = \"mlx\", precision = \"fast\");",
    "    fac <- svd_factor(",
    "      x,",
    "      k = payload$k,",
    "      method = payload$method,",
    "      n_oversamples = payload$n_oversamples,",
    "      n_iter = payload$n_iter",
    "    );",
    "    err <- max(abs(fac@d - payload$reference_d) / pmax(abs(payload$reference_d), 1e-12));",
    "  })[[\"elapsed\"]];",
    "};",
    "cat(sprintf(\"RESULT\\t%s\\t%s\\t%s\\t%.6f\\t%.6f\\n\", fac@method, fac@engine, fac@backend, median(timings), err));"
  )

  res <- run_rscript_expr(code, args = payload_path)
  parsed <- parse_worker_result(res$output)
  if (!is.null(res$status) || is.null(parsed)) {
    warning(
      sprintf(
        "MLX worker failed for %s; last output: %s",
        method,
        paste(tail(res$output, 4L), collapse = " | ")
      ),
      call. = FALSE
    )
    return(data.frame(
      implementation = sprintf("svd_factor(%s)", method),
      preferred_backend = "mlx",
      phase = "steady_state",
      selected_method = "unavailable",
      selected_engine = "unavailable",
      selected_backend = "unavailable",
      elapsed = NA_real_,
      rel_sv_err = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    implementation = sprintf("svd_factor(%s)", method),
    preferred_backend = "mlx",
    phase = "steady_state",
    selected_method = parsed[[1L]],
    selected_engine = parsed[[2L]],
    selected_backend = parsed[[3L]],
    elapsed = as.numeric(parsed[[4L]]),
    rel_sv_err = as.numeric(parsed[[5L]]),
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
    rows[[length(rows) + 1L]] <- benchmark_mlx_factor_path(
      host = host,
      k = k,
      n_oversamples = n_oversamples,
      n_iter = n_iter,
      method = "rsvd",
      reference_d = reference_d
    )
    rows[[length(rows) + 1L]] <- benchmark_mlx_factor_path(
      host = host,
      k = k,
      n_oversamples = n_oversamples,
      n_iter = n_iter,
      method = "auto",
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
    (dir.exists(file.path(repo_root, "backends", "amatrix.mlx")) ||
       nzchar(system.file(package = "amatrix.mlx")))
) {
  cat("- MLX rows are skipped here by default because nested/file-entry MLX benchmarking is unstable on this machine.\n")
  cat("- Use the direct `Rscript -e` spot-benchmark commands documented in `docs/gpu-svd-analysis.md` for MLX steady-state numbers.\n")
}
if (nzchar(bench_lib)) {
  cat(sprintf("- Prepended library path from AMATRIX_BENCH_LIB=%s\n", bench_lib))
}
cat("\n")
print(rows)
