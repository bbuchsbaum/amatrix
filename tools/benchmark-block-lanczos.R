#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else if (requireNamespace("amatrix", quietly = TRUE)) {
    library(amatrix)
  } else {
    stop("Either the installed 'amatrix' package or 'pkgload' is required", call. = FALSE)
  }
  if (!requireNamespace("irlba", quietly = TRUE)) {
    stop("Package 'irlba' is required for this benchmark", call. = FALSE)
  }
  if (requireNamespace("pkgload", quietly = TRUE) && dir.exists("backends/amatrix.mlx")) {
    pkgload::load_all("backends/amatrix.mlx", quiet = TRUE)
  } else if (requireNamespace("amatrix.mlx", quietly = TRUE)) {
    library(amatrix.mlx)
  }
})

if (exists("amatrix_mlx_is_available", mode = "function")) {
  options(amatrix.mlx.available = TRUE)
}

benchmark_elapsed <- function(fn, reps = 3L, warmup = NULL) {
  if (is.function(warmup)) {
    warmup()
  }
  timings <- numeric(reps)
  last <- NULL
  for (idx in seq_len(reps)) {
    gc()
    timings[[idx]] <- system.time(last <- fn())[["elapsed"]]
  }
  list(elapsed = median(timings), result = last)
}

benchmark_case <- function(n, p, k, block_size, n_steps, seed) {
  set.seed(seed)
  host <- matrix(rnorm(n * p), nrow = n, ncol = p)
  x_mlx <- adgeMatrix(host, preferred_backend = "mlx", precision = "fast")
  ref_d <- La.svd(host, nu = k, nv = k)$d[seq_len(k)]
  block_size_used <- min(as.integer(block_size), n, p)

  cpu_bench <- benchmark_elapsed(function() {
    irlba::irlba(host, nv = k, nu = k, work = max(3L * k, k + 20L))
  })
  compat_bench <- benchmark_elapsed(
    function() {
      am_irlba(x_mlx, nv = k, nu = k, work = max(3L * k, k + 20L))
    },
    warmup = function() {
      invisible(am_irlba(x_mlx, nv = k, nu = k, work = max(3L * k, k + 20L)))
    }
  )
  block_bench <- benchmark_elapsed(
    function() {
      am_block_lanczos(
        x_mlx,
        nv = k,
        nu = k,
        block_size = block_size,
        n_steps = n_steps
      )
    },
    warmup = function() {
      invisible(am_block_lanczos(
        x_mlx,
        nv = k,
        nu = k,
        block_size = block_size,
        n_steps = n_steps
      ))
    }
  )

  cpu <- cpu_bench$result
  compat <- compat_bench$result
  block <- block_bench$result

  data.frame(
    case = sprintf("%dx%d", n, p),
    k = k,
    block_size = block_size_used,
    n_steps = n_steps,
    implementation = c("irlba::irlba", "am_irlba", "am_block_lanczos"),
    elapsed = c(cpu_bench$elapsed, compat_bench$elapsed, block_bench$elapsed),
    max_rel_sv_err = c(
      max(abs(cpu$d - ref_d) / pmax(abs(ref_d), 1e-12)),
      max(abs(compat$d - ref_d) / pmax(abs(ref_d), 1e-12)),
      max(abs(block$d - ref_d) / pmax(abs(ref_d), 1e-12))
    ),
    iter = c(cpu$iter, compat$iter, block$iter),
    mprod = c(cpu$mprod, compat$mprod, block$mprod),
    stringsAsFactors = FALSE
  )
}

benchmark_block_lanczos_cases <- function(cases = NULL) {
  if (is.null(cases)) {
    cases <- list(
      list(n = 1200L, p = 600L, k = 10L, block_size = 11L, n_steps = 4L, seed = 20260406L),
      list(n = 3000L, p = 1200L, k = 20L, block_size = 20L, n_steps = 4L, seed = 20260407L)
    )
  }

  rows <- do.call(
    rbind,
    lapply(cases, function(case) {
      benchmark_case(
        n = case$n,
        p = case$p,
        k = case$k,
        block_size = case$block_size,
        n_steps = case$n_steps,
        seed = case$seed
      )
    })
  )

  rows$elapsed <- sprintf("%.3f", rows$elapsed)
  rows$max_rel_sv_err <- sprintf("%.4f", rows$max_rel_sv_err)
  rows
}

print_block_lanczos_benchmark <- function(cases = NULL) {
  rows <- benchmark_block_lanczos_cases(cases = cases)
  cat("Notes:\n")
  cat("- `am_irlba` is the compatibility wrapper around `irlba::irlba(..., fastpath = FALSE)`.\n")
  cat("- `am_block_lanczos` is the GEMM-oriented block Krylov surface; the default cases here use the current mildly oversampled block heuristic.\n")
  cat("- `max_rel_sv_err` is measured against `La.svd()` on the same host matrix.\n\n")
  print(rows, row.names = FALSE)
  invisible(rows)
}

if (isTRUE(getOption("amatrix.block_lanczos.benchmark.autorun", FALSE))) {
  print_block_lanczos_benchmark()
}
