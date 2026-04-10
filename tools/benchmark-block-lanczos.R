#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else if (requireNamespace("amatrix", quietly = TRUE)) {
    library(amatrix)
  } else {
    stop("Either the installed 'amatrix' package or 'pkgload' is required", call. = FALSE)
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for sparse block Lanczos benchmarks", call. = FALSE)
  }
  if (!requireNamespace("irlba", quietly = TRUE)) {
    stop("Package 'irlba' is required for this benchmark", call. = FALSE)
  }
})

source(file.path("tools", "benchmark-helpers.R"), local = TRUE)

available_spectral_backends <- function() {
  available_benchmark_backends(include_cpu = FALSE, include_mlx = TRUE, include_metal = TRUE, include_opencl = TRUE)
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

benchmark_dense_case <- function(n, p, k, block_size, n_steps, seed) {
  set.seed(seed)
  host <- matrix(rnorm(n * p), nrow = n, ncol = p)
  ref_d <- La.svd(host, nu = k, nv = k)$d[seq_len(k)]
  block_size_used <- min(as.integer(block_size), n, p)
  gpu_backends <- available_spectral_backends()
  cpu_x <- adgeMatrix(host, preferred_backend = "cpu", precision = "strict")

  rows <- list()

  cpu_native <- benchmark_elapsed(function() {
    irlba::irlba(host, nv = k, nu = k, work = max(3L * k, k + 20L))
  })
  cpu_compat <- benchmark_elapsed(function() {
    irlba(cpu_x, nv = k, nu = k, work = max(3L * k, k + 20L))
  })
  cpu_block <- benchmark_elapsed(function() {
    block_lanczos(cpu_x, nv = k, nu = k, block_size = block_size, n_steps = n_steps)
  })

  cpu_results <- list(
    native = cpu_native$result,
    compat = cpu_compat$result,
    block = cpu_block$result
  )

  rows[[length(rows) + 1L]] <- data.frame(
    case = sprintf("%dx%d", n, p),
    backend = "cpu",
    precision = "strict",
    k = k,
    block_size = block_size_used,
    n_steps = n_steps,
    implementation = c("irlba::irlba", "irlba", "block_lanczos"),
    elapsed = c(cpu_native$elapsed, cpu_compat$elapsed, cpu_block$elapsed),
    max_rel_sv_err = c(
      max(abs(cpu_results$native$d - ref_d) / pmax(abs(ref_d), 1e-12)),
      max(abs(cpu_results$compat$d - ref_d) / pmax(abs(ref_d), 1e-12)),
      max(abs(cpu_results$block$d - ref_d) / pmax(abs(ref_d), 1e-12))
    ),
    iter = c(cpu_results$native$iter, cpu_results$compat$iter, cpu_results$block$iter),
    mprod = c(cpu_results$native$mprod, cpu_results$compat$mprod, cpu_results$block$mprod),
    stringsAsFactors = FALSE
  )

  for (backend_spec in gpu_backends) {
    x_gpu <- adgeMatrix(host, preferred_backend = backend_spec$name, precision = backend_spec$precision)
    compat_bench <- benchmark_elapsed(
      function() {
        irlba(x_gpu, nv = k, nu = k, work = max(3L * k, k + 20L))
      },
      warmup = function() {
        invisible(irlba(x_gpu, nv = k, nu = k, work = max(3L * k, k + 20L)))
      }
    )
    block_bench <- benchmark_elapsed(
      function() {
        block_lanczos(
          x_gpu,
          nv = k,
          nu = k,
          block_size = block_size,
          n_steps = n_steps
        )
      },
      warmup = function() {
        invisible(block_lanczos(
          x_gpu,
          nv = k,
          nu = k,
          block_size = block_size,
          n_steps = n_steps
        ))
      }
    )

    compat <- compat_bench$result
    block <- block_bench$result

    rows[[length(rows) + 1L]] <- data.frame(
      case = sprintf("%dx%d", n, p),
      backend = backend_spec$name,
      precision = backend_spec$precision,
      k = k,
      block_size = block_size_used,
      n_steps = n_steps,
      implementation = c("irlba", "block_lanczos"),
      elapsed = c(compat_bench$elapsed, block_bench$elapsed),
      max_rel_sv_err = c(
        max(abs(compat$d - ref_d) / pmax(abs(ref_d), 1e-12)),
        max(abs(block$d - ref_d) / pmax(abs(ref_d), 1e-12))
      ),
      iter = c(compat$iter, block$iter),
      mprod = c(compat$mprod, block$mprod),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

benchmark_sparse_case <- function(n, p, k, block_size, n_steps, density, seed) {
  set.seed(seed)
  host <- Matrix::rsparsematrix(n, p, density = density)
  ref_d <- La.svd(as.matrix(host), nu = k, nv = k)$d[seq_len(k)]
  block_size_used <- min(as.integer(block_size), n, p)
  gpu_backends <- available_spectral_backends()
  cpu_x <- adgCMatrix(host, preferred_backend = "cpu", precision = "strict")

  rows <- list()

  cpu_block <- benchmark_elapsed(function() {
    block_lanczos(cpu_x, nv = k, nu = k, block_size = block_size, n_steps = n_steps)
  })

  rows[[length(rows) + 1L]] <- data.frame(
    case = sprintf("%dx%d", n, p),
    matrix_type = "sparse",
    density = density,
    backend = "cpu",
    precision = "strict",
    k = k,
    block_size = block_size_used,
    n_steps = n_steps,
    implementation = "block_lanczos",
    elapsed = cpu_block$elapsed,
    max_rel_sv_err = max(abs(cpu_block$result$d - ref_d) / pmax(abs(ref_d), 1e-12)),
    iter = cpu_block$result$iter,
    mprod = cpu_block$result$mprod,
    stringsAsFactors = FALSE
  )

  for (backend_spec in gpu_backends) {
    x_gpu <- adgCMatrix(host, preferred_backend = backend_spec$name, precision = backend_spec$precision)
    block_bench <- benchmark_elapsed(
      function() {
        block_lanczos(
          x_gpu,
          nv = k,
          nu = k,
          block_size = block_size,
          n_steps = n_steps
        )
      },
      warmup = function() {
        invisible(block_lanczos(
          x_gpu,
          nv = k,
          nu = k,
          block_size = block_size,
          n_steps = n_steps
        ))
      }
    )

    rows[[length(rows) + 1L]] <- data.frame(
      case = sprintf("%dx%d", n, p),
      matrix_type = "sparse",
      density = density,
      backend = backend_spec$name,
      precision = backend_spec$precision,
      k = k,
      block_size = block_size_used,
      n_steps = n_steps,
      implementation = "block_lanczos",
      elapsed = block_bench$elapsed,
      max_rel_sv_err = max(abs(block_bench$result$d - ref_d) / pmax(abs(ref_d), 1e-12)),
      iter = block_bench$result$iter,
      mprod = block_bench$result$mprod,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

benchmark_case <- function(n, p, k, block_size, n_steps, seed, matrix_type = c("dense", "sparse"), density = 0.05) {
  matrix_type <- match.arg(matrix_type)

  if (identical(matrix_type, "sparse")) {
    return(benchmark_sparse_case(
      n = n,
      p = p,
      k = k,
      block_size = block_size,
      n_steps = n_steps,
      density = density,
      seed = seed
    ))
  }

  rows <- benchmark_dense_case(
    n = n,
    p = p,
    k = k,
    block_size = block_size,
    n_steps = n_steps,
    seed = seed
  )
  rows$matrix_type <- "dense"
  rows$density <- NA_real_
  rows[, c("case", "matrix_type", "density", "backend", "precision", "k", "block_size", "n_steps", "implementation", "elapsed", "max_rel_sv_err", "iter", "mprod")]
}

benchmark_block_lanczos_cases <- function(cases = NULL) {
  if (is.null(cases)) {
    cases <- list(
      list(n = 1200L, p = 600L, k = 10L, block_size = 11L, n_steps = 4L, seed = 20260406L, matrix_type = "dense"),
      list(n = 3000L, p = 1200L, k = 20L, block_size = 20L, n_steps = 4L, seed = 20260407L, matrix_type = "dense"),
      list(n = 4000L, p = 1000L, k = 8L, block_size = 16L, n_steps = 4L, density = 0.05, seed = 20260409L, matrix_type = "sparse")
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
        seed = case$seed,
        matrix_type = if (is.null(case$matrix_type)) "dense" else case$matrix_type,
        density = if (is.null(case$density)) 0.05 else case$density
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
  cat("- `matrix_type` distinguishes the legacy dense benchmark path from the new sparse block-Krylov path.\n")
  cat("- `backend` indicates where the amatrix wrapper/block path was run; `irlba::irlba` is the host baseline.\n")
  cat("- `irlba()` is the compatibility wrapper around `irlba::irlba(..., fastpath = FALSE)`.\n")
  cat("- `block_lanczos()` is the GEMM-oriented block Krylov surface; the default cases here use the current mildly oversampled block heuristic.\n")
  cat("- Sparse cases benchmark `block_lanczos()` directly, because that is the path that now uses compiled repeated-product plans.\n")
  cat("- `max_rel_sv_err` is measured against `La.svd()` on the same host matrix.\n\n")
  print(rows, row.names = FALSE)
  invisible(rows)
}

if (isTRUE(getOption("amatrix.block_lanczos.benchmark.autorun", FALSE))) {
  print_block_lanczos_benchmark()
}
