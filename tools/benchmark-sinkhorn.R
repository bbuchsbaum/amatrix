#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(amatrix)
  }
})

source(file.path("tools", "benchmark-helpers.R"), local = TRUE)

available_sinkhorn_backends <- function() {
  available_benchmark_backends(include_cpu = TRUE, include_mlx = TRUE, include_opencl = TRUE)
}

make_positive_matrix <- function(n, seed = 1L) {
  set.seed(seed)
  x <- exp(matrix(rnorm(n * n), nrow = n, ncol = n))
  storage.mode(x) <- "double"
  x
}

sinkhorn_base <- function(A, max_iter = 50L, tol = 0, check_every = 5L, eps = 1e-15) {
  work <- A
  converged <- FALSE
  row_error <- Inf
  col_error <- Inf

  for (iter in seq_len(max_iter)) {
    work <- work / pmax(base::rowSums(work), eps)
    work <- t(t(work) / pmax(base::colSums(work), eps))

    if ((iter %% check_every) == 0L || iter == max_iter) {
      row_error <- max(abs(base::rowSums(work) - 1))
      col_error <- max(abs(base::colSums(work) - 1))
      if (max(row_error, col_error) < tol) {
        converged <- TRUE
        break
      }
    }
  }

  list(
    matrix = work,
    iterations = iter,
    converged = converged,
    row_error = row_error,
    col_error = col_error
  )
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

benchmark_sinkhorn_case <- function(n, fixed_iters = 50L, reps = 3L, seed = 1L) {
  A <- make_positive_matrix(n, seed = seed)
  rows <- list()

  base_bench <- benchmark_elapsed(function() {
    sinkhorn_base(A, max_iter = fixed_iters, tol = 0)
  }, reps = reps)

  rows[[length(rows) + 1L]] <- data.frame(
    case = sprintf("%dx%d", n, n),
    backend = "base",
    precision = "double",
    method = "host",
    elapsed = base_bench$elapsed,
    iterations = base_bench$result$iterations,
    row_error = base_bench$result$row_error,
    col_error = base_bench$result$col_error,
    stringsAsFactors = FALSE
  )

  for (backend_spec in available_sinkhorn_backends()) {
    backend_name <- backend_spec$name
    precision <- backend_spec$precision

    if (identical(backend_name, "cpu")) {
      bench <- benchmark_elapsed(function() {
        sinkhorn(A, max_iter = fixed_iters, tol = 0, mode = "fast", return_info = TRUE)
      }, reps = reps)
    } else {
      x_gpu <- adgeMatrix(A, preferred_backend = backend_name, precision = precision)
      bench <- benchmark_elapsed(
        function() {
          sinkhorn(x_gpu, max_iter = fixed_iters, tol = 0, return_info = TRUE)
        },
        reps = reps,
        warmup = function() {
          invisible(sinkhorn(x_gpu, max_iter = 5L, tol = 0, return_info = TRUE))
        }
      )
    }

    rows[[length(rows) + 1L]] <- data.frame(
      case = sprintf("%dx%d", n, n),
      backend = backend_name,
      precision = precision,
      method = bench$result$method,
      elapsed = bench$elapsed,
      iterations = bench$result$iterations,
      row_error = bench$result$row_error,
      col_error = bench$result$col_error,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  out$elapsed <- sprintf("%.3f", out$elapsed)
  out$row_error <- formatC(out$row_error, format = "e", digits = 2L)
  out$col_error <- formatC(out$col_error, format = "e", digits = 2L)
  out
}

print_sinkhorn_benchmark <- function(cases = c(200L, 500L, 1000L), fixed_iters = 50L, reps = 3L) {
  rows <- do.call(
    rbind,
    lapply(seq_along(cases), function(idx) {
      benchmark_sinkhorn_case(
        n = cases[[idx]],
        fixed_iters = fixed_iters,
        reps = reps,
        seed = 1000L + idx
      )
    })
  )

  cat("Notes:\n")
  cat("- `base` is a plain matrix reference loop.\n")
  cat("- `cpu` is the exported amatrix surface on a host matrix.\n")
  cat("- GPU rows use fast adgeMatrix inputs and report whether the algorithm stayed on the resident path.\n")
  cat("- Set `AMATRIX_BENCHMARK_ARRAYFIRE=1` or `AMATRIX_ARRAYFIRE_PROBE_GPU=1` to include ArrayFire in this sweep.\n\n")
  print(rows, row.names = FALSE)
  invisible(rows)
}

if (sys.nframe() == 0L) {
  print_sinkhorn_benchmark()
}
