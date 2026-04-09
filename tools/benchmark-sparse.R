#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(amatrix)
  }
  library(Matrix)
})

source(file.path("tools", "benchmark-helpers.R"), local = TRUE)

available_sparse_backends <- function() {
  benchmark_backend_names(include_cpu = TRUE, include_mlx = TRUE, include_metal = TRUE, include_opencl = FALSE)
}

.time_ms <- function(fn, reps = 5L, warmup = 1L, batch_reps = 10L) {
  for (i in seq_len(warmup)) {
    tryCatch(fn(), error = function(e) NULL)
  }

  times <- vapply(seq_len(reps), function(i) {
    t0 <- unclass(Sys.time())
    for (j in seq_len(batch_reps)) {
      tryCatch(fn(), error = function(e) stop(e), silent = FALSE)
    }
    ((unclass(Sys.time()) - t0) * 1000) / batch_reps
  }, numeric(1))

  median(times)
}

make_sparse_case <- function(nrow, ncol, density, rhs_width, seed) {
  set.seed(seed)
  X_host <- rsparsematrix(nrow, ncol, density = density)
  rhs_host <- matrix(rnorm(ncol * rhs_width), nrow = ncol, ncol = rhs_width)
  list(X_host = X_host, rhs_host = rhs_host)
}

make_sparse_operand <- function(x_host, backend) {
  precision <- if (identical(backend, "cpu")) "strict" else "fast"
  as_adgCMatrix(x_host, preferred_backend = backend, precision = precision)
}

.release_residency <- function(x) {
  if (!inherits(x, "aMatrix")) {
    return(invisible(FALSE))
  }

  entry <- amatrix:::.amatrix_resident_entry(x)
  if (is.null(entry)) {
    return(invisible(FALSE))
  }

  backend <- tryCatch(amatrix:::.amatrix_get_backend(entry$backend), error = function(e) NULL)
  if (is.null(backend)) {
    return(invisible(FALSE))
  }

  if (isTRUE(entry$sparse) && is.function(backend$sparse_resident_drop)) {
    try(backend$sparse_resident_drop(entry$resident_key), silent = TRUE)
  } else if (is.function(backend$resident_drop)) {
    try(backend$resident_drop(entry$resident_key), silent = TRUE)
  }
  amatrix:::.amatrix_drop_resident_binding(x)
  invisible(TRUE)
}

benchmark_sparse_cell <- function(backend, nrow, ncol, density, rhs_width, reps = 5L) {
  case <- make_sparse_case(nrow = nrow, ncol = ncol, density = density,
                           rhs_width = rhs_width, seed = 1000L + rhs_width + as.integer(density * 1000))
  op_label <- if (rhs_width == 1L) "spmv" else "spmm"

  x_plan <- make_sparse_operand(case$X_host, backend)
  plan <- amatrix_backend_plan(x_plan, "matmul", y = case$rhs_host)
  chosen <- plan$chosen

  cold_ms <- .time_ms(function() {
    x <- make_sparse_operand(case$X_host, backend)
    on.exit(.release_residency(x), add = TRUE)
    invisible(x %*% case$rhs_host)
  }, reps = reps)

  resident_x <- make_sparse_operand(case$X_host, backend)
  on.exit(.release_residency(resident_x), add = TRUE)
  invisible(resident_x %*% case$rhs_host)
  resident_ms <- .time_ms(function() {
    invisible(resident_x %*% case$rhs_host)
  }, reps = reps)

  cpu_ref_ms <- .time_ms(function() {
    invisible(case$X_host %*% case$rhs_host)
  }, reps = reps)

  data.frame(
    backend = backend,
    chosen = chosen,
    op = op_label,
    nrow = nrow,
    ncol = ncol,
    nnz = length(case$X_host@x),
    density = density,
    rhs_width = rhs_width,
    cold_ms = cold_ms,
    resident_ms = resident_ms,
    cpu_ref_ms = cpu_ref_ms,
    cold_speedup_vs_cpu = cpu_ref_ms / cold_ms,
    resident_speedup_vs_cpu = cpu_ref_ms / resident_ms,
    stringsAsFactors = FALSE
  )
}

run_sparse_benchmarks <- function(
  backends = available_sparse_backends(),
  sizes = list(c(4000L, 1000L), c(8000L, 2000L)),
  densities = c(0.001, 0.005, 0.01, 0.05),
  rhs_widths = c(1L, 8L, 32L),
  reps = 5L
) {
  rows <- list()

  for (backend in backends) {
    for (size in sizes) {
      nr <- size[[1L]]
      nc <- size[[2L]]

      for (density in densities) {
        for (rhs_width in rhs_widths) {
          rows[[length(rows) + 1L]] <- benchmark_sparse_cell(
            backend = backend,
            nrow = nr,
            ncol = nc,
            density = density,
            rhs_width = rhs_width,
            reps = reps
          )
        }
      }
    }
  }

  do.call(rbind, rows)
}

print_sparse_benchmark_summary <- function(results) {
  cat("=== Sparse Product Benchmark ===\n")
  print(results, row.names = FALSE)

  cat("\n=== Cases where non-CPU backend wins vs CPU reference ===\n")
  wins <- subset(results, backend != "cpu" & chosen == backend & resident_speedup_vs_cpu > 1)
  if (nrow(wins) == 0L) {
    cat("No non-CPU sparse wins recorded in this run.\n")
  } else {
    print(wins, row.names = FALSE)
  }
}

if (!identical(Sys.getenv("AMATRIX_BENCHMARK_NO_AUTORUN", unset = ""), "1")) {
  results <- run_sparse_benchmarks()
  print_sparse_benchmark_summary(results)
}
