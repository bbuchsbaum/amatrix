# Profile-guided dispatch calibration.
#
# GPU backends win only above a machine-specific matrix-size threshold. Below
# that threshold the upload + JIT overhead dominates and GPU is slower than CPU.
#
# amatrix_calibrate() benchmarks the requested backends, derives per-(op,
# backend) element-count thresholds, and stores them in the session and
# optionally on disk. amatrix_backend_plan() then uses these thresholds to
# gate cold-path GPU dispatch automatically.
#
# The resident path is never gated by calibration: if a matrix is already on
# the device the upload cost has been paid, so GPU is always preferred.

# ── Public API ────────────────────────────────────────────────────────────────

#' Calibrate GPU dispatch thresholds for this machine
#'
#' Runs micro-benchmarks for each (op, backend, size) combination and derives
#' the minimum matrix size at which each GPU backend reliably beats CPU. The
#' results are stored in the current session and optionally persisted to disk
#' for automatic reuse in future sessions.
#'
#' @param backend Character vector of backend names to benchmark. Defaults to
#'   all available non-CPU backends.
#' @param ops Character vector of operations to benchmark. Supported values:
#'   \code{"matmul"} (alias for \code{"gemm"}), \code{"gemm"},
#'   \code{"gemv"}, \code{"spmv"}, \code{"spmm"}, \code{"crossprod"},
#'   \code{"rowSums"}, \code{"colSums"}, \code{"qr"}, \code{"chol"},
#'   \code{"solve"}, \code{"svd"}.
#' @param sizes List of integer(2) vectors giving (nrow, ncol) test sizes.
#' @param sparse_densities Numeric vector of target sparse densities used when
#'   benchmarking \code{"spmv"} and \code{"spmm"}.
#' @param n_reps Number of timed repetitions per cell (after 2 warm-up reps).
#' @param margin Fraction by which GPU must beat CPU to count as a win (default
#'   0.10 = GPU must be at least 10\% faster).
#' @param persist Logical. Save calibration to the user cache directory so
#'   future sessions load it automatically.
#' @param quiet Logical. Suppress progress messages.
#' @return Invisibly, the calibration list (thresholds + full results table).
amatrix_calibrate <- function(
  backend  = NULL,
  ops      = c("gemm", "gemv", "crossprod", "rowSums", "colSums", "qr", "chol", "solve", "svd"),
  sizes    = list(c(64L, 32L), c(128L, 64L), c(256L, 128L), c(512L, 256L), c(1024L, 512L)),
  sparse_densities = c(0.01, 0.05, 0.20),
  n_reps   = 10L,
  margin   = 0.10,
  persist  = TRUE,
  quiet    = FALSE
) {
  stopifnot(
    is.character(ops), length(ops) >= 1L,
    is.list(sizes), length(sizes) >= 1L,
    is.numeric(n_reps), n_reps >= 1L,
    is.numeric(margin), margin >= 0, margin < 1
  )
  n_reps <- as.integer(n_reps)
  ops <- .amatrix_normalize_calibration_ops(ops)
  sparse_densities <- as.numeric(sparse_densities)
  stopifnot(length(sparse_densities) >= 1L, all(is.finite(sparse_densities)), all(sparse_densities > 0), all(sparse_densities <= 1))

  if (is.null(backend)) {
    backend <- setdiff(amatrix_backend_names(), "cpu")
  }
  backend <- as.character(backend)

  all_rows <- list()
  thresholds <- list()

  for (be in backend) {
    be_obj <- tryCatch(.amatrix_get_backend(be), error = function(e) NULL)
    if (is.null(be_obj) || !isTRUE(be_obj$available())) {
      if (!quiet) message(sprintf("amatrix_calibrate: '%s' not available, skipping", be))
      next
    }

    precision <- if ("fast" %in% be_obj$precision_modes()) "fast" else "strict"
    caps      <- be_obj$capabilities()
    be_rows   <- list()

    for (op in ops) {
      if (!.amatrix_backend_supports_calibration_op(be_obj, op)) next

      density_values <- if (op %in% c("spmv", "spmm")) sparse_densities else NA_real_
      for (density in density_values) {
      for (sz in sizes) {
        nr <- as.integer(sz[[1L]])
        nc <- as.integer(sz[[2L]])

        if (!quiet) message(sprintf(
          "  calibrating %s / %s / %dx%d%s ...", be, op, nr, nc,
          if (is.na(density)) "" else sprintf(" @ density %.3f", density)
        ), appendLF = FALSE)

        row <- tryCatch(
          .amatrix_benchmark_op(be, be_obj, op, nr, nc, precision, n_reps, sparse_density = density),
          error = function(e) NULL
        )

        if (is.null(row)) {
          if (!quiet) message(" failed, skipped")
          next
        }

        row$margin <- margin
        row$gpu_wins <- row$gpu_ms < row$cpu_ms * (1 - margin)

        if (!quiet) message(sprintf(
          " cpu=%.1fms  gpu=%.1fms  speedup=%.2fx  %s",
          row$cpu_ms, row$gpu_ms,
          row$cpu_ms / row$gpu_ms,
          if (row$gpu_wins) "GPU wins" else "CPU wins"
        ))

        be_rows[[length(be_rows) + 1L]] <- row
      }
      }
    }

    if (length(be_rows) == 0L) next

    results_be <- do.call(rbind, lapply(be_rows, as.data.frame, stringsAsFactors = FALSE))
    all_rows   <- c(all_rows, be_rows)
    thresholds[[be]] <- .amatrix_derive_thresholds(results_be, unique(results_be$op))
  }

  results_df <- if (length(all_rows) > 0L) {
    do.call(rbind, lapply(all_rows, as.data.frame, stringsAsFactors = FALSE))
  } else {
    data.frame(
      backend = character(), op = character(), op_base = character(),
      nrow = integer(), ncol = integer(), elements = integer(),
      nnz = integer(), density = numeric(), density_bucket = character(),
      cpu_ms = numeric(), gpu_ms = numeric(),
      margin = numeric(), gpu_wins = logical(),
      stringsAsFactors = FALSE
    )
  }

  calibration <- list(
    version       = "1",
    calibrated_at = Sys.time(),
    thresholds    = thresholds,
    results       = results_df
  )

  .amatrix_state$calibration <- calibration

  if (persist) {
    path <- .amatrix_calibration_path()
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tryCatch(
      saveRDS(calibration, path),
      error = function(e) {
        if (!quiet) message(sprintf("amatrix_calibrate: could not persist to %s: %s", path, e$message))
      }
    )
    if (!quiet) message(sprintf("amatrix_calibrate: saved to %s", path))
  }

  if (!quiet) {
    message("\nThresholds (min elements for GPU dispatch on cold path):")
    for (be in names(thresholds)) {
      for (op in names(thresholds[[be]])) {
        thresh <- thresholds[[be]][[op]]
        message(sprintf("  %-12s %-10s  %s",
          be, op,
          if (is.infinite(thresh)) "never (CPU always faster)"
          else if (thresh == 0L)   "always (GPU always faster)"
          else                     sprintf(">= %d elements", thresh)
        ))
      }
    }
  }

  invisible(calibration)
}

#' Inspect the current calibration
#'
#' Returns the calibration list stored in the session (loaded from disk if not
#' yet loaded). Returns \code{NULL} if no calibration is available.
amatrix_calibration_info <- function() {
  .amatrix_load_calibration()
  .amatrix_state$calibration
}

# ── Internals ─────────────────────────────────────────────────────────────────

.amatrix_calibration_path <- function() {
  cache_dir <- tryCatch(
    tools::R_user_dir("amatrix", "cache"),
    error = function(e) file.path(Sys.getenv("HOME"), ".amatrix", "cache")
  )
  file.path(cache_dir, "calibration.rds")
}

# Lazy-load calibration from disk on first use.
.amatrix_load_calibration <- function() {
  if (!is.null(.amatrix_state$calibration)) return(invisible(NULL))
  path <- .amatrix_calibration_path()
  if (!file.exists(path)) return(invisible(NULL))
  cal <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.null(cal) && identical(cal$version, "1")) {
    .amatrix_state$calibration <- cal
  }
  invisible(NULL)
}

# Returns TRUE when the workload clears the calibrated threshold for
# (op class, backend). Always returns TRUE when no calibration data exists
# (backward compatible) or for the CPU backend.
.amatrix_calibration_ok <- function(x, op, backend_name, y = NULL) {
  if (backend_name == "cpu") return(TRUE)
  .amatrix_load_calibration()
  cal <- .amatrix_state$calibration
  if (is.null(cal)) return(TRUE)

  op_key <- .amatrix_dispatch_signature(x, op, y = y)
  thresh <- cal$thresholds[[backend_name]][[op_key]]
  if (is.null(thresh)) return(TRUE)     # op not benchmarked
  if (is.infinite(thresh)) return(FALSE) # GPU never wins

  .amatrix_dispatch_workload(x, op, y = y) >= thresh
}

# Benchmark one (backend, op, size) cell. Returns a one-row list.
.amatrix_benchmark_op <- function(be_name, be_obj, op, nr, nc, precision, n_reps, sparse_density = NA_real_) {
  X_host  <- matrix(seq_len(nr * nc) / (nr * nc + 1L), nr, nc)
  Y_host <- matrix(seq_len(nc * nc) / (nc * nc + 1L), nc, nc)
  v_host <- matrix(seq_len(nc) / (nc + 1L), nc, 1L)
  spd_host <- crossprod(X_host) + diag(nc)  # nc×nc SPD for chol
  solve_rhs_width <- min(max(2L, getOption("amatrix.calibration.solve_rhs_width", 8L)), max(2L, nc))
  solve_rhs_host <- matrix(seq_len(nc * solve_rhs_width) / (nc * solve_rhs_width + 1L), nc, solve_rhs_width)
  svd_rank <- min(nr, nc, max(1L, as.integer(getOption("amatrix.calibration.svd_rank", 16L))))

  X_adge  <- as_adgeMatrix(X_host,  preferred_backend = be_name, precision = precision)
  Y_adge <- as_adgeMatrix(Y_host, preferred_backend = be_name, precision = precision)
  v_adge <- as_adgeMatrix(v_host, preferred_backend = be_name, precision = precision)
  SPD_adge <- as_adgeMatrix(spd_host, preferred_backend = be_name, precision = precision)
  solve_rhs_adge <- as_adgeMatrix(solve_rhs_host, preferred_backend = be_name, precision = precision)

  if (op %in% c("spmv", "spmm")) {
    X_sparse_host <- .amatrix_sparse_benchmark_matrix(nr, nc, density = sparse_density)
    rhs_width <- if (identical(op, "spmv")) 1L else min(max(2L, getOption("amatrix.calibration.spmm_rhs_width", 8L)), max(2L, nc))
    rhs_host <- matrix(seq_len(nc * rhs_width) / (nc * rhs_width + 1L), nc, rhs_width)
    X_sparse <- as_adgCMatrix(X_sparse_host, preferred_backend = be_name, precision = precision)
    rhs_adge <- as_adgeMatrix(rhs_host, preferred_backend = be_name, precision = precision)

    if (!isTRUE(be_obj$supports("matmul", X_sparse, y = rhs_adge))) {
      return(NULL)
    }

    cpu_fn <- function() X_sparse_host %*% rhs_host
    gpu_fn <- function() be_obj$matmul(X_sparse, rhs_adge)
    nnz <- length(X_sparse_host@x)
    density_value <- .amatrix_sparse_density(X_sparse_host)
    density_bucket <- .amatrix_sparse_density_bucket(density_value)
    op_label <- paste0(op, ":", density_bucket)
    workload <- as.integer(nnz * rhs_width)
  } else {
    probe_y <- switch(op,
      gemm = Y_adge,
      gemv = v_adge,
      solve = solve_rhs_adge,
      NULL
    )
    if (!isTRUE(.amatrix_backend_supports_benchmark_cell(be_obj, op, X_adge, probe_y, SPD_adge))) {
      return(NULL)
    }

    cpu_fn <- switch(op,
      gemm      = function() X_host %*% Y_host,
      gemv      = function() X_host %*% v_host,
      crossprod = function() base::crossprod(X_host),
      rowSums   = function() base::rowSums(X_host),
      colSums   = function() base::colSums(X_host),
      qr        = function() base::qr(X_host),
      chol      = function() base::chol(spd_host),
      solve     = function() base::solve(spd_host, solve_rhs_host),
      svd       = function() base::svd(X_host, nu = svd_rank, nv = svd_rank)
    )
    gpu_fn <- switch(op,
      gemm      = function() be_obj$matmul(X_adge, Y_adge),
      gemv      = function() be_obj$matmul(X_adge, v_adge),
      crossprod = function() be_obj$crossprod(X_adge),
      rowSums   = function() be_obj$rowSums(X_adge),
      colSums   = function() be_obj$colSums(X_adge),
      qr        = function() be_obj$qr(X_adge),
      chol      = function() be_obj$chol(SPD_adge),
      solve     = function() be_obj$solve(SPD_adge, solve_rhs_adge),
      svd       = function() be_obj$svd(X_adge, nu = svd_rank, nv = svd_rank)
    )
    nnz <- NA_integer_
    density_value <- NA_real_
    density_bucket <- NA_character_
    op_label <- op
    workload <- .amatrix_benchmark_workload(op, nr, nc)
  }

  # Two warm-up reps each (discard compilation cost for GPU, cache for CPU)
  for (i in seq_len(2L)) tryCatch(cpu_fn(), error = function(e) NULL)
  for (i in seq_len(2L)) tryCatch(gpu_fn(), error = function(e) NULL)

  .time_ms <- function(fn, reps) {
    times <- vapply(seq_len(reps), function(i) {
      t0 <- proc.time()[["elapsed"]]
      tryCatch(fn(), error = function(e) NULL)
      (proc.time()[["elapsed"]] - t0) * 1000
    }, numeric(1))
    median(times)
  }

  cpu_ms <- .time_ms(cpu_fn, n_reps)
  gpu_ms <- .time_ms(gpu_fn, n_reps)

  list(
    backend  = be_name,
    op       = op_label,
    op_base  = op,
    nrow     = nr,
    ncol     = nc,
    nnz      = nnz,
    density  = density_value,
    density_bucket = density_bucket,
    elements = workload,
    cpu_ms   = cpu_ms,
    gpu_ms   = gpu_ms
  )
}

# For each op, find the smallest element count where gpu_wins is TRUE.
# Returns Inf if GPU never wins, 0L if GPU wins at every tested size.
.amatrix_derive_thresholds <- function(results, ops) {
  thresholds <- list()
  for (op in ops) {
    sub <- results[results$op == op, , drop = FALSE]
    if (nrow(sub) == 0L) next
    sub <- sub[order(sub$elements), ]
    winning <- sub$elements[sub$gpu_wins]
    thresholds[[op]] <- if (length(winning) == 0L) Inf else min(winning)
  }
  thresholds
}

.amatrix_normalize_calibration_ops <- function(ops) {
  stopifnot(is.character(ops), length(ops) >= 1L)
  ops <- vapply(ops, function(op) {
    if (identical(op, "matmul")) "gemm" else op
  }, character(1))

  allowed <- c("gemm", "gemv", "crossprod", "rowSums", "colSums", "qr", "chol", "solve", "svd", "spmv", "spmm")
  invalid <- setdiff(ops, allowed)
  if (length(invalid) > 0L) {
    stop(sprintf(
      "unsupported calibration ops: %s",
      paste(invalid, collapse = ", ")
    ), call. = FALSE)
  }

  unique(ops)
}

.amatrix_backend_supports_calibration_op <- function(backend, op) {
  caps <- backend$capabilities()
  if (op %in% c("gemm", "gemv", "spmv", "spmm")) {
    return("matmul" %in% caps)
  }
  op %in% caps
}

.amatrix_backend_supports_benchmark_cell <- function(backend, op, X_adge, y_arg = NULL, SPD_adge = NULL) {
  mapped_op <- if (op %in% c("gemm", "gemv")) "matmul" else op
  x_arg <- if (op %in% c("chol", "solve")) SPD_adge else X_adge
  isTRUE(backend$supports(mapped_op, x_arg, y = y_arg))
}

.amatrix_rhs_dims <- function(y) {
  if (is.null(y)) {
    return(NULL)
  }

  dims <- dim(y)
  if (is.null(dims)) {
    return(c(length(y), 1L))
  }

  if (length(dims) == 1L) {
    return(c(dims[[1L]], 1L))
  }

  dims
}

.amatrix_rhs_width <- function(y) {
  dims <- .amatrix_rhs_dims(y)
  if (is.null(dims)) {
    return(NA_integer_)
  }
  as.integer(dims[[2L]])
}

.amatrix_dispatch_signature <- function(x, op, y = NULL) {
  if (!identical(op, "matmul")) {
    return(op)
  }

  width <- .amatrix_rhs_width(y)
  if (inherits(x, "adgCMatrix")) {
    bucket <- .amatrix_sparse_density_bucket(x)
    if (!is.na(width) && width <= 1L) {
      return(paste0("spmv:", bucket))
    }
    return(paste0("spmm:", bucket))
  }

  if (!is.na(width) && width <= 1L) {
    return("gemv")
  }

  "gemm"
}

.amatrix_dispatch_workload <- function(x, op, y = NULL) {
  sig <- .amatrix_dispatch_signature(x, op, y = y)
  width <- .amatrix_rhs_width(y)
  width_eff <- if (is.na(width)) 1L else max(1L, width)

  if (inherits(x, "adgCMatrix")) {
    nnz <- length(x@x)
    return(as.integer(nnz * if (startsWith(sig, "spmm:")) width_eff else 1L))
  }

  base_work <- as.integer(nrow(x) * ncol(x))
  if (identical(sig, "solve")) {
    return(as.integer(base_work * width_eff))
  }
  if (identical(sig, "svd")) {
    return(as.integer(base_work * min(nrow(x), ncol(x))))
  }
  if (sig %in% c("gemm")) {
    return(as.integer(base_work * width_eff))
  }

  base_work
}

.amatrix_benchmark_workload <- function(op, nr, nc) {
  switch(op,
    gemm = as.integer(nr * nc * nc),
    gemv = as.integer(nr * nc),
    crossprod = as.integer(nr * nc),
    rowSums = as.integer(nr * nc),
    colSums = as.integer(nr * nc),
    qr = as.integer(nr * nc),
    chol = as.integer(nc * nc),
    solve = as.integer(nc * nc * min(max(2L, getOption("amatrix.calibration.solve_rhs_width", 8L)), max(2L, nc))),
    svd = as.integer(nr * nc * min(nr, nc)),
    spmv = as.integer(nr * nc),
    spmm = as.integer(nr * nc * nc),
    as.integer(nr * nc)
  )
}

.amatrix_sparse_density <- function(x) {
  stopifnot(inherits(x, "sparseMatrix"))
  total <- nrow(x) * ncol(x)
  if (total <= 0L) {
    return(0)
  }
  length(x@x) / total
}

.amatrix_sparse_density_bucket <- function(x) {
  density <- if (inherits(x, "sparseMatrix")) .amatrix_sparse_density(x) else as.double(x)

  if (!is.finite(density) || density <= 0.01) {
    return("ultra_sparse")
  }
  if (density <= 0.05) {
    return("sparse")
  }
  if (density <= 0.20) {
    return("semi_dense")
  }
  "denseish"
}

.amatrix_sparse_benchmark_matrix <- function(nr, nc, density = 0.05) {
  stopifnot(nr >= 1L, nc >= 1L, is.finite(density), density > 0, density <= 1)

  total <- nr * nc
  target <- max(1L, min(total, as.integer(ceiling(total * density))))
  positions <- sample.int(total, target, replace = FALSE)
  rows <- ((positions - 1L) %% nr) + 1L
  cols <- ((positions - 1L) %/% nr) + 1L

  missing_cols <- setdiff(seq_len(nc), unique(cols))
  if (length(missing_cols) > 0L) {
    rows <- c(rows, sample.int(nr, length(missing_cols), replace = TRUE))
    cols <- c(cols, missing_cols)
  }

  vals <- seq_along(rows) / (length(rows) + 1L)
  Matrix::sparseMatrix(
    i = rows,
    j = cols,
    x = vals,
    dims = c(nr, nc),
    giveCsparse = TRUE
  )
}
