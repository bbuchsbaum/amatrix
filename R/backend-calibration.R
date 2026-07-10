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
#' Runs micro-benchmarks for each (op, backend, size) combination and
#' derives the minimum matrix element count at which each GPU backend
#' reliably outperforms CPU. Results are stored in the current session
#' and optionally persisted to disk for reuse in future sessions.
#'
#' @param backend Character vector of backend names to benchmark.
#'   Defaults to all registered non-CPU backends that report
#'   \code{available = TRUE}.
#' @param ops Character vector of operations to benchmark. Supported
#'   values: \code{"matmul"} (alias for \code{"gemm"}), \code{"gemm"},
#'   \code{"gemv"}, \code{"spmv"}, \code{"spmm"}, \code{"crossprod"},
#'   \code{"rowSums"}, \code{"colSums"}, \code{"qr"}, \code{"chol"},
#'   \code{"solve"}, \code{"svd"}.
#' @param sizes List of integer vectors of length 2 giving
#'   \code{c(nrow, ncol)} test matrix dimensions.
#' @param sparse_densities Numeric vector of target fill densities used
#'   when benchmarking \code{"spmv"} and \code{"spmm"}.
#' @param n_reps Positive integer. Number of timed repetitions per
#'   benchmark cell, after 2 warm-up repetitions.
#' @param margin Non-negative numeric less than 1. Fraction by which
#'   GPU median time must beat CPU to count as a GPU win (default
#'   \code{0.10} means GPU must be at least 10\% faster).
#' @param persist Logical. If \code{TRUE} (default), save calibration
#'   to the user cache directory so future sessions load it
#'   automatically.
#' @param quiet Logical. Suppress progress messages.
#'
#' @return Invisibly, a list with elements \code{version},
#'   \code{calibrated_at} (POSIXct), \code{thresholds} (nested list
#'   keyed by backend then op), and \code{results} (data.frame of all
#'   benchmark measurements).
#'
#' @seealso \code{\link{amatrix_calibration_info}},
#'   \code{\link{amatrix_backend_plan}}
#' @export
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
    if (is.null(be_obj) || !.amatrix_backend_available_safe(be_obj)) {
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
    version       = "2",
    calibrated_at = Sys.time(),
    sys_hash      = .amatrix_sys_hash(),
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

#' Retrieve the current calibration state
#'
#' Returns the calibration data stored in the current session. If no
#' calibration has been run yet, the function attempts to load a
#' previously persisted calibration from the user cache directory.
#' Returns \code{NULL} when no calibration is available.
#'
#' @return A list as returned by \code{\link{amatrix_calibrate}}, or
#'   \code{NULL} if no calibration data exists for this session.
#'
#' @examples
#' cal <- amatrix_calibration_info()
#' is.null(cal) # TRUE when no calibration has been run
#'
#' @seealso \code{\link{amatrix_calibrate}}
#' @export
amatrix_calibration_info <- function() {
  .amatrix_load_calibration()
  .amatrix_state$calibration
}

#' Report amatrix benchmark status across ops and backends
#'
#' Reads a machine-local benchmark baseline CSV (a table of recorded
#' per-op cold and warm timings), if one is present, together with the
#' cached calibration in the user cache directory, and returns a
#' structured data.frame surfacing per-op cold vs warm timings and the
#' currently-calibrated dispatch thresholds.
#'
#' This is the user-facing honesty surface for Track 4's speed contract:
#' users can see (a) which backends are calibrated on their machine,
#' (b) cold-start vs warm-run ratios per op, and (c) where the
#' dispatcher will currently route.
#'
#' @param baseline_path Path to a benchmark baseline CSV of recorded
#'   per-op timings. Defaults to a \code{baseline.csv} looked up relative
#'   to the current working directory. Pass \code{NULL} to skip baseline
#'   reading entirely and return only calibration data.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{baseline}{data.frame with columns \code{op}, \code{size},
#'       \code{backend}, \code{cold_ms}, \code{warm_ms},
#'       \code{warm_vs_cold_ratio}, \code{speedup_vs_cpu}. Rows with
#'       missing cold OR warm data use \code{NA} for the missing variant.
#'       Empty when the baseline file is absent.}
#'     \item{calibration}{data.frame with columns \code{backend}, \code{op},
#'       \code{threshold_elements}, \code{gpu_wins}. Rows come from the
#'       cached calibration; empty when no calibration is available.}
#'   }
#'
#' @examples
#' \dontrun{
#' rep <- amatrix_benchmark_report()
#' head(rep$baseline)
#' head(rep$calibration)
#' }
#'
#' @seealso \code{\link{amatrix_calibrate}},
#'   \code{\link{amatrix_calibration_info}}
#' @export
amatrix_benchmark_report <- function(baseline_path = file.path("tools", "baseline.csv")) {
  baseline_df <- data.frame(
    op = character(0), size = character(0), backend = character(0),
    cold_ms = numeric(0), warm_ms = numeric(0),
    warm_vs_cold_ratio = numeric(0), speedup_vs_cpu = numeric(0),
    stringsAsFactors = FALSE
  )

  if (!is.null(baseline_path) && file.exists(baseline_path)) {
    raw <- tryCatch(
      utils::read.csv(baseline_path, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (!is.null(raw) && all(c("op", "size", "backend", "variant", "median_ms") %in% names(raw))) {
      # Pivot cold/warm rows into side-by-side columns. Speedup is backend-local.
      key <- paste(raw$op, raw$size, raw$backend, sep = "\001")
      cold_idx <- raw$variant == "cold"
      warm_idx <- raw$variant == "warm"
      cold_map <- stats::setNames(raw$median_ms[cold_idx], key[cold_idx])
      warm_map <- stats::setNames(raw$median_ms[warm_idx], key[warm_idx])
      all_keys <- unique(c(names(cold_map), names(warm_map)))

      # Carry speedup_vs_cpu from whichever variant is present (prefer cold).
      speedup_cold <- stats::setNames(
        raw$speedup_vs_cpu[cold_idx] %||% NA_real_,
        key[cold_idx]
      )
      speedup_warm <- stats::setNames(
        raw$speedup_vs_cpu[warm_idx] %||% NA_real_,
        key[warm_idx]
      )

      parts <- strsplit(all_keys, "\001", fixed = TRUE)
      ops      <- vapply(parts, `[[`, character(1), 1L)
      sizes    <- vapply(parts, `[[`, character(1), 2L)
      backends <- vapply(parts, `[[`, character(1), 3L)

      cold_vec <- cold_map[all_keys]
      warm_vec <- warm_map[all_keys]
      speedup  <- ifelse(is.na(speedup_cold[all_keys]), speedup_warm[all_keys], speedup_cold[all_keys])
      ratio    <- ifelse(
        is.finite(cold_vec) & is.finite(warm_vec) & warm_vec > 0,
        cold_vec / warm_vec,
        NA_real_
      )

      baseline_df <- data.frame(
        op = ops,
        size = sizes,
        backend = backends,
        cold_ms = unname(cold_vec),
        warm_ms = unname(warm_vec),
        warm_vs_cold_ratio = unname(ratio),
        speedup_vs_cpu = unname(speedup),
        stringsAsFactors = FALSE
      )
      rownames(baseline_df) <- NULL
      # Stable ordering for readability.
      baseline_df <- baseline_df[order(baseline_df$backend, baseline_df$op, baseline_df$size), , drop = FALSE]
      rownames(baseline_df) <- NULL
    }
  }

  .amatrix_load_calibration()
  calibration_df <- data.frame(
    backend = character(0), op = character(0),
    threshold_elements = numeric(0), gpu_wins = logical(0),
    stringsAsFactors = FALSE
  )
  cal <- .amatrix_state$calibration
  if (!is.null(cal) && !is.null(cal$thresholds)) {
    rows <- list()
    for (be in names(cal$thresholds)) {
      for (op in names(cal$thresholds[[be]])) {
        thresh <- cal$thresholds[[be]][[op]]
        rows[[length(rows) + 1L]] <- data.frame(
          backend = be,
          op = op,
          threshold_elements = if (is.infinite(thresh)) Inf else as.numeric(thresh),
          gpu_wins = !is.infinite(thresh),
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(rows) > 0L) {
      calibration_df <- do.call(rbind, rows)
      rownames(calibration_df) <- NULL
    }
  }

  list(
    baseline = baseline_df,
    calibration = calibration_df
  )
}

# ── Internals ─────────────────────────────────────────────────────────────────

.amatrix_calibration_path <- function() {
  cache_dir <- tryCatch(
    tools::R_user_dir("amatrix", "cache"),
    error = function(e) file.path(Sys.getenv("HOME"), ".amatrix", "cache")
  )
  file.path(cache_dir, "calibration.rds")
}

# Hardware / runtime identity hash. When this changes we invalidate cached
# calibration data — thresholds measured on one CPU class are not portable
# to another. Track 4.
.amatrix_sys_hash <- function() {
  # Stable identity string for cache invalidation. We don't need a
  # cryptographic hash — just a short, deterministic fingerprint of the
  # fields that would change if the calibration target moves. Using the
  # fields directly (joined by pipe) keeps the implementation base-R-only
  # and makes hash mismatches self-documenting in the cache file.
  paste0(
    "v1|",
    Sys.info()[["sysname"]], "|",
    Sys.info()[["release"]], "|",
    Sys.info()[["machine"]], "|",
    R.version$platform, "|",
    R.version$major, ".", R.version$minor
  )
}

# Lazy-load calibration from disk on first use. Invalidates the cached file
# when the sys_hash changes (hardware or R version).
.amatrix_load_calibration <- function() {
  if (!is.null(.amatrix_state$calibration)) return(invisible(NULL))
  path <- .amatrix_calibration_path()
  if (!file.exists(path)) return(invisible(NULL))
  cal <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(cal)) return(invisible(NULL))

  # Version gate — v1 is grandfathered with NULL hash (will re-calibrate on
  # next explicit amatrix_calibrate() call). v2+ carries sys_hash.
  if (!identical(cal$version, "1") && !identical(cal$version, "2")) {
    return(invisible(NULL))
  }

  # Hardware invalidation for v2+: skip the cache if the hash does not match.
  if (identical(cal$version, "2")) {
    current_hash <- tryCatch(.amatrix_sys_hash(), error = function(e) NULL)
    if (is.null(current_hash) || !identical(cal$sys_hash, current_hash)) {
      # Stale — do NOT populate state. Caller must re-calibrate.
      return(invisible(NULL))
    }
  }

  .amatrix_state$calibration <- cal
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

  workload <- .amatrix_dispatch_workload(x, op, y = y)
  if (identical(workload, 0L)) return(TRUE)

  workload >= thresh
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

# For each op, find the smallest element count from which gpu_wins remains
# TRUE for all larger tested workloads. Returns Inf if GPU never wins.
.amatrix_derive_thresholds <- function(results, ops) {
  thresholds <- list()
  for (op in ops) {
    sub <- results[results$op == op, , drop = FALSE]
    if (nrow(sub) == 0L) next
    sub <- sub[order(sub$elements), ]
    winners <- which(sub$gpu_wins)
    if (length(winners) == 0L) {
      thresholds[[op]] <- Inf
      next
    }

    threshold_idx <- NA_integer_
    for (idx in winners) {
      if (all(sub$gpu_wins[idx:nrow(sub)])) {
        threshold_idx <- idx
        break
      }
    }

    thresholds[[op]] <- if (is.na(threshold_idx)) Inf else sub$elements[[threshold_idx]]
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
