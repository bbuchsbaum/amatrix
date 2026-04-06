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
#'   \code{"matmul"}, \code{"crossprod"}, \code{"qr"}, \code{"chol"}.
#' @param sizes List of integer(2) vectors giving (nrow, ncol) test sizes.
#' @param n_reps Number of timed repetitions per cell (after 2 warm-up reps).
#' @param margin Fraction by which GPU must beat CPU to count as a win (default
#'   0.10 = GPU must be at least 10\% faster).
#' @param persist Logical. Save calibration to the user cache directory so
#'   future sessions load it automatically.
#' @param quiet Logical. Suppress progress messages.
#' @return Invisibly, the calibration list (thresholds + full results table).
amatrix_calibrate <- function(
  backend  = NULL,
  ops      = c("matmul", "crossprod", "qr", "chol"),
  sizes    = list(c(64L, 32L), c(128L, 64L), c(256L, 128L), c(512L, 256L), c(1024L, 512L)),
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
      if (!op %in% caps) next

      for (sz in sizes) {
        nr <- as.integer(sz[[1L]])
        nc <- as.integer(sz[[2L]])

        if (!quiet) message(sprintf(
          "  calibrating %s / %s / %dx%d ...", be, op, nr, nc
        ), appendLF = FALSE)

        row <- tryCatch(
          .amatrix_benchmark_op(be, be_obj, op, nr, nc, precision, n_reps),
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

    if (length(be_rows) == 0L) next

    results_be <- do.call(rbind, lapply(be_rows, as.data.frame, stringsAsFactors = FALSE))
    all_rows   <- c(all_rows, be_rows)
    thresholds[[be]] <- .amatrix_derive_thresholds(results_be, ops)
  }

  results_df <- if (length(all_rows) > 0L) {
    do.call(rbind, lapply(all_rows, as.data.frame, stringsAsFactors = FALSE))
  } else {
    data.frame(
      backend = character(), op = character(),
      nrow = integer(), ncol = integer(), elements = integer(),
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

# Returns TRUE when the matrix size clears the calibrated threshold for
# (op, backend). Always returns TRUE when no calibration data exists (backward
# compatible) or for the CPU backend.
.amatrix_calibration_ok <- function(x, op, backend_name) {
  if (backend_name == "cpu") return(TRUE)
  .amatrix_load_calibration()
  cal <- .amatrix_state$calibration
  if (is.null(cal)) return(TRUE)

  thresh <- cal$thresholds[[backend_name]][[op]]
  if (is.null(thresh)) return(TRUE)     # op not benchmarked
  if (is.infinite(thresh)) return(FALSE) # GPU never wins

  nrow(x) * ncol(x) >= thresh
}

# Benchmark one (backend, op, size) cell. Returns a one-row list.
.amatrix_benchmark_op <- function(be_name, be_obj, op, nr, nc, precision, n_reps) {
  X_host  <- matrix(seq_len(nr * nc) / (nr * nc + 1L), nr, nc)
  spd_host <- crossprod(X_host) + diag(nc)  # nc×nc SPD for chol

  X_adge  <- as_adgeMatrix(X_host,  preferred_backend = be_name, precision = precision)
  SPD_adge <- as_adgeMatrix(spd_host, preferred_backend = be_name, precision = precision)

  cpu_fn <- switch(op,
    matmul    = function() base::crossprod(X_host),  # X^T X same FLOP count as X X^T
    crossprod = function() base::crossprod(X_host),
    qr        = function() base::qr(X_host),
    chol      = function() base::chol(spd_host)
  )
  gpu_fn <- switch(op,
    matmul    = function() be_obj$crossprod(X_adge),
    crossprod = function() be_obj$crossprod(X_adge),
    qr        = function() be_obj$qr(X_adge),
    chol      = function() be_obj$chol(SPD_adge)
  )

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
    op       = op,
    nrow     = nr,
    ncol     = nc,
    elements = nr * nc,
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
