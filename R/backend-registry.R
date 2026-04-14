.amatrix_optional_backends_enabled <- function() {
  !identical(getOption("amatrix.optional_backends", TRUE), FALSE)
}

.amatrix_optional_backend_specs <- function() {
  list(
    mlx = list(
      package = "amatrix.mlx",
      register_fun = "amatrix_mlx_register",
      enabled = function() TRUE
    ),
    metal = list(
      package = "amatrix.metal",
      register_fun = "amatrix_metal_register",
      enabled = function() isTRUE(getOption("amatrix.enable_metal", FALSE))
    ),
    opencl = list(
      package = "amatrix.opencl",
      register_fun = "amatrix_opencl_register",
      enabled = function() isTRUE(getOption("amatrix.enable_opencl", FALSE))
    ),
    arrayfire = list(
      package = "amatrix.arrayfire",
      register_fun = "amatrix_arrayfire_register",
      enabled = function() isTRUE(getOption("amatrix.enable_arrayfire", FALSE))
    )
  )
}

.amatrix_optional_backend_enabled <- function(spec) {
  enabled <- spec$enabled
  if (is.null(enabled)) {
    return(TRUE)
  }
  isTRUE(tryCatch(enabled(), error = function(e) FALSE))
}

.amatrix_try_register_optional_backend <- function(name) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))

  if (!.amatrix_optional_backends_enabled()) {
    return(FALSE)
  }

  if (exists(name, envir = .amatrix_state$backends, inherits = FALSE)) {
    return(TRUE)
  }

  spec <- .amatrix_optional_backend_specs()[[name]]
  if (is.null(spec)) {
    return(FALSE)
  }
  if (!.amatrix_optional_backend_enabled(spec)) {
    return(FALSE)
  }

  lib_locs <- unique(c(.libPaths(), .Library, .Library.site))
  ns <- tryCatch(
    loadNamespace(spec$package, lib.loc = lib_locs),
    error = function(e) NULL
  )
  if (is.null(ns)) {
    return(FALSE)
  }

  register_backend <- get0(spec$register_fun, envir = ns, inherits = FALSE)
  if (!is.function(register_backend)) {
    return(FALSE)
  }

  isTRUE(tryCatch({
    register_backend(overwrite = TRUE)
    exists(name, envir = .amatrix_state$backends, inherits = FALSE)
  }, error = function(e) FALSE))
}

#' Register a backend with the amatrix dispatch system
#'
#' Adds a named backend to the session backend registry. The backend
#' must be a named list containing all required callable fields. Once
#' registered, the backend is available for dispatch by any
#' \code{aMatrix} object whose \code{preferred_backend} or \code{policy}
#' slot names it.
#'
#' @param name Character string. Unique identifier for the backend
#'   (e.g. \code{"mlx"}, \code{"opencl"}).
#' @param backend Named list implementing the backend contract. Required
#'   fields: \code{capabilities}, \code{features}, \code{precision_modes}
#'   (each a zero-argument function returning a character vector),
#'   \code{available} (zero-argument logical function), \code{supports},
#'   \code{matmul}, \code{crossprod}, \code{tcrossprod}, \code{ewise},
#'   \code{rowSums}, \code{colSums}.
#' @param overwrite Logical. Allow replacement of an existing registration
#'   with the same \code{name}. Default \code{FALSE}.
#'
#' @return Invisibly, \code{name}.
#'
#' @examples
#' # Minimal no-op backend for illustration only
#' noop <- list(
#'   capabilities   = function() character(),
#'   features       = function() character(),
#'   precision_modes = function() "strict",
#'   available      = function() FALSE,
#'   supports       = function(op, x, y = NULL) FALSE,
#'   matmul         = function(x, y) x,
#'   crossprod      = function(x, y = NULL) x,
#'   tcrossprod     = function(x, y = NULL) x,
#'   ewise          = function(x, y, op) x,
#'   rowSums        = function(x) numeric(nrow(x)),
#'   colSums        = function(x) numeric(ncol(x))
#' )
#' amatrix_register_backend("noop_test", noop, overwrite = TRUE)
#'
#' @seealso \code{\link{amatrix_backend_names}},
#'   \code{\link{amatrix_backend_status}}
#' @export
amatrix_register_backend <- function(name, backend, overwrite = FALSE) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))

  .bad_backend <- function(msg) {
    stop(errorCondition(msg, class = "amatrix_bad_backend", call = sys.call(-1L)))
  }

  if (!is.list(backend)) {
    .bad_backend("backend must be a named list")
  }

  required_fields <- c(
    "capabilities",
    "features",
    "precision_modes",
    "available",
    "supports",
    "matmul",
    "crossprod",
    "tcrossprod",
    "ewise",
    "rowSums",
    "colSums"
  )
  missing_fields <- setdiff(required_fields, names(backend))
  if (length(missing_fields) > 0L) {
    .bad_backend(sprintf("backend is missing required fields: %s", paste(missing_fields, collapse = ", ")))
  }

  if (!is.function(backend$capabilities)) {
    stop("backend$capabilities must be a function")
  }
  if (!is.function(backend$features)) {
    stop("backend$features must be a function")
  }
  if (!is.function(backend$precision_modes)) {
    stop("backend$precision_modes must be a function")
  }

  backend_capabilities <- backend$capabilities()
  backend_features <- backend$features()
  backend_precision_modes <- backend$precision_modes()
  if (!is.character(backend_capabilities)) {
    stop("backend$capabilities() must return a character vector")
  }
  if (!is.character(backend_features)) {
    stop("backend$features() must return a character vector")
  }
  if (!is.character(backend_precision_modes)) {
    stop("backend$precision_modes() must return a character vector")
  }
  if (length(backend_precision_modes) == 0L) {
    stop("backend$precision_modes() must return at least one precision mode")
  }
  if (!all(backend_precision_modes %in% .amatrix_valid_precisions)) {
    stop(sprintf(
      "backend$precision_modes() must be a subset of: %s",
      paste(.amatrix_valid_precisions, collapse = ", ")
    ))
  }

  exists_already <- exists(name, envir = .amatrix_state$backends, inherits = FALSE)
  if (exists_already && !overwrite) {
    stop(errorCondition(
      sprintf("backend '%s' is already registered", name),
      class = c("amatrix_backend_exists", "amatrix_bad_backend"),
      call = NULL
    ))
  }

  assign(name, backend, envir = .amatrix_state$backends)
  if (isTRUE(exists_already)) {
    calibration <- .amatrix_state$calibration
    if (!is.null(calibration)) {
      if (!is.null(calibration$thresholds)) {
        calibration$thresholds[[name]] <- NULL
      }
      if (is.data.frame(calibration$results) && "backend" %in% names(calibration$results)) {
        calibration$results <- calibration$results[calibration$results$backend != name, , drop = FALSE]
      }
      .amatrix_state$calibration <- calibration
    }

    if (!is.null(.amatrix_state$backend_health)) {
      .amatrix_state$backend_health[[name]] <- NULL
    }

    .amatrix_cache_clear()
  }
  invisible(name)
}

.amatrix_get_backend <- function(name) {
  backend <- get0(name, envir = .amatrix_state$backends, inherits = FALSE)
  if (is.null(backend) && isTRUE(.amatrix_try_register_optional_backend(name))) {
    backend <- get0(name, envir = .amatrix_state$backends, inherits = FALSE)
  }
  if (is.null(backend)) {
    stop(sprintf("backend '%s' is not registered", name))
  }
  backend
}

#' List names of all registered backends
#'
#' Returns the names of every backend currently in the session registry.
#' When optional backends are enabled (the default), this also attempts
#' to auto-register any installed optional backend packages before
#' returning the list.
#'
#' @return Character vector of registered backend names, sorted
#'   alphabetically. Always includes at least \code{"cpu"}.
#'
#' @examples
#' amatrix_backend_names()
#'
#' @seealso \code{\link{amatrix_backend_status}},
#'   \code{\link{amatrix_register_backend}}
#' @export
amatrix_backend_names <- function() {
  if (.amatrix_optional_backends_enabled()) {
    invisible(lapply(names(.amatrix_optional_backend_specs()), .amatrix_try_register_optional_backend))
  }
  sort(ls(envir = .amatrix_state$backends, all.names = FALSE))
}

#' Query the capabilities of a registered backend
#'
#' Returns the unique capability strings advertised by the named backend,
#' as reported by its \code{capabilities()} function.
#'
#' @param name Character string. Name of a registered backend.
#'
#' @return Character vector of capability identifiers (e.g.
#'   \code{"matmul"}, \code{"svd"}).
#'
#' @examples
#' amatrix_backend_capabilities("cpu")
#'
#' @seealso \code{\link{amatrix_backend_features}},
#'   \code{\link{amatrix_backend_status}}
#' @export
amatrix_backend_capabilities <- function(name) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  backend <- .amatrix_get_backend(name)
  unique(backend$capabilities())
}

#' Query the features of a registered backend
#'
#' Returns the unique feature strings advertised by the named backend,
#' as reported by its \code{features()} function. Features describe
#' optional capabilities such as sparse residency or deferred execution.
#'
#' @param name Character string. Name of a registered backend.
#'
#' @return Character vector of feature identifiers.
#'
#' @examples
#' amatrix_backend_features("cpu")
#'
#' @seealso \code{\link{amatrix_backend_capabilities}},
#'   \code{\link{amatrix_backend_status}}
#' @export
amatrix_backend_features <- function(name) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  backend <- .amatrix_get_backend(name)
  unique(backend$features())
}

#' Query the precision modes supported by a registered backend
#'
#' Returns the precision mode strings advertised by the named backend.
#' Valid values are \code{"strict"} (double precision) and \code{"fast"}
#' (single/mixed precision).
#'
#' @param name Character string. Name of a registered backend.
#'
#' @return Character vector of precision mode identifiers, a subset of
#'   \code{c("strict", "fast")}.
#'
#' @examples
#' amatrix_backend_precision_modes("cpu")
#'
#' @seealso \code{\link{amatrix_backend_capabilities}},
#'   \code{\link{amatrix_backend_status}}
#' @export
amatrix_backend_precision_modes <- function(name) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  backend <- .amatrix_get_backend(name)
  unique(backend$precision_modes())
}

#' Summarise the status of registered backends
#'
#' Returns a data.frame with one row per backend describing its
#' availability, supported precision modes, features, capabilities,
#' and whether it supports GPU residency.
#'
#' @param names Character vector of backend names to query. When
#'   \code{NULL} (default) all registered backends are included,
#'   with optional backends auto-registered first if possible.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{name}{Character. Backend identifier.}
#'     \item{available}{Logical. Whether the backend reports itself
#'       as available on this machine.}
#'     \item{precision_modes}{Character. Comma-separated precision
#'       modes (\code{"strict"}, \code{"fast"}).}
#'     \item{features}{Character. Comma-separated feature strings.}
#'     \item{residency_capable}{Logical. Whether the backend
#'       supports GPU-resident matrix storage.}
#'     \item{capabilities}{Character. Comma-separated operation
#'       capability strings.}
#'   }
#'
#' @examples
#' amatrix_backend_status()
#' amatrix_backend_status("cpu")
#'
#' @seealso \code{\link{amatrix_backend_names}},
#'   \code{\link{amatrix_register_backend}}
#' @export
amatrix_backend_status <- function(names = NULL) {
  if (is.null(names)) {
    if (.amatrix_optional_backends_enabled()) {
      invisible(lapply(names(.amatrix_optional_backend_specs()), .amatrix_try_register_optional_backend))
    }
    names <- amatrix_backend_names()
  }

  stopifnot(is.character(names))

  rows <- lapply(names, function(name) {
    backend <- tryCatch(.amatrix_get_backend(name), error = function(e) NULL)
    if (is.null(backend)) {
      stop(sprintf("backend '%s' is not registered", name))
    }
    health <- .amatrix_backend_health_get(name)
    data.frame(
      name = name,
      available = isTRUE(backend$available()),
      health = health$status,
      health_reason = health$reason %||% NA_character_,
      precision_modes = paste(amatrix_backend_precision_modes(name), collapse = ","),
      features = paste(amatrix_backend_features(name), collapse = ","),
      residency_capable = .amatrix_backend_residency_capable(backend),
      capabilities = paste(amatrix_backend_capabilities(name), collapse = ","),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

# ============================================================================
# Track 5: Backend health, fallback telemetry, and honest tiering
# ============================================================================
#
# Infrastructure for the three Track 5 contracts:
#
#   1. Health state:        Every backend has a health record. Default is
#                           "unprobed". A runtime failure in dispatch marks
#                           it "unhealthy:<reason>"; a successful canary
#                           probe or a clean dispatch marks it "healthy".
#                           amatrix_backend_status() surfaces this.
#
#   2. Fallback telemetry:  When a backend call fails and dispatch falls
#                           back to CPU, a structured event is appended to
#                           .amatrix_state$fallback_log. amatrix_fallback_log()
#                           returns the log as a data.frame. The conformance
#                           suite asserts the log is empty after a clean run —
#                           a non-empty log means a backend claimed support
#                           for an op it cannot actually execute.
#
#   3. Canary probe:        amatrix_backend_health_probe(name) runs a small
#                           10x10 matmul + residual check against base R
#                           and marks the backend healthy/unhealthy.
#
# See planning_docs/quality-tracking.md §7 (stop-ship rule 7) and §8.

.amatrix_backend_health_init <- function() {
  if (is.null(.amatrix_state$backend_health)) {
    .amatrix_state$backend_health <- new.env(parent = emptyenv())
  }
  if (is.null(.amatrix_state$fallback_log)) {
    .amatrix_state$fallback_log <- list()
  }
  invisible()
}

.amatrix_backend_health_mark <- function(name, status, reason = NULL) {
  .amatrix_backend_health_init()
  .amatrix_state$backend_health[[name]] <- list(
    status = status,
    reason = reason,
    timestamp = Sys.time()
  )
  invisible()
}

.amatrix_backend_health_get <- function(name) {
  .amatrix_backend_health_init()
  rec <- .amatrix_state$backend_health[[name]]
  if (is.null(rec)) {
    return(list(status = "unprobed", reason = NA_character_, timestamp = NA))
  }
  rec
}

#' Run a canary health probe against a registered backend
#'
#' Executes a small matmul round-trip against the named backend and
#' compares the result to the base R reference. On success the backend
#' is marked \code{healthy}; on failure it is marked
#' \code{unhealthy:<reason>}. Subsequent calls to
#' \code{amatrix_backend_status()} reflect the recorded health.
#'
#' The probe is intentionally tiny (10x10 double-precision matmul) so it
#' completes in milliseconds even on cold GPU. It is not a benchmark; it
#' is a liveness check.
#'
#' @param name Character string. Name of a registered backend.
#' @param tol Numeric. Residual tolerance for the probe, default
#'   \code{1e-8} (float64) or \code{1e-4} (if the backend only supports
#'   fast precision).
#'
#' @return Invisibly, the health record as a list with elements
#'   \code{status}, \code{reason}, \code{timestamp}.
#'
#' @examples
#' amatrix_backend_health_probe("cpu")
#'
#' @seealso \code{\link{amatrix_backend_status}},
#'   \code{\link{amatrix_fallback_log}}
#' @export
amatrix_backend_health_probe <- function(name, tol = NULL) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))

  backend <- tryCatch(.amatrix_get_backend(name), error = function(e) NULL)
  if (is.null(backend)) {
    .amatrix_backend_health_mark(name, "unhealthy", "not registered")
    return(invisible(.amatrix_backend_health_get(name)))
  }

  if (!isTRUE(backend$available())) {
    .amatrix_backend_health_mark(name, "unhealthy", "backend$available() returned FALSE")
    return(invisible(.amatrix_backend_health_get(name)))
  }

  if (!is.function(backend$matmul)) {
    .amatrix_backend_health_mark(name, "unhealthy", "backend$matmul is not a function")
    return(invisible(.amatrix_backend_health_get(name)))
  }

  # Determine tolerance from precision modes.
  modes <- tryCatch(unique(backend$precision_modes()), error = function(e) character())
  if (is.null(tol)) {
    tol <- if ("strict" %in% modes) 1e-8 else 1e-4
  }

  canary <- tryCatch({
    set.seed(2026041500L)
    x_host <- matrix(rnorm(100L), nrow = 10L, ncol = 10L)
    y_host <- matrix(rnorm(100L), nrow = 10L, ncol = 10L)
    x_a <- new_adgeMatrix(x_host, preferred_backend = name,
                          policy = "cpu", precision = if ("strict" %in% modes) "strict" else "fast")
    y_a <- new_adgeMatrix(y_host, preferred_backend = name,
                          policy = "cpu", precision = if ("strict" %in% modes) "strict" else "fast")
    result <- backend$matmul(x_a, y_a)
    host_result <- if (inherits(result, "aMatrix")) as.matrix(result)
      else if (inherits(result, "Matrix")) as.matrix(result)
      else result
    reference <- x_host %*% y_host
    max_err <- max(abs(host_result - reference))
    if (!is.finite(max_err) || max_err > tol) {
      list(ok = FALSE, reason = sprintf("canary residual %.2e exceeds tol %.2e", max_err, tol))
    } else {
      list(ok = TRUE, reason = NULL)
    }
  }, error = function(e) {
    list(ok = FALSE, reason = sprintf("canary error: %s", conditionMessage(e)))
  })

  if (isTRUE(canary$ok)) {
    .amatrix_backend_health_mark(name, "healthy", NULL)
  } else {
    .amatrix_backend_health_mark(name, "unhealthy", canary$reason)
  }

  invisible(.amatrix_backend_health_get(name))
}

# ---------------------------------------------------------------------------
# Fallback telemetry
# ---------------------------------------------------------------------------

.amatrix_log_fallback <- function(op, backend, reason,
                                  from_backend = NULL, to_backend = "cpu") {
  .amatrix_backend_health_init()
  idx <- length(.amatrix_state$fallback_log) + 1L
  .amatrix_state$fallback_log[[idx]] <- list(
    timestamp = Sys.time(),
    op = op,
    from_backend = from_backend %||% backend,
    to_backend = to_backend,
    reason = reason
  )
  invisible()
}

#' Return the amatrix backend fallback log
#'
#' The fallback log records every runtime fall-through from a preferred
#' backend to the CPU reference path. A non-empty log after a clean
#' conformance run is a stop-ship condition (planning_docs/quality-tracking.md
#' §7 rule 7): it means a backend claimed support for an op it cannot
#' actually execute.
#'
#' @return A data.frame with columns \code{timestamp}, \code{op},
#'   \code{from_backend}, \code{to_backend}, \code{reason}. Zero rows
#'   means no fallbacks have been recorded.
#'
#' @examples
#' amatrix_fallback_log()
#' amatrix_fallback_log_reset()
#'
#' @seealso \code{\link{amatrix_fallback_log_reset}},
#'   \code{\link{amatrix_backend_health_probe}}
#' @export
amatrix_fallback_log <- function() {
  .amatrix_backend_health_init()
  events <- .amatrix_state$fallback_log
  if (length(events) == 0L) {
    return(data.frame(
      timestamp = as.POSIXct(character(0)),
      op = character(0),
      from_backend = character(0),
      to_backend = character(0),
      reason = character(0),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, lapply(events, function(e) {
    data.frame(
      timestamp = e$timestamp,
      op = e$op %||% NA_character_,
      from_backend = e$from_backend %||% NA_character_,
      to_backend = e$to_backend %||% NA_character_,
      reason = e$reason %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }))
}

#' Clear the amatrix backend fallback log
#'
#' Resets the fallback log to empty. Typically called at the start of a
#' test block to isolate the assertion that the log is empty after a
#' clean run.
#'
#' @return Invisibly, \code{NULL}.
#'
#' @seealso \code{\link{amatrix_fallback_log}}
#' @export
amatrix_fallback_log_reset <- function() {
  .amatrix_backend_health_init()
  .amatrix_state$fallback_log <- list()
  invisible()
}
