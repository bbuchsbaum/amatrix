# amatrix_memory_stats() -- snapshot of GPU residency and model cache usage.
# amatrix_gc()          -- free dead residency entries and optionally clear cache.

#' Report GPU residency and model cache usage
#'
#' Returns a snapshot of GPU device memory usage and model cache
#' occupancy for the current session. Backends that expose a
#' \code{memory_usage()} method contribute device byte counts;
#' backends without that method show \code{NA} for byte fields.
#'
#' @return An object of class \code{amatrix_memory_stats}, which is a
#'   list with two components:
#'   \describe{
#'     \item{residency}{data.frame with one row per registered backend
#'       and columns \code{backend} (character), \code{resident_objects}
#'       (integer count of GPU-resident R objects), \code{bytes_used}
#'       (numeric, device bytes in use, or \code{NA}), and
#'       \code{bytes_total} (numeric, total device capacity, or
#'       \code{NA}).}
#'     \item{model_cache}{List with \code{n_entries} (integer, number
#'       of cached matrix factors) and \code{max_size} (integer or
#'       \code{Inf}, the cache size limit).}
#'   }
#'
#' @examples
#' stats <- amatrix_memory_stats()
#' print(stats)
#'
#' @seealso \code{\link{amatrix_gc}}, \code{\link{amatrix_residency_info}}
#' @export
amatrix_memory_stats <- function() {
  res_env <- .amatrix_state$residency
  keys    <- ls(envir = res_env, all.names = FALSE)

  entries <- lapply(keys, function(k) {
    get0(k, envir = res_env, inherits = FALSE)
  })

  all_backends <- tryCatch(amatrix_backend_names(), error = function(e) "cpu")
  # Also include any backend that has a live residency entry (e.g. unloaded)
  resident_backends <- unique(vapply(entries, function(e) e$backend, character(1)))
  all_backends <- unique(c(all_backends, resident_backends))

  rows <- lapply(all_backends, function(be) {
    n_res <- sum(vapply(entries, function(e) identical(e$backend, be), logical(1)))

    be_obj <- tryCatch(.amatrix_get_backend(be), error = function(e) NULL)
    bytes_used  <- NA_real_
    bytes_total <- NA_real_

    if (!is.null(be_obj) && is.function(be_obj[["memory_usage"]])) {
      mem <- tryCatch(be_obj$memory_usage(), error = function(e) NULL)
      if (is.list(mem)) {
        bytes_used  <- as.numeric(mem$used)
        bytes_total <- as.numeric(mem$total)
      }
    }

    data.frame(
      backend          = be,
      resident_objects = n_res,
      bytes_used       = bytes_used,
      bytes_total      = bytes_total,
      stringsAsFactors = FALSE
    )
  })

  cache_n <- length(ls(envir = .amatrix_state$model_cache, all.names = FALSE))

  structure(
    list(
      residency   = do.call(rbind, rows),
      model_cache = list(n_entries = cache_n, max_size = .amatrix_cache_max_size())
    ),
    class = "amatrix_memory_stats"
  )
}

#' @export
print.amatrix_memory_stats <- function(x, ...) {
  cat("-- amatrix memory stats ----------------------------------------\n")

  max_s <- x$model_cache$max_size
  max_label <- if (is.infinite(max_s)) "unlimited" else as.character(max_s)
  cat(sprintf("  model cache: %d entries (max: %s)\n",
    x$model_cache$n_entries, max_label))

  cat("  residency:\n")
  for (i in seq_len(nrow(x$residency))) {
    r <- x$residency[i, ]
    bytes_str <- if (is.na(r$bytes_used)) "" else {
      used_mb  <- r$bytes_used  / 1048576
      total_mb <- r$bytes_total / 1048576
      sprintf("  %.1f / %.1f MB", used_mb, total_mb)
    }
    cat(sprintf("    %-12s  %d resident object(s)%s\n",
      r$backend, r$resident_objects, bytes_str))
  }
  cat("----------------------------------------------------------------\n")
  invisible(x)
}

#' Free stale GPU residency entries and optionally flush the model cache
#'
#' Scans the residency registry and removes entries whose backend no
#' longer reports the associated device buffer as present. Such stale
#' entries arise when a backend is unloaded or the device is reset
#' between sessions. Optionally flushes all cached matrix factors (QR,
#' Cholesky, SVD) from the session model cache.
#'
#' @param cache Logical. If \code{TRUE}, also remove all model-cache
#'   entries. Default \code{FALSE}.
#'
#' @return Invisibly, a list with two integer elements:
#'   \describe{
#'     \item{dead_entries}{Number of stale residency slots removed.}
#'     \item{cache_entries_cleared}{Number of model-cache entries
#'       flushed (0 when \code{cache = FALSE}).}
#'   }
#'
#' @examples
#' amatrix_gc()
#' amatrix_gc(cache = TRUE)
#'
#' @seealso \code{\link{amatrix_memory_stats}},
#'   \code{\link{amatrix_residency_info}}
#' @export
amatrix_gc <- function(cache = FALSE) {
  res_env <- .amatrix_state$residency
  keys    <- ls(envir = res_env, all.names = FALSE)
  dead    <- 0L

  for (k in keys) {
    entry <- get0(k, envir = res_env, inherits = FALSE)
    if (is.null(entry)) next

    be_obj <- tryCatch(.amatrix_get_backend(entry$backend), error = function(e) NULL)
    alive  <- !is.null(be_obj) &&
              .amatrix_backend_residency_capable(be_obj) &&
              isTRUE(be_obj$resident_has(entry$resident_key))

    if (!alive) {
      rm(list = k, envir = res_env)
      dead <- dead + 1L
    }
  }

  cleared <- 0L
  if (cache) {
    cache_keys <- ls(envir = .amatrix_state$model_cache, all.names = FALSE)
    cleared <- length(cache_keys)
    if (cleared > 0L) {
      rm(list = cache_keys, envir = .amatrix_state$model_cache)
      # Also clear LRU access-time tracking
      atime_env <- .amatrix_state$cache_atime
      if (!is.null(atime_env)) {
        rm(list = ls(envir = atime_env), envir = atime_env)
      }
    }
  }

  invisible(list(dead_entries = dead, cache_entries_cleared = cleared))
}
