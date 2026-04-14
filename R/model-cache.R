# Centralised model-cache helpers with LRU eviction.
#
# The model cache stores amQR, amChol, and amSVD factors keyed by a string
# derived from the source matrix's object_id. It is an R environment
# (.amatrix_state$model_cache) shared by chol-factor.R, svd-factor.R, and
# models-lm.R.
#
# By default the cache is unbounded (max_size = Inf). When a finite max_size
# is configured, the least-recently-used entry is evicted whenever a new entry
# would exceed the limit. Access order is tracked via a monotonically
# increasing integer counter stored in .amatrix_state$cache_atime.

# ── LRU state initialisation (called from zzz.R .onLoad) ──────────────────
.amatrix_cache_init <- function() {
  if (is.null(.amatrix_state$cache_max_size)) {
    .amatrix_state$cache_max_size    <- Inf
  }
  if (is.null(.amatrix_state$cache_atime)) {
    .amatrix_state$cache_atime       <- new.env(parent = emptyenv())
  }
  if (is.null(.amatrix_state$cache_atime_counter)) {
    .amatrix_state$cache_atime_counter <- 0L
  }
}

# ── Public getters/setters ─────────────────────────────────────────────────

#' Get or set the model cache maximum size
#'
#' \code{amatrix_cache_max_size} returns the current limit.
#' \code{amatrix_set_cache_max_size} changes the limit and immediately
#' evicts the least-recently-used entries if the cache exceeds the new
#' bound. When \code{max_size} is \code{Inf} (the default) the cache
#' grows without bound.
#'
#' @param max_size Positive numeric scalar or \code{Inf}; the maximum
#'   number of factorizations to retain in the model cache.
#'
#' @return \code{amatrix_cache_max_size} returns a length-1 numeric
#'   giving the current limit. \code{amatrix_set_cache_max_size}
#'   returns the new limit invisibly.
#'
#' @examples
#' old <- amatrix_cache_max_size()
#' amatrix_set_cache_max_size(10)
#' amatrix_cache_max_size()
#' amatrix_set_cache_max_size(old)
#'
#' @export
amatrix_cache_max_size <- function() {
  .amatrix_cache_max_size()
}

.amatrix_cache_max_size <- function() {
  sz <- .amatrix_state$cache_max_size
  if (is.null(sz)) Inf else sz
}

#' @rdname amatrix_cache_max_size
#' @export
amatrix_set_cache_max_size <- function(max_size) {
  stopifnot(is.numeric(max_size), length(max_size) == 1L, max_size >= 1L)
  .amatrix_state$cache_max_size <- max_size
  # If currently over limit, evict until within budget
  .amatrix_cache_trim()
  invisible(max_size)
}

# ── Internal cache primitives ──────────────────────────────────────────────

.amatrix_cache_get <- function(cache_key) {
  val <- get0(cache_key, envir = .amatrix_state$model_cache, inherits = FALSE)
  if (!is.null(val)) {
    .amatrix_cache_touch(cache_key)
  }
  val
}

.amatrix_cache_set <- function(cache_key, value) {
  max_size <- .amatrix_cache_max_size()

  if (!is.infinite(max_size)) {
    # Check if key already present (update, not new entry)
    already_present <- exists(cache_key,
      envir = .amatrix_state$model_cache, inherits = FALSE)
    if (!already_present) {
      n_current <- length(ls(envir = .amatrix_state$model_cache,
                             all.names = FALSE))
      if (n_current >= max_size) {
        .amatrix_cache_evict_lru()
      }
    }
  }

  assign(cache_key, value, envir = .amatrix_state$model_cache)
  .amatrix_cache_touch(cache_key)
  invisible(value)
}

# Record / refresh the access time for a key using a monotonic counter.
.amatrix_cache_touch <- function(cache_key) {
  atime_env <- .amatrix_state$cache_atime
  if (is.null(atime_env)) return(invisible(NULL))
  ctr <- .amatrix_state$cache_atime_counter + 1L
  .amatrix_state$cache_atime_counter <- ctr
  assign(cache_key, ctr, envir = atime_env)
  invisible(NULL)
}

.amatrix_cache_clear <- function() {
  cache_keys <- ls(envir = .amatrix_state$model_cache, all.names = FALSE)
  if (length(cache_keys) > 0L) {
    rm(list = cache_keys, envir = .amatrix_state$model_cache)
  }

  atime_env <- .amatrix_state$cache_atime
  if (!is.null(atime_env)) {
    atime_keys <- ls(envir = atime_env, all.names = FALSE)
    if (length(atime_keys) > 0L) {
      rm(list = atime_keys, envir = atime_env)
    }
  }

  .amatrix_state$cache_atime_counter <- 0L
  invisible(NULL)
}

# Evict the single least-recently-used cache entry.
.amatrix_cache_evict_lru <- function() {
  atime_env <- .amatrix_state$cache_atime
  cache_env <- .amatrix_state$model_cache
  if (is.null(atime_env)) return(invisible(NULL))

  akeys <- ls(envir = atime_env, all.names = FALSE)
  if (length(akeys) == 0L) return(invisible(NULL))

  atimes <- vapply(akeys, function(k) {
    get0(k, envir = atime_env, inherits = FALSE) %||% 0L
  }, integer(1))

  lru_key <- akeys[[which.min(atimes)]]

  if (exists(lru_key, envir = cache_env, inherits = FALSE)) {
    rm(list = lru_key, envir = cache_env)
  }
  if (exists(lru_key, envir = atime_env, inherits = FALSE)) {
    rm(list = lru_key, envir = atime_env)
  }
  invisible(lru_key)
}

# Remove entries until cache is within the current max_size limit.
.amatrix_cache_trim <- function() {
  max_size <- .amatrix_cache_max_size()
  if (is.infinite(max_size)) return(invisible(NULL))
  while (length(ls(envir = .amatrix_state$model_cache,
                   all.names = FALSE)) > max_size) {
    .amatrix_cache_evict_lru()
  }
  invisible(NULL)
}

# NULL-coalescing operator (private)
`%||%` <- function(x, y) if (is.null(x)) y else x
