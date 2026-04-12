.amatrix_valid_policies   <- c("auto", "cpu", "mlx", "metal", "arrayfire", "torch")
.amatrix_valid_precisions <- c("strict", "fast")
.amatrix_valid_modes      <- c("exact", "balanced", "fast")

# Resolve mode= + backend= into the (preferred_backend, policy, precision) triple
# used by the internal constructors.
#
# mode="exact"    — strict float64, CPU-pinned. No GPU, no silent downcast.
# mode="balanced" — strict float64, auto routing. GPU where numerically safe.
# mode="fast"     — fast (float32-oriented), auto routing. Full GPU throughput.
#
# backend= is an escape hatch that overrides the mode-derived preferred_backend.
# The legacy preferred_backend/policy/precision params still work and take
# precedence over mode-derived values when explicitly supplied.
.amatrix_resolve_mode <- function(mode, backend, preferred_backend, policy, precision) {
  if (!is.null(mode)) {
    mode <- match.arg(mode, .amatrix_valid_modes)
    derived_precision <- switch(mode,
      exact    = "strict",
      balanced = "strict",
      fast     = "fast"
    )
    derived_backend <- switch(mode,
      exact    = "cpu",             # hard-pinned: exact means CPU semantics
      balanced = "cpu",             # M8: auto-route to GPU where float64-safe; currently CPU-pinned
      fast     = "cpu"              # backend= needed to route to GPU
    )
    pb  <- if (!is.null(preferred_backend)) preferred_backend else if (!is.null(backend)) backend else derived_backend
    pol <- if (!is.null(policy)) policy else amatrix_default_policy()
    pre <- if (!is.null(precision)) precision else derived_precision
    return(list(preferred_backend = pb, policy = pol, precision = pre))
  }
  # No mode: honour explicit args or fall back to package defaults.
  list(
    preferred_backend = if (!is.null(preferred_backend)) preferred_backend else if (!is.null(backend)) backend else "cpu",
    policy            = if (!is.null(policy)) policy else amatrix_default_policy(),
    precision         = if (!is.null(precision)) precision else amatrix_default_precision()
  )
}

#' Get the session-level default dispatch policy
#'
#' Returns the dispatch policy used when an \code{aMatrix} object does
#' not specify its own policy. The policy controls which backend is
#' preferred for operations on new matrices.
#'
#' @return Character string, one of \code{"auto"}, \code{"cpu"},
#'   \code{"mlx"}, \code{"metal"}, \code{"arrayfire"}, or
#'   \code{"torch"}.
#'
#' @examples
#' amatrix_default_policy()
#'
#' @seealso \code{\link{amatrix_set_default_policy}},
#'   \code{\link{amatrix_default_precision}}
#' @export
amatrix_default_policy <- function() {
  .amatrix_state$default_policy
}

#' Set the session-level default dispatch policy
#'
#' Sets the dispatch policy applied to new \code{aMatrix} objects that
#' do not specify their own policy. The change affects all subsequent
#' matrix constructions in the current session.
#'
#' @param policy Character string. Must be one of \code{"auto"},
#'   \code{"cpu"}, \code{"mlx"}, \code{"metal"}, \code{"arrayfire"},
#'   or \code{"torch"}.
#'
#' @return Invisibly, \code{policy}.
#'
#' @examples
#' old <- amatrix_default_policy()
#' amatrix_set_default_policy("auto")
#' amatrix_set_default_policy(old) # restore
#'
#' @seealso \code{\link{amatrix_default_policy}},
#'   \code{\link{amatrix_set_default_precision}}
#' @export
amatrix_set_default_policy <- function(policy) {
  stopifnot(is.character(policy), length(policy) == 1L, nzchar(policy))
  if (!(policy %in% .amatrix_valid_policies)) {
    stop(sprintf("policy must be one of: %s", paste(.amatrix_valid_policies, collapse = ", ")))
  }
  .amatrix_state$default_policy <- policy
  invisible(policy)
}

#' Get the session-level default precision mode
#'
#' Returns the precision mode used when constructing new \code{aMatrix}
#' objects that do not specify their own precision.
#'
#' @return Character string, either \code{"strict"} (double precision)
#'   or \code{"fast"} (single/mixed precision).
#'
#' @examples
#' amatrix_default_precision()
#'
#' @seealso \code{\link{amatrix_set_default_precision}},
#'   \code{\link{amatrix_default_policy}}
#' @export
amatrix_default_precision <- function() {
  .amatrix_state$default_precision
}

.amatrix_backend_precision_modes_safe <- function(backend_name) {
  if (is.null(backend_name) ||
      !nzchar(backend_name) ||
      backend_name %in% c("auto", "cpu")) {
    return(NULL)
  }

  backend <- tryCatch(.amatrix_get_backend(backend_name), error = function(e) NULL)
  if (is.null(backend) || !is.function(backend$precision_modes)) {
    return(NULL)
  }

  modes <- tryCatch(unique(backend$precision_modes()), error = function(e) character())
  modes[modes %in% .amatrix_valid_precisions]
}

.amatrix_resolve_backend_precision <- function(backend_name, precision, precision_missing = FALSE) {
  if (!isTRUE(precision_missing)) {
    return(precision)
  }

  modes <- .amatrix_backend_precision_modes_safe(backend_name)
  if (is.null(modes) || length(modes) == 0L || precision %in% modes) {
    return(precision)
  }

  if ("fast" %in% modes) {
    return("fast")
  }

  modes[[1L]]
}

.amatrix_check_backend_precision <- function(backend_name, precision) {
  modes <- .amatrix_backend_precision_modes_safe(backend_name)
  if (is.null(modes) || length(modes) == 0L || precision %in% modes) {
    return(invisible(TRUE))
  }

  stop(sprintf(
    "backend '%s' does not support precision '%s' (supports: %s)",
    backend_name,
    precision,
    paste(modes, collapse = ", ")
  ), call. = FALSE)
}

#' Set the session-level default precision mode
#'
#' Sets the precision mode applied to new \code{aMatrix} objects that
#' do not specify their own precision. Use \code{"strict"} for
#' reproducible double-precision results and \code{"fast"} for maximum
#' GPU throughput with single/mixed precision.
#'
#' @param precision Character string. Must be one of \code{"strict"}
#'   or \code{"fast"}.
#'
#' @return Invisibly, \code{precision}.
#'
#' @examples
#' old <- amatrix_default_precision()
#' amatrix_set_default_precision("strict")
#' amatrix_set_default_precision(old) # restore
#'
#' @seealso \code{\link{amatrix_default_precision}},
#'   \code{\link{amatrix_set_default_policy}}
#' @export
amatrix_set_default_precision <- function(precision) {
  stopifnot(is.character(precision), length(precision) == 1L, nzchar(precision))
  if (!(precision %in% .amatrix_valid_precisions)) {
    stop(sprintf("precision must be one of: %s", paste(.amatrix_valid_precisions, collapse = ", ")))
  }
  .amatrix_state$default_precision <- precision
  invisible(precision)
}

.amatrix_backend_preference <- function(x, op = NULL) {
  pinned_backend <- .amatrix_live_resident_backend(x)
  if (!is.null(pinned_backend)) {
    return(unique(c(pinned_backend, "cpu")))
  }

  unique(c(x@preferred_backend, x@policy, amatrix_default_policy(), "cpu"))
}
