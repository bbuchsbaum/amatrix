.amatrix_valid_policies   <- c("auto", "cpu", "mlx", "metal", "arrayfire", "opencl", "torch")
.amatrix_valid_precisions <- c("strict", "fast")
.amatrix_valid_modes      <- c("exact", "balanced", "fast")

.amatrix_default_preferred_backend <- function(
  policy = amatrix_default_policy(),
  precision = amatrix_default_precision()
) {
  if (is.null(policy) || !nzchar(policy) || identical(policy, "auto")) {
    if (identical(precision, "fast")) {
      return(.amatrix_default_fast_backend())
    }

    return("cpu")
  }

  policy
}

.amatrix_auto_fast_backend_order <- function() {
  c("mlx", "metal", "arrayfire", "opencl", "torch")
}

.amatrix_default_fast_backend <- function() {
  registered <- setdiff(amatrix_backend_names(), "cpu")
  ordered <- unique(c(
    intersect(.amatrix_auto_fast_backend_order(), registered),
    setdiff(registered, .amatrix_auto_fast_backend_order())
  ))

  for (backend_name in ordered) {
    backend <- tryCatch(.amatrix_get_backend(backend_name), error = function(e) NULL)
    if (is.null(backend) || !isTRUE(backend$available())) {
      next
    }

    modes <- tryCatch(unique(backend$precision_modes()), error = function(e) character())
    if ("fast" %in% modes) {
      return(backend_name)
    }
  }

  "cpu"
}

# Resolve mode= + backend= into the (preferred_backend, policy, precision) triple
# used by the internal constructors.
#
# mode="exact"    — strict float64, CPU-pinned. No GPU, no silent downcast.
# mode="balanced" — DEPRECATED. Historical "strict float64, auto-route to GPU
#                   where numerically safe" mode. The float64 GPU routing was
#                   never implemented (routes to CPU just like "exact"), so
#                   the mode is a source of user-expectation mismatch. Track 5
#                   maps `balanced` -> `exact` with a one-time deprecation
#                   warning per session. See planning_docs/quality-tracking.md
#                   §8. The `balanced` string is still accepted by
#                   .amatrix_valid_modes for backward compatibility.
# mode="fast"     — fast (float32-oriented), auto routing. Full GPU throughput.
#
# backend= is an escape hatch that overrides the mode-derived preferred_backend.
# The legacy preferred_backend/policy/precision params still work and take
# precedence over mode-derived values when explicitly supplied.

.amatrix_balanced_deprecation_warned <- function() {
  isTRUE(.amatrix_state$balanced_deprecation_warned)
}

.amatrix_warn_balanced_deprecation_once <- function() {
  if (.amatrix_balanced_deprecation_warned()) return(invisible())
  .amatrix_state$balanced_deprecation_warned <- TRUE
  warning(
    "mode = \"balanced\" is deprecated and will be removed in a future ",
    "release. It currently behaves identically to mode = \"exact\" ",
    "(CPU-pinned, strict float64); the float64 GPU routing it implied ",
    "was never implemented. Use mode = \"exact\" for strict float64 ",
    "CPU semantics or mode = \"fast\" for GPU throughput. See ",
    "planning_docs/quality-tracking.md \u00a78.",
    call. = FALSE
  )
}

.amatrix_resolve_mode <- function(mode, backend, preferred_backend, policy, precision) {
  if (!is.null(mode)) {
    mode <- match.arg(mode, .amatrix_valid_modes)
    if (identical(mode, "balanced")) {
      .amatrix_warn_balanced_deprecation_once()
      mode <- "exact"
    }
    pol <- if (!is.null(policy)) policy else amatrix_default_policy()
    derived_precision <- switch(mode,
      exact    = "strict",
      fast     = "fast"
    )
    derived_backend <- switch(mode,
      exact    = "cpu",             # hard-pinned: exact means CPU semantics
      fast     = .amatrix_default_preferred_backend(pol, "fast")
    )
    pb  <- if (!is.null(preferred_backend)) preferred_backend else if (!is.null(backend)) backend else derived_backend
    pre <- if (!is.null(precision)) precision else derived_precision
    return(list(preferred_backend = pb, policy = pol, precision = pre))
  }
  # No mode: honour explicit args or fall back to package defaults.
  pol <- if (!is.null(policy)) policy else amatrix_default_policy()
  pre <- if (!is.null(precision)) precision else amatrix_default_precision()
  list(
    preferred_backend = if (!is.null(preferred_backend)) preferred_backend else if (!is.null(backend)) backend else .amatrix_default_preferred_backend(pol, pre),
    policy            = pol,
    precision         = pre
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
    stop(errorCondition(
      sprintf("policy must be one of: %s", paste(.amatrix_valid_policies, collapse = ", ")),
      class = "amatrix_bad_arg",
      call = NULL
    ))
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
    stop(errorCondition(
      sprintf("precision must be one of: %s", paste(.amatrix_valid_precisions, collapse = ", ")),
      class = "amatrix_bad_arg",
      call = NULL
    ))
  }
  .amatrix_state$default_precision <- precision
  invisible(precision)
}

#' Evaluate code with temporary amatrix defaults
#'
#' Temporarily overrides the session-default dispatch policy and/or
#' precision mode for the duration of \code{code}, then restores the
#' previous values on exit, even when \code{code} errors.
#'
#' @param policy Optional temporary policy. Must be one of
#'   \code{"auto"}, \code{"cpu"}, \code{"mlx"}, \code{"metal"},
#'   \code{"arrayfire"}, or \code{"torch"}.
#' @param precision Optional temporary precision. Must be either
#'   \code{"strict"} or \code{"fast"}.
#' @param code Expression to evaluate under the temporary defaults.
#'
#' @return The result of evaluating \code{code}.
#'
#' @examples
#' with_amatrix(policy = "auto", precision = "fast", {
#'   adgeMatrix(matrix(1:4, nrow = 2))
#' })
#'
#' @seealso \code{\link{adgeMatrix}},
#'   \code{\link{amatrix_set_default_policy}},
#'   \code{\link{amatrix_set_default_precision}}
#' @export
with_amatrix <- function(policy = NULL, precision = NULL, code) {
  old_policy <- amatrix_default_policy()
  old_precision <- amatrix_default_precision()

  on.exit({
    amatrix_set_default_policy(old_policy)
    amatrix_set_default_precision(old_precision)
  }, add = TRUE)

  if (!is.null(policy)) {
    amatrix_set_default_policy(policy)
  }
  if (!is.null(precision)) {
    amatrix_set_default_precision(precision)
  }

  force(code)
}

.amatrix_backend_preference <- function(x, op = NULL) {
  pinned_backend <- .amatrix_live_resident_backend(x)
  explicit_policy <- nzchar(x@policy) && !identical(x@policy, "auto")
  force_cpu <- explicit_policy && identical(x@policy, "cpu")
  preferred_candidates <- function(...) {
    candidates <- c(...)
    candidates <- candidates[!is.na(candidates) & nzchar(candidates) & candidates != "auto"]
    unique(candidates)
  }

  if (!is.null(pinned_backend)) {
    pinned_available <- tryCatch(
      isTRUE(.amatrix_get_backend(pinned_backend)$available()),
      error = function(e) FALSE
    )

    if (isTRUE(pinned_available)) {
      return(preferred_candidates(
        pinned_backend,
        "cpu"
      ))
    }

    return(preferred_candidates(
      pinned_backend,
      if (force_cpu) "cpu" else x@preferred_backend,
      if (explicit_policy && !force_cpu) x@policy else NULL,
      if (force_cpu) x@preferred_backend else NULL,
      amatrix_default_policy(),
      "cpu"
    ))
  }

  preferred_candidates(
    if (force_cpu) "cpu" else x@preferred_backend,
    if (explicit_policy && !force_cpu) x@policy else NULL,
    if (force_cpu) x@preferred_backend else NULL,
    amatrix_default_policy(),
    "cpu"
  )
}
