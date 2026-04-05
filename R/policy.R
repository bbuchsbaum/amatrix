.amatrix_valid_policies   <- c("auto", "cpu", "mlx", "arrayfire", "torch")
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
      balanced = "cpu",             # auto routing not implemented yet; defaults to cpu
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

amatrix_default_policy <- function() {
  .amatrix_state$default_policy
}

amatrix_set_default_policy <- function(policy) {
  stopifnot(is.character(policy), length(policy) == 1L, nzchar(policy))
  if (!(policy %in% .amatrix_valid_policies)) {
    stop(sprintf("policy must be one of: %s", paste(.amatrix_valid_policies, collapse = ", ")))
  }
  .amatrix_state$default_policy <- policy
  invisible(policy)
}

amatrix_default_precision <- function() {
  .amatrix_state$default_precision
}

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
