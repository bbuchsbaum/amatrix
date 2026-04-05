.amatrix_valid_policies <- c("auto", "cpu", "mlx", "arrayfire", "torch")
.amatrix_valid_precisions <- c("strict", "fast")

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
