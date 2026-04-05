.validate_amatrix_slots <- function(object) {
  if (!identical(length(object@preferred_backend), 1L) || is.na(object@preferred_backend)) {
    return("preferred_backend must be a single non-missing string")
  }

  if (!identical(length(object@policy), 1L) || is.na(object@policy)) {
    return("policy must be a single non-missing string")
  }

  if (!(object@policy %in% .amatrix_valid_policies)) {
    return(sprintf("policy must be one of: %s", paste(.amatrix_valid_policies, collapse = ", ")))
  }

  if (!identical(length(object@precision), 1L) || is.na(object@precision)) {
    return("precision must be a single non-missing string")
  }

  if (!(object@precision %in% .amatrix_valid_precisions)) {
    return(sprintf("precision must be one of: %s", paste(.amatrix_valid_precisions, collapse = ", ")))
  }

  if (!identical(length(object@object_id), 1L) || is.na(object@object_id) || !nzchar(object@object_id)) {
    return("object_id must be a single non-missing string")
  }

  TRUE
}

setClass(
  "aMatrix",
  contains = "VIRTUAL",
  slots = c(
    preferred_backend = "character",
    policy = "character",
    precision = "character",
    object_id = "character",
    finalizer_env = "environment"
  ),
  prototype = list(
    preferred_backend = "cpu",
    policy = "auto",
    precision = "strict",
    object_id = "",
    finalizer_env = new.env(parent = emptyenv())
  ),
  validity = .validate_amatrix_slots
)

setClass(
  "adgeMatrix",
  contains = c("aMatrix", "dgeMatrix")
)

setClass(
  "adgCMatrix",
  contains = c("aMatrix", "dgCMatrix")
)

setMethod("show", "adgeMatrix", function(object) {
  cat(sprintf(
    "An amatrix dense matrix [%s|policy=%s|precision=%s]\n",
    object@preferred_backend,
    object@policy,
    object@precision
  ))
  callNextMethod()
})

setMethod("show", "adgCMatrix", function(object) {
  cat(sprintf(
    "An amatrix sparse matrix [%s|policy=%s|precision=%s]\n",
    object@preferred_backend,
    object@policy,
    object@precision
  ))
  callNextMethod()
})
