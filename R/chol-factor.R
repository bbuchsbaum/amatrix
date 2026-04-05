.validate_amChol <- function(object) {
  R <- object@factor
  if (!is.matrix(R)) {
    return("factor must be a matrix")
  }
  if (nrow(R) > 0L && ncol(R) > 0L) {
    if (nrow(R) != ncol(R)) {
      return("factor must be square")
    }
    if (any(abs(R[lower.tri(R)]) > 0)) {
      return("factor must be upper triangular")
    }
  }
  TRUE
}

setClass(
  "amChol",
  slots = list(
    factor = "matrix",
    source_id = "character",
    precision = "character",
    backend = "character"
  ),
  prototype = list(
    factor = matrix(numeric(0), 0L, 0L),
    source_id = NA_character_,
    precision = NA_character_,
    backend = NA_character_
  ),
  validity = .validate_amChol
)

setMethod("show", "amChol", function(object) {
  cat(sprintf(
    "amChol [%dx%d | %s | source: %s]\n",
    nrow(object@factor),
    ncol(object@factor),
    object@precision,
    object@source_id
  ))
  invisible(object)
})

as.matrix.amChol <- function(x, ...) x@factor

am_chol_factor <- function(X) {
  if (!inherits(X, "adgeMatrix")) {
    stop("X must be an adgeMatrix (symmetric positive definite)", call. = FALSE)
  }

  cache_key <- paste0("chol:", X@object_id)
  cached <- get0(cache_key, envir = .amatrix_state$model_cache, inherits = FALSE)
  if (!is.null(cached)) {
    return(cached)
  }

  R <- as.matrix(am_chol(X))

  factor_obj <- new(
    "amChol",
    factor = R,
    source_id = X@object_id,
    precision = X@precision,
    backend = X@preferred_backend
  )

  assign(cache_key, factor_obj, envir = .amatrix_state$model_cache)
  factor_obj
}

am_chol_solve <- function(factor, B) {
  if (!inherits(factor, "amChol")) {
    stop("factor must be an amChol object", call. = FALSE)
  }

  R <- factor@factor
  B_in <- B
  if (is.vector(B)) {
    B <- as.matrix(B)
  } else {
    B <- as.matrix(B)
  }

  # A = R' R, so A x = B solved via forwardsolve(R', B) then backsolve(R, .)
  z <- forwardsolve(t(R), B)
  x <- backsolve(R, z)

  if (is.vector(B_in) && ncol(x) == 1L) {
    x <- as.matrix(x)
  }
  x
}

chol_solve <- function(factor, B) {
  am_chol_solve(factor, B)
}
