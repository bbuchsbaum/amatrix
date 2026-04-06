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

chol_factor <- function(X) {
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

chol_solve <- function(factor, B) {
  if (!inherits(factor, "amChol")) {
    stop("factor must be an amChol object", call. = FALSE)
  }

  R <- factor@factor
  B_in <- B
  B_mat <- if (is.vector(B)) matrix(B, ncol = 1L) else as.matrix(B)

  # GPU path: dispatch through the backend's chol_solve_factor when the factor
  # was computed in fast mode on a GPU-capable backend.
  x <- if (isTRUE(factor@precision == "fast") &&
           nzchar(factor@backend) && factor@backend != "cpu") {
    backend <- tryCatch(
      amatrix:::.amatrix_get_backend(factor@backend),
      error = function(e) NULL
    )
    if (!is.null(backend) && is.function(backend$chol_solve_factor)) {
      tryCatch(
        backend$chol_solve_factor(R, B_mat),
        error = function(e) {
          z <- forwardsolve(t(R), B_mat)
          backsolve(R, z)
        }
      )
    } else {
      z <- forwardsolve(t(R), B_mat)
      backsolve(R, z)
    }
  } else {
    # CPU path: standard triangular solve
    z <- forwardsolve(t(R), B_mat)
    backsolve(R, z)
  }

  if (is.vector(B_in) && ncol(x) == 1L) x <- drop(x)
  x
}

chol_diag <- function(factor) {
  if (!inherits(factor, "amChol")) {
    stop("factor must be an amChol object", call. = FALSE)
  }
  diag(factor@factor)
}

chol_logdet <- function(factor) {
  if (!inherits(factor, "amChol")) {
    stop("factor must be an amChol object", call. = FALSE)
  }
  2 * sum(log(diag(factor@factor)))
}

quad_form <- function(factor, v) {
  if (!inherits(factor, "amChol")) {
    stop("factor must be an amChol object", call. = FALSE)
  }
  z <- chol_solve(factor, v)
  if (is.vector(v)) {
    as.double(crossprod(v, z))
  } else {
    crossprod(as.matrix(v), as.matrix(z))
  }
}
