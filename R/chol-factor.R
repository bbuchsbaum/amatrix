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
  cached <- .amatrix_cache_get(cache_key)
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

  .amatrix_cache_set(cache_key, factor_obj)
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
      .amatrix_get_backend(factor@backend),
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

solve_triangular <- function(R, B, lower = FALSE) {
  R_mat <- if (inherits(R, "amChol")) R@factor else as.matrix(R)
  scalar_out <- is.vector(B) || (is.matrix(B) && ncol(B) == 1L)
  B_mat <- if (is.vector(B)) matrix(B, ncol = 1L) else as.matrix(B)
  x <- if (isTRUE(lower)) {
    forwardsolve(R_mat, B_mat)
  } else {
    backsolve(R_mat, B_mat)
  }
  if (scalar_out && ncol(x) == 1L) drop(x) else x
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

# ---------------------------------------------------------------------------
# LU factorization  (general square systems; mirrors the amChol pattern)
# ---------------------------------------------------------------------------

setClass(
  "amLU",
  slots = list(
    A          = "matrix",    # original square matrix; LAPACK DGESV factorises on solve
    source_id  = "character",
    precision  = "character",
    backend    = "character"
  ),
  prototype = list(
    A          = matrix(numeric(0), 0L, 0L),
    source_id  = NA_character_,
    precision  = NA_character_,
    backend    = NA_character_
  )
)

setMethod("show", "amLU", function(object) {
  cat(sprintf(
    "amLU [%dx%d | %s | source: %s]\n",
    nrow(object@A), ncol(object@A),
    object@precision, object@source_id
  ))
  invisible(object)
})

lu_factor <- function(A) {
  A_mat <- if (inherits(A, "adgeMatrix")) {
    m <- as.matrix(amatrix_materialize_host(A))
    storage.mode(m) <- "double"
    m
  } else {
    m <- as.matrix(A)
    storage.mode(m) <- "double"
    m
  }
  if (nrow(A_mat) != ncol(A_mat)) {
    stop("A must be a square matrix", call. = FALSE)
  }
  src  <- if (inherits(A, "adgeMatrix")) A@object_id else NA_character_
  prec <- if (inherits(A, "adgeMatrix")) A@precision  else NA_character_
  be   <- if (inherits(A, "adgeMatrix")) A@preferred_backend else NA_character_
  new("amLU", A = A_mat, source_id = src, precision = prec, backend = be)
}

lu_solve <- function(factor, B) {
  if (!inherits(factor, "amLU")) {
    stop("factor must be an amLU object", call. = FALSE)
  }
  scalar_out <- is.vector(B) || (is.matrix(B) && ncol(B) == 1L)
  B_mat <- if (is.vector(B)) matrix(B, ncol = 1L) else as.matrix(B)
  x <- base::solve(factor@A, B_mat)
  if (scalar_out && ncol(x) == 1L) drop(x) else x
}
