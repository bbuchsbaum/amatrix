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
    factor_obj = "ANY",
    source_id = "character",
    precision = "character",
    backend = "character"
  ),
  prototype = list(
    factor = matrix(numeric(0), 0L, 0L),
    factor_obj = NULL,
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

  factor_value <- am_chol(X)
  backend_name <- amatrix_backend_plan(X, "chol")$chosen
  R <- as.matrix(factor_value)

  factor_obj <- new(
    "amChol",
    factor = R,
    factor_obj = if (inherits(factor_value, "aMatrix")) factor_value else NULL,
    source_id = X@object_id,
    precision = X@precision,
    backend = backend_name
  )

  .amatrix_cache_set(cache_key, factor_obj)
  factor_obj
}

.amatrix_resident_triangular_solve <- function(R_obj, B, backend_name, lower = FALSE, transpose = FALSE) {
  backend <- tryCatch(
    .amatrix_get_backend(backend_name),
    error = function(e) NULL
  )
  if (is.null(backend) || !is.function(backend$solve_triangular_resident)) {
    return(NULL)
  }

  factor_arg <- .amatrix_prepare_resident_arg(R_obj, backend_name, promote_amatrix = FALSE)
  rhs_arg <- .amatrix_prepare_resident_arg(B, backend_name, promote_amatrix = FALSE)
  if (is.null(factor_arg) || is.null(rhs_arg)) {
    .amatrix_cleanup_temp_resident(list(rhs_arg), backend_name)
    return(NULL)
  }

  out_key <- .amatrix_next_resident_key(backend_name)
  result <- tryCatch(
    backend$solve_triangular_resident(
      factor_arg$key,
      rhs_arg$key,
      out_key,
      lower = lower,
      transpose = transpose,
      defer = FALSE
    ),
    error = function(e) {
      try(backend$resident_drop(out_key), silent = TRUE)
      NULL
    }
  )
  if (isTRUE(backend$resident_has(out_key))) {
    try(backend$resident_drop(out_key), silent = TRUE)
  }
  .amatrix_cleanup_temp_resident(list(rhs_arg), backend_name)
  result
}

.amatrix_amchol_backend <- function(factor) {
  if (!inherits(factor, "amChol") || !nzchar(factor@backend) || identical(factor@backend, "cpu")) {
    return(NULL)
  }
  tryCatch(.amatrix_get_backend(factor@backend), error = function(e) NULL)
}

.amatrix_amchol_resident_triangular_solve <- function(factor, B_mat, lower = FALSE, transpose = FALSE) {
  factor_obj <- factor@factor_obj
  if (!inherits(factor_obj, "aMatrix")) {
    return(NULL)
  }

  backend_name <- .amatrix_live_resident_backend(factor_obj)
  if (is.null(backend_name)) {
    return(NULL)
  }

  .amatrix_resident_triangular_solve(
    factor_obj,
    B_mat,
    backend_name,
    lower = lower,
    transpose = transpose
  )
}

.amatrix_amchol_resident_solve <- function(factor, B_mat) {
  z <- .amatrix_amchol_resident_triangular_solve(factor, B_mat, lower = FALSE, transpose = TRUE)
  if (is.null(z)) {
    return(NULL)
  }

  x <- .amatrix_amchol_resident_triangular_solve(factor, z, lower = FALSE, transpose = FALSE)
  if (is.null(x)) {
    return(NULL)
  }

  x
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
    backend <- .amatrix_amchol_backend(factor)
    if (!is.null(backend) && is.function(backend$chol_solve_factor)) {
      resident_value <- .amatrix_amchol_resident_solve(factor, B_mat)
      if (!is.null(resident_value)) {
        resident_value
      } else {
        tryCatch(
          backend$chol_solve_factor(R, B_mat),
          error = function(e) {
            if (is.function(backend$solve_triangular_factor)) {
              z <- backend$solve_triangular_factor(R, B_mat, lower = FALSE, transpose = TRUE)
              backend$solve_triangular_factor(R, z, lower = FALSE, transpose = FALSE)
            } else {
              z <- forwardsolve(t(R), B_mat)
              backsolve(R, z)
            }
          }
        )
      }
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
  scalar_out <- is.vector(B) || (is.matrix(B) && ncol(B) == 1L)
  B_mat <- if (is.vector(B)) matrix(B, ncol = 1L) else as.matrix(B)
  x <- NULL

  if (inherits(R, "amChol")) {
    x <- .amatrix_amchol_resident_triangular_solve(R, B_mat, lower = lower, transpose = FALSE)
    if (is.null(x)) {
      backend <- .amatrix_amchol_backend(R)
      if (!is.null(backend) && is.function(backend$solve_triangular_factor)) {
        x <- tryCatch(
          backend$solve_triangular_factor(R@factor, B_mat, lower = lower, transpose = FALSE),
          error = function(e) NULL
        )
      }
    }
    R_mat <- R@factor
  } else {
    R_mat <- as.matrix(R)
    if (inherits(R, "adgeMatrix")) {
      backend_name <- .amatrix_live_resident_backend(R)
      if (!is.null(backend_name)) {
        x <- .amatrix_resident_triangular_solve(R, B_mat, backend_name, lower = lower, transpose = FALSE)
      }
      if (is.null(x) && isTRUE(R@precision == "fast") && nzchar(R@preferred_backend) && R@preferred_backend != "cpu") {
        backend <- tryCatch(.amatrix_get_backend(R@preferred_backend), error = function(e) NULL)
        if (!is.null(backend) && is.function(backend$solve_triangular_factor)) {
          x <- tryCatch(
            backend$solve_triangular_factor(R_mat, B_mat, lower = lower, transpose = FALSE),
            error = function(e) NULL
          )
        }
      }
    }
  }

  if (is.null(x)) {
    x <- if (isTRUE(lower)) {
      forwardsolve(R_mat, B_mat)
    } else {
      backsolve(R_mat, B_mat)
    }
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
