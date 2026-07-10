# Track 6 classed-condition helper. Uses base R errorCondition() to attach a
# condition class so callers can expect_error(..., class = "amatrix_bad_arg").
.amatrix_abort_bad_arg <- function(msg) {
  stop(errorCondition(msg, class = "amatrix_bad_arg", call = sys.call(-1L)))
}

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

.amatrix_amchol_dim <- function(factor) {
  stopifnot(inherits(factor, "amChol"))

  if (length(factor@factor) > 0L) {
    return(dim(factor@factor))
  }
  if (inherits(factor@factor_obj, "aMatrix")) {
    return(dim(factor@factor_obj))
  }
  c(0L, 0L)
}

.amatrix_amchol_factor_matrix <- function(factor) {
  stopifnot(inherits(factor, "amChol"))

  if (length(factor@factor) > 0L) {
    return(factor@factor)
  }
  if (!inherits(factor@factor_obj, "aMatrix")) {
    return(matrix(numeric(0), 0L, 0L))
  }

  mat <- as.matrix(factor@factor_obj)
  cache_key <- paste0("chol:", factor@source_id)
  cached <- .amatrix_cache_get(cache_key)
  if (inherits(cached, "amChol") && identical(cached@source_id, factor@source_id)) {
    cached@factor <- mat
    .amatrix_cache_set(cache_key, cached)
  }
  mat
}

#' Cholesky factorization result
#'
#' Stores the upper-triangular Cholesky factor \code{R} of a
#' symmetric positive-definite \code{adgeMatrix}, as returned by
#' \code{\link{chol_factor}}. When the factor is resident on a GPU
#' backend the host-side \code{@factor} slot may be an empty
#' zero-row matrix; use \code{\link{chol_solve}} rather than
#' accessing slots directly.
#'
#' @slot factor Numeric matrix; the upper-triangular factor \code{R}
#'   such that \code{t(R) \%*\% R} equals the source matrix.
#'   May be \code{matrix(numeric(0), 0, 0)} when the factor lives
#'   only on the device.
#' @slot factor_obj Either an \code{adgeMatrix} holding the
#'   GPU-resident factor, or \code{NULL} for CPU-only factors.
#' @slot source_id Character string; the \code{object_id} of the
#'   source \code{adgeMatrix}.
#' @slot precision Character string; \code{"strict"} or \code{"fast"}.
#' @slot backend Character string; the backend that computed the
#'   factorization.
#'
#' @exportClass amChol
#' @seealso \code{\link{chol_factor}}, \code{\link{chol_solve}},
#'   \code{\link{chol_logdet}}
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

#' @noRd
setMethod("show", "amChol", function(object) {
  dims <- .amatrix_amchol_dim(object)
  cat(sprintf(
    "amChol [%dx%d | %s | source: %s]\n",
    dims[[1L]],
    dims[[2L]],
    object@precision,
    object@source_id
  ))
  invisible(object)
})

#' Materialize the dense upper Cholesky factor of an amChol
#'
#' Returns the factor as a base R matrix regardless of where it currently
#' lives: the host copy when present, otherwise the resident/deferred factor
#' object is downloaded and cached.
#'
#' @param x An \code{\linkS4class{amChol}} object.
#' @param ... Ignored.
#' @return A base R numeric matrix containing the upper-triangular factor.
#' @method as.matrix amChol
#' @export
as.matrix.amChol <- function(x, ...) .amatrix_amchol_factor_matrix(x)

#' Compute the Cholesky factorization of an adgeMatrix
#'
#' Computes the upper-triangular Cholesky factor \code{R} of a
#' symmetric positive-definite \code{adgeMatrix} \code{X} such that
#' \code{t(R) \%*\% R == X}. Results are cached by \code{object_id};
#' repeated calls with the same object return the cached factor.
#'
#' @param X An \code{adgeMatrix} that is symmetric positive definite.
#'
#' @return An \code{\linkS4class{amChol}} object.
#'
#' @examples
#' m <- crossprod(matrix(rnorm(16), 4, 4)) + diag(4)
#' A <- adgeMatrix(m)
#' fac <- chol_factor(A)
#' fac
#'
#' @export
chol_factor <- function(X) {
  if (!inherits(X, "adgeMatrix")) {
    .amatrix_abort_bad_arg("X must be an adgeMatrix (symmetric positive definite)")
  }

  # Reject NA/NaN inputs before dispatch. LAPACK dpotrf is platform-fragile
  # on non-finite input: some builds (e.g. Linux reference LAPACK) return
  # info == 0 with a silent NaN factor instead of failing, while others
  # (e.g. macOS Accelerate) error. anyNA() is cheap and is TRUE for both NA
  # and NaN but FALSE for +/-Inf, so genuine Inf factorizations still flow
  # through and return a non-finite factor as before. Skip the check for
  # deferred objects, whose @x is a placeholder until first host access.
  if (!isTRUE(X@finalizer_env$host_deferred) && anyNA(X@x)) {
    .amatrix_abort_bad_arg(
      "X contains NA/NaN values; chol_factor() requires a finite symmetric positive-definite matrix"
    )
  }

  cache_key <- paste0("chol:", X@object_id)
  cached <- .amatrix_cache_get(cache_key)
  if (!is.null(cached)) {
    return(cached)
  }

  factor_value <- am_chol(X)
  backend_name <- amatrix_backend_plan(X, "chol")$chosen
  factor_backend <- if (inherits(factor_value, "aMatrix")) .amatrix_live_resident_backend(factor_value) else NULL
  R <- if (!is.null(factor_backend)) {
    matrix(numeric(0), 0L, 0L)
  } else {
    as.matrix(factor_value)
  }

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
  rhs_arg <- .amatrix_prepare_resident_arg(B, backend_name, promote_amatrix = TRUE)
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

.amatrix_amchol_result_dimnames <- function(factor_obj, rhs_template, out_ncol) {
  row_names <- NULL
  col_names <- NULL

  if (inherits(factor_obj, "aMatrix") && length(factor_obj@Dimnames) >= 1L) {
    row_names <- factor_obj@Dimnames[[1L]]
  }

  if (inherits(rhs_template, "aMatrix") && length(rhs_template@Dimnames) >= 2L) {
    col_names <- rhs_template@Dimnames[[2L]]
  } else if (is.matrix(rhs_template)) {
    col_names <- colnames(rhs_template)
  }

  if (!is.null(col_names) && length(col_names) != out_ncol) {
    col_names <- NULL
  }

  list(row_names, col_names)
}

.amatrix_amchol_wrap_resident_result <- function(factor_obj, rhs_template, backend_name, resident_key, out_ncol) {
  stopifnot(inherits(factor_obj, "aMatrix"))

  template <- if (inherits(rhs_template, "adgeMatrix")) rhs_template else factor_obj
  out_dim <- c(as.integer(nrow(factor_obj)), as.integer(out_ncol))

  obj <- new_adgeMatrix_deferred(
    dim = out_dim,
    dimnames = .amatrix_amchol_result_dimnames(factor_obj, rhs_template, out_ncol),
    preferred_backend = backend_name,
    policy = template@policy,
    precision = template@precision
  )

  .amatrix_bind_resident(obj, backend_name, resident_key)
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

  backend <- tryCatch(.amatrix_get_backend(backend_name), error = function(e) NULL)
  if (is.null(backend)) {
    return(NULL)
  }
  if (is.function(backend$prefer_solve_triangular_resident) &&
      !isTRUE(backend$prefer_solve_triangular_resident(factor_obj, B_mat, lower = lower, transpose = transpose))) {
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

.amatrix_amchol_resident_solve <- function(factor, B_mat, rhs_template = B_mat) {
  factor_obj <- factor@factor_obj
  if (!inherits(factor_obj, "aMatrix")) {
    return(NULL)
  }

  backend_name <- .amatrix_live_resident_backend(factor_obj)
  if (is.null(backend_name)) {
    return(NULL)
  }

  backend <- tryCatch(.amatrix_get_backend(backend_name), error = function(e) NULL)
  if (is.null(backend) || !is.function(backend$resident_materialize)) {
    return(NULL)
  }
  if (is.function(backend$prefer_chol_solve_resident) &&
      !isTRUE(backend$prefer_chol_solve_resident(factor_obj, B_mat))) {
    return(NULL)
  }

  factor_arg <- .amatrix_prepare_resident_arg(factor_obj, backend_name, promote_amatrix = FALSE)
  rhs_arg <- .amatrix_prepare_resident_arg(B_mat, backend_name, promote_amatrix = TRUE)
  if (is.null(factor_arg) || is.null(rhs_arg)) {
    .amatrix_cleanup_temp_resident(list(rhs_arg), backend_name)
    return(NULL)
  }

  x_key <- .amatrix_next_resident_key(backend_name)
  result <- NULL
  keep_result <- FALSE
  return_resident <- inherits(rhs_template, "adgeMatrix") &&
    identical(rhs_template@preferred_backend, backend_name)

  on.exit({
    .amatrix_cleanup_temp_resident(list(rhs_arg), backend_name)
    if (!isTRUE(keep_result) && isTRUE(backend$resident_has(x_key))) {
      try(backend$resident_drop(x_key), silent = TRUE)
    }
  }, add = TRUE)

  if (is.function(backend$chol_solve_resident)) {
    ok <- tryCatch({
      backend$chol_solve_resident(
        factor_arg$key,
        rhs_arg$key,
        x_key,
        defer = TRUE
      )
      if (isTRUE(return_resident)) {
        keep_result <- TRUE
        result <- .amatrix_amchol_wrap_resident_result(
          factor_obj,
          rhs_template,
          backend_name,
          x_key,
          out_ncol = if (is.null(dim(B_mat))) 1L else ncol(B_mat)
        )
      } else {
        result <- backend$resident_materialize(x_key)
      }
      TRUE
    }, error = function(e) FALSE)
  } else {
    if (!is.function(backend$solve_triangular_resident)) {
      return(NULL)
    }

    z_key <- .amatrix_next_resident_key(backend_name)
    on.exit({
      if (isTRUE(backend$resident_has(z_key))) {
        try(backend$resident_drop(z_key), silent = TRUE)
      }
    }, add = TRUE)

    ok <- tryCatch({
      backend$solve_triangular_resident(
        factor_arg$key,
        rhs_arg$key,
        z_key,
        lower = FALSE,
        transpose = TRUE,
        defer = TRUE
      )
      backend$solve_triangular_resident(
        factor_arg$key,
        z_key,
        x_key,
        lower = FALSE,
        transpose = FALSE,
        defer = TRUE
      )
      if (isTRUE(return_resident)) {
        keep_result <- TRUE
        result <- .amatrix_amchol_wrap_resident_result(
          factor_obj,
          rhs_template,
          backend_name,
          x_key,
          out_ncol = if (is.null(dim(B_mat))) 1L else ncol(B_mat)
        )
      } else {
        result <- backend$resident_materialize(x_key)
      }
      TRUE
    }, error = function(e) FALSE)
  }

  if (!ok) {
    return(NULL)
  }

  result
}

.amatrix_rhs_batches_to_list <- function(B, arg_name = "B") {
  if (is.list(B)) {
    return(B)
  }

  if (is.array(B) && length(dim(B)) == 3L) {
    nb <- dim(B)[[3L]]
    return(lapply(seq_len(nb), function(idx) {
      out <- B[, , idx, drop = FALSE]
      dim(out) <- dim(out)[1:2]
      out
    }))
  }

  stop(arg_name, " must be a list of RHS objects or a 3-D array [n, k, batch]", call. = FALSE)
}

.amatrix_triangular_rhs_arg <- function(B) {
  if (is.vector(B)) {
    return(matrix(B, ncol = 1L))
  }
  B
}

#' Solve a linear system using a Cholesky factor
#'
#' Solves \code{X \%*\% x = B} where \code{X} is the symmetric
#' positive-definite matrix whose Cholesky factorization is stored in
#' \code{factor}. Dispatches to a GPU backend when the factor was
#' computed in \code{"fast"} precision and a device-resident factor
#' is available.
#'
#' @param factor An \code{\linkS4class{amChol}} object from
#'   \code{\link{chol_factor}}.
#' @param B Numeric vector or matrix; the right-hand side. The number
#'   of rows must equal the dimension of the factor.
#'
#' @return Numeric vector or matrix \code{x} satisfying
#'   \code{X \%*\% x == B}. Returns a vector when \code{B} is a
#'   vector.
#'
#' @examples
#' m <- crossprod(matrix(rnorm(16), 4, 4)) + diag(4)
#' A <- adgeMatrix(m)
#' fac <- chol_factor(A)
#' b <- rnorm(4)
#' x <- chol_solve(fac, b)
#'
#' @export
chol_solve <- function(factor, B) {
  if (!inherits(factor, "amChol")) {
    .amatrix_abort_bad_arg("factor must be an amChol object")
  }

  B_in <- B
  B_arg <- .amatrix_triangular_rhs_arg(B)
  B_mat <- NULL

  # GPU path: dispatch through the backend's chol_solve_factor when the factor
  # was computed in fast mode on a GPU-capable backend.
  x <- if (isTRUE(factor@precision == "fast") &&
           nzchar(factor@backend) && factor@backend != "cpu") {
    backend <- .amatrix_amchol_backend(factor)
    if (!is.null(backend) && is.function(backend$chol_solve_factor)) {
      resident_value <- .amatrix_amchol_resident_solve(factor, B_arg)
      if (!is.null(resident_value)) {
        resident_value
      } else {
        R <- .amatrix_amchol_factor_matrix(factor)
        B_mat <- as.matrix(.amatrix_host_arg(B_arg))
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
      R <- .amatrix_amchol_factor_matrix(factor)
      B_mat <- as.matrix(.amatrix_host_arg(B_arg))
      z <- forwardsolve(t(R), B_mat)
      backsolve(R, z)
    }
  } else {
    # CPU path: standard triangular solve
    R <- .amatrix_amchol_factor_matrix(factor)
    B_mat <- as.matrix(.amatrix_host_arg(B_arg))
    z <- forwardsolve(t(R), B_mat)
    backsolve(R, z)
  }

  if (is.vector(B_in) && ncol(x) == 1L) x <- drop(x)
  x
}

#' Solve many right-hand-side batches with one Cholesky factor
#'
#' This is the same operation as repeatedly calling
#' \code{chol_solve(factor, B[[i]])}, but it packs all RHS batches into one
#' wide solve and then splits the result.  BLAS/GPU backends generally amortize
#' launch and dispatch overhead much better on one wide RHS than on many small
#' independent solves.
#'
#' @param factor An \code{amChol} object from \code{\link{chol_factor}}.
#' @param B A list of RHS vectors/matrices, or a 3-D array \code{[n, k, batch]}.
#' @return A list of solution vectors/matrices.
#' @export
chol_solve_batches <- function(factor, B) {
  if (!inherits(factor, "amChol")) {
    .amatrix_abort_bad_arg("factor must be an amChol object")
  }

  rhs_list <- .amatrix_rhs_batches_to_list(B, "B")
  dims <- .amatrix_amchol_dim(factor)
  if (length(rhs_list) == 0L) {
    return(list())
  }

  rhs_mats <- vector("list", length(rhs_list))
  rhs_vector <- logical(length(rhs_list))
  rhs_ncols <- integer(length(rhs_list))

  for (idx in seq_along(rhs_list)) {
    rhs_vector[[idx]] <- is.vector(rhs_list[[idx]])
    rhs_arg <- .amatrix_triangular_rhs_arg(rhs_list[[idx]])
    if (nrow(rhs_arg) != dims[[1L]]) {
      stop("all RHS batches must have nrow equal to the factor dimension", call. = FALSE)
    }
    rhs_mat <- as.matrix(.amatrix_host_arg(rhs_arg))
    if (!is.double(rhs_mat)) {
      storage.mode(rhs_mat) <- "double"
    }
    rhs_mats[[idx]] <- rhs_mat
    rhs_ncols[[idx]] <- ncol(rhs_mat)
  }

  # Pack RHS batches into one wide solve.  BLAS/GPU triangular solves amortize
  # launch and dispatch overhead much better on one n x sum(k_i) RHS than on
  # many small independent calls.
  rhs_packed <- do.call(cbind, rhs_mats)
  solved <- chol_solve(factor, rhs_packed)
  solved_mat <- as.matrix(.amatrix_host_arg(solved))

  starts <- cumsum(c(1L, rhs_ncols[-length(rhs_ncols)]))
  Map(function(start, width, was_vector) {
    value <- solved_mat[, seq.int(start, length.out = width), drop = FALSE]
    if (isTRUE(was_vector) && width == 1L) {
      drop(value)
    } else {
      value
    }
  }, starts, rhs_ncols, rhs_vector)
}

#' Extract the diagonal of a Cholesky factor
#'
#' Returns the diagonal of the upper-triangular matrix \code{R} stored
#' in an \code{\linkS4class{amChol}} object.
#'
#' @param factor An \code{\linkS4class{amChol}} object from
#'   \code{\link{chol_factor}}.
#'
#' @return Numeric vector of length equal to the matrix dimension.
#'
#' @examples
#' m <- crossprod(matrix(rnorm(16), 4, 4)) + diag(4)
#' A <- adgeMatrix(m)
#' fac <- chol_factor(A)
#' chol_diag(fac)
#'
#' @export
chol_diag <- function(factor) {
  if (!inherits(factor, "amChol")) {
    .amatrix_abort_bad_arg("factor must be an amChol object")
  }
  diag(.amatrix_amchol_factor_matrix(factor))
}

#' Log-determinant from a Cholesky factor
#'
#' Computes \code{log(det(X))} from the Cholesky factor of a
#' symmetric positive-definite matrix \code{X} as
#' \code{2 * sum(log(diag(R)))}, which avoids forming the full
#' determinant and is numerically stable.
#'
#' @param factor An \code{\linkS4class{amChol}} object from
#'   \code{\link{chol_factor}}.
#'
#' @return Scalar double; the log-determinant of the source matrix.
#'
#' @examples
#' m <- crossprod(matrix(rnorm(16), 4, 4)) + diag(4)
#' A <- adgeMatrix(m)
#' fac <- chol_factor(A)
#' chol_logdet(fac)
#'
#' @export
chol_logdet <- function(factor) {
  if (!inherits(factor, "amChol")) {
    .amatrix_abort_bad_arg("factor must be an amChol object")
  }
  2 * sum(log(diag(.amatrix_amchol_factor_matrix(factor))))
}

#' Solve a triangular linear system
#'
#' Solves \code{R \%*\% x = B} (or \code{t(R) \%*\% x = B} when
#' \code{lower = TRUE}) for \code{x}, where \code{R} is a triangular
#' matrix. Dispatches to a GPU backend when \code{R} is an
#' \code{\linkS4class{amChol}} or \code{adgeMatrix} with a live
#' resident key and a capable backend.
#'
#' @param R An \code{\linkS4class{amChol}}, \code{adgeMatrix}, or
#'   numeric matrix holding the triangular factor. Upper triangular
#'   by default.
#' @param B Numeric vector or matrix; the right-hand side.
#' @param lower Logical scalar; \code{FALSE} (default) treats
#'   \code{R} as upper triangular, \code{TRUE} treats it as lower
#'   triangular.
#'
#' @return Numeric vector or matrix \code{x} satisfying
#'   \code{R \%*\% x == B} (or its transpose variant).
#'   Returns a vector when \code{B} is a vector or single-column
#'   matrix.
#'
#' @examples
#' R <- chol(crossprod(matrix(rnorm(16), 4, 4)) + diag(4))
#' b <- rnorm(4)
#' x <- solve_triangular(R, b)
#'
#' @export
solve_triangular <- function(R, B, lower = FALSE) {
  scalar_out <- is.vector(B) || (is.matrix(B) && ncol(B) == 1L)
  B_arg <- .amatrix_triangular_rhs_arg(B)
  B_mat <- NULL
  x <- NULL

  if (inherits(R, "amChol")) {
    x <- .amatrix_amchol_resident_triangular_solve(R, B_arg, lower = lower, transpose = FALSE)
    if (is.null(x)) {
      backend <- .amatrix_amchol_backend(R)
      if (!is.null(backend) && is.function(backend$solve_triangular_factor)) {
        R_mat <- .amatrix_amchol_factor_matrix(R)
        B_mat <- as.matrix(.amatrix_host_arg(B_arg))
        x <- tryCatch(
          backend$solve_triangular_factor(R_mat, B_mat, lower = lower, transpose = FALSE),
          error = function(e) NULL
        )
      }
    }
    R_mat <- .amatrix_amchol_factor_matrix(R)
  } else {
    R_mat <- as.matrix(R)
    if (inherits(R, "adgeMatrix")) {
      backend_name <- .amatrix_live_resident_backend(R)
      if (!is.null(backend_name)) {
        x <- .amatrix_resident_triangular_solve(R, B_arg, backend_name, lower = lower, transpose = FALSE)
      }
      if (is.null(x) && isTRUE(R@precision == "fast") && nzchar(R@preferred_backend) && R@preferred_backend != "cpu") {
        backend <- tryCatch(.amatrix_get_backend(R@preferred_backend), error = function(e) NULL)
        if (!is.null(backend) && is.function(backend$solve_triangular_factor)) {
          B_mat <- as.matrix(.amatrix_host_arg(B_arg))
          x <- tryCatch(
            backend$solve_triangular_factor(R_mat, B_mat, lower = lower, transpose = FALSE),
            error = function(e) NULL
          )
        }
      }
    }
  }

  if (is.null(x)) {
    if (is.null(B_mat)) {
      B_mat <- as.matrix(.amatrix_host_arg(B_arg))
    }
    x <- if (isTRUE(lower)) {
      forwardsolve(R_mat, B_mat)
    } else {
      backsolve(R_mat, B_mat)
    }
  }
  if (scalar_out && ncol(x) == 1L) drop(x) else x
}

#' Evaluate a quadratic form using a Cholesky factor
#'
#' Computes \code{t(v) \%*\% solve(X) \%*\% v} (for a vector
#' \code{v}) or \code{t(V) \%*\% solve(X) \%*\% V} (for a matrix
#' \code{V}) efficiently via the Cholesky factor of \code{X}, without
#' forming the inverse.
#'
#' @param factor An \code{\linkS4class{amChol}} object from
#'   \code{\link{chol_factor}}.
#' @param v Numeric vector or matrix. For a vector, the result is a
#'   scalar; for a matrix with \code{p} columns, the result is a
#'   \code{p x p} matrix.
#'
#' @return Scalar double (when \code{v} is a vector) or numeric matrix
#'   of dimensions \code{ncol(v) x ncol(v)} containing the quadratic
#'   form.
#'
#' @examples
#' m <- crossprod(matrix(rnorm(16), 4, 4)) + diag(4)
#' A <- adgeMatrix(m)
#' fac <- chol_factor(A)
#' v <- rnorm(4)
#' quad_form(fac, v)
#'
#' @export
quad_form <- function(factor, v) {
  if (!inherits(factor, "amChol")) {
    .amatrix_abort_bad_arg("factor must be an amChol object")
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

#' LU factorization result for general square matrices
#'
#' Stores the original square matrix for use with LAPACK's
#' \code{DGESV} routine. Unlike \code{\linkS4class{amChol}}, which
#' caches the explicit triangular factor, \code{amLU} retains
#' \code{A} and delegates factorization to \code{base::solve} at
#' solve time.
#'
#' @slot A Numeric square matrix; the original matrix passed to
#'   \code{\link{lu_factor}}.
#' @slot source_id Character string; the \code{object_id} of the
#'   source \code{adgeMatrix}, or \code{NA} for base matrices.
#' @slot precision Character string; \code{"strict"} or \code{"fast"},
#'   or \code{NA} for base matrices.
#' @slot backend Character string; the preferred backend of the source
#'   object, or \code{NA} for base matrices.
#'
#' @exportClass amLU
#' @seealso \code{\link{lu_factor}}, \code{\link{lu_solve}}
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

#' @noRd
setMethod("show", "amLU", function(object) {
  cat(sprintf(
    "amLU [%dx%d | %s | source: %s]\n",
    nrow(object@A), ncol(object@A),
    object@precision, object@source_id
  ))
  invisible(object)
})

#' Store a general square matrix for LU-based solving
#'
#' Wraps a square numeric matrix or \code{adgeMatrix} in an
#' \code{\linkS4class{amLU}} object. The actual LU decomposition is
#' performed by \code{base::solve} at solve time via LAPACK's
#' \code{DGESV}.
#'
#' @param A A square numeric matrix or \code{adgeMatrix}.
#'
#' @return An \code{\linkS4class{amLU}} object.
#'
#' @examples
#' m <- matrix(c(2, 1, 5, 3), nrow = 2)
#' fac <- lu_factor(m)
#' fac
#'
#' @export
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
    .amatrix_abort_bad_arg("A must be a square matrix")
  }
  src  <- if (inherits(A, "adgeMatrix")) A@object_id else NA_character_
  prec <- if (inherits(A, "adgeMatrix")) A@precision  else NA_character_
  be   <- if (inherits(A, "adgeMatrix")) A@preferred_backend else NA_character_
  new("amLU", A = A_mat, source_id = src, precision = prec, backend = be)
}

#' Solve a linear system using an LU factor
#'
#' Solves \code{A \%*\% x = B} where \code{A} is the square matrix
#' stored in \code{factor}, delegating to \code{base::solve} (LAPACK
#' \code{DGESV}).
#'
#' @param factor An \code{\linkS4class{amLU}} object from
#'   \code{\link{lu_factor}}.
#' @param B Numeric vector or matrix; the right-hand side. The number
#'   of rows must equal the dimension of \code{factor@A}.
#'
#' @return Numeric vector or matrix \code{x} satisfying
#'   \code{A \%*\% x == B}. Returns a vector when \code{B} is a
#'   vector or single-column matrix.
#'
#' @examples
#' m <- matrix(c(2, 1, 5, 3), nrow = 2)
#' fac <- lu_factor(m)
#' b <- c(1, 2)
#' lu_solve(fac, b)
#'
#' @export
lu_solve <- function(factor, B) {
  if (!inherits(factor, "amLU")) {
    .amatrix_abort_bad_arg("factor must be an amLU object")
  }
  scalar_out <- is.vector(B) || (is.matrix(B) && ncol(B) == 1L)
  B_mat <- if (is.vector(B)) matrix(B, ncol = 1L) else as.matrix(B)
  x <- base::solve(factor@A, B_mat)
  if (scalar_out && ncol(x) == 1L) drop(x) else x
}
