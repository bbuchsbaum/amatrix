setOldClass(c("amQR", "amDenseQR"))

.amatrix_qr_or <- function(lhs, rhs) {
  if (!is.null(lhs)) {
    return(lhs)
  }
  rhs
}

.amatrix_qr_state <- function(factor = NULL, factor_source = NULL, q = NULL, q_key = NULL, backend_ops = NULL, factor_builder = NULL) {
  state <- new.env(parent = emptyenv())
  state$factor <- factor
  state$factor_source <- factor_source
  state$factor_builder <- factor_builder
  state$q <- q
  state$q_key <- q_key
  state$backend_ops <- backend_ops
  if (!is.null(backend_ops)) {
    reg.finalizer(
      state,
      function(e) {
        backend_name <- get0("backend_ops", envir = e, inherits = FALSE)
        if (is.null(backend_name) || !nzchar(backend_name)) {
          return(invisible(NULL))
        }
        backend <- tryCatch(.amatrix_get_backend(backend_name), error = function(err) NULL)
        if (is.null(backend) || !is.function(backend$resident_drop)) {
          return(invisible(NULL))
        }
        keys <- character()
        key <- get0("q_key", envir = e, inherits = FALSE)
        if (!is.null(key)) {
          keys <- c(keys, as.character(key))
        }
        factor <- get0("factor", envir = e, inherits = FALSE)
        if (is.list(factor) && !is.null(factor$block_q_keys)) {
          keys <- c(keys, as.character(factor$block_q_keys))
        }
        if (is.list(factor) && !is.null(factor$top_q_key)) {
          keys <- c(keys, as.character(factor$top_q_key))
        }
        if (is.list(factor) && !is.null(factor$top_r_key)) {
          keys <- c(keys, as.character(factor$top_r_key))
        }
        if (is.list(factor) && !is.null(factor$r_stack_key)) {
          keys <- c(keys, as.character(factor$r_stack_key))
        }
        keys <- unique(keys[nzchar(keys)])
        for (key in keys) {
          if (is.function(backend$resident_has) && !isTRUE(backend$resident_has(key))) {
            next
          }
          try(backend$resident_drop(key), silent = TRUE)
        }
        invisible(NULL)
      },
      onexit = TRUE
    )
  }
  state
}

.amatrix_wrap_qr <- function(qr_obj, x, method = "cpu") {
  stopifnot(is.list(qr_obj))
  source_dim <- dim(.amatrix_host_arg(x))
  qr_factor <- qr_obj[["factor", exact = TRUE]]
  qr_factor_builder <- qr_obj[["factor_builder", exact = TRUE]]
  qr_factor_source <- qr_obj[["factor_source", exact = TRUE]]
  qr_payload <- qr_obj[["qr", exact = TRUE]]
  qr_representation <- qr_obj[["representation", exact = TRUE]]
  qr_rank <- qr_obj[["rank", exact = TRUE]]
  qr_pivot <- qr_obj[["pivot", exact = TRUE]]
  qr_r <- qr_obj[["r", exact = TRUE]]
  qr_q <- qr_obj[["q", exact = TRUE]]
  qr_q_key <- qr_obj[["q_key", exact = TRUE]]

  if (!is.null(qr_payload)) {
    qr_base <- if (inherits(qr_obj, "qr")) {
      qr_obj
    } else if (inherits(qr_payload, "qr")) {
      qr_payload
    } else {
      qr_obj
    }
    representation <- "base_qr"
    rank <- as.integer(.amatrix_qr_or(qr_rank, qr_base$rank))
    q_materialized <- FALSE
    r_materialized <- FALSE
    thin <- TRUE
    pivot <- .amatrix_qr_or(qr_pivot, qr_base$pivot)
    pivoted <- !is.null(pivot) && !identical(as.integer(pivot), seq_along(pivot))
    state <- .amatrix_qr_state(
      factor = .amatrix_qr_or(qr_factor, qr_base),
      factor_source = .amatrix_qr_or(qr_factor_source, "native")
    )
    wrapped_payload <- qr_base
  } else if (identical(qr_representation, "mlx_compact_qr")) {
    representation <- "mlx_compact_qr"
    rank <- as.integer(.amatrix_qr_or(qr_rank, .amatrix_qr_or(qr_factor$rank, if (!is.null(source_dim)) min(source_dim) else NA_integer_)))
    q_materialized <- FALSE
    r_materialized <- !is.null(qr_r)
    r_dim <- if (!is.null(qr_r)) dim(as.matrix(qr_r)) else if (!is.null(source_dim)) c(min(source_dim), source_dim[[2]]) else NULL
    thin <- !is.null(source_dim) && !is.null(r_dim) && identical(r_dim[[1]], min(source_dim))
    pivot <- if (!is.null(qr_pivot)) as.integer(qr_pivot) else NULL
    pivoted <- !is.null(pivot) && !identical(as.integer(pivot), seq_along(pivot))
    state <- .amatrix_qr_state(
      factor = qr_factor,
      factor_source = qr_factor_source,
      factor_builder = qr_factor_builder,
      q = NULL,
      q_key = qr_q_key,
      backend_ops = qr_obj[["backend_ops", exact = TRUE]]
    )
    wrapped_payload <- qr_obj
  } else if (!is.null(qr_q) && !is.null(qr_r)) {
    representation <- .amatrix_qr_or(qr_representation, "explicit_qr")
    rank <- .amatrix_explicit_qr_rank(qr_obj)
    q_materialized <- TRUE
    r_materialized <- TRUE
    q_dim <- dim(as.matrix(qr_q))
    thin <- !is.null(source_dim) && !is.null(q_dim) && identical(q_dim[[2]], min(source_dim))
    pivot <- if (!is.null(qr_pivot)) as.integer(qr_pivot) else NULL
    pivoted <- !is.null(pivot) && !identical(as.integer(pivot), seq_along(pivot))
    state <- .amatrix_qr_state(
      factor = qr_factor,
      factor_source = qr_factor_source,
      factor_builder = qr_factor_builder,
      q = qr_q,
      q_key = qr_q_key,
      backend_ops = qr_obj[["backend_ops", exact = TRUE]]
    )
    wrapped_payload <- qr_obj
  } else if (!is.null(qr_q_key) && !is.null(qr_r)) {
    representation <- .amatrix_qr_or(qr_representation, "explicit_qr")
    rank <- .amatrix_explicit_qr_rank(qr_obj)
    q_materialized <- FALSE
    r_materialized <- TRUE
    r_dim <- dim(as.matrix(qr_r))
    thin <- !is.null(source_dim) && !is.null(r_dim) && identical(r_dim[[1]], min(source_dim))
    pivot <- if (!is.null(qr_pivot)) as.integer(qr_pivot) else NULL
    pivoted <- !is.null(pivot) && !identical(as.integer(pivot), seq_along(pivot))
    state <- .amatrix_qr_state(
      factor = qr_factor,
      factor_source = qr_factor_source,
      factor_builder = qr_factor_builder,
      q = NULL,
      q_key = qr_q_key,
      backend_ops = qr_obj[["backend_ops", exact = TRUE]]
    )
    wrapped_payload <- qr_obj
  } else {
    stop("unsupported QR payload", call. = FALSE)
  }

  structure(
    list(
      qr = wrapped_payload,
      representation = representation,
      rank = rank,
      backend = if (inherits(x, "aMatrix")) x@preferred_backend else "cpu",
      precision = if (inherits(x, "aMatrix")) x@precision else amatrix_default_precision(),
      method = method,
      source_class = class(x)[[1]],
      source_dim = source_dim,
      thin = thin,
      pivoted = pivoted,
      pivot = pivot,
      state = state,
      q_materialized = q_materialized,
      r_materialized = r_materialized
    ),
    class = c("amQR", "amDenseQR")
  )
}

.amatrix_unwrap_qr <- function(qr) {
  if (inherits(qr, "amQR") && !is.null(qr$qr)) {
    return(qr$qr)
  }
  qr
}

.amatrix_qr_kind <- function(qr) {
  if (inherits(qr, "amQR")) {
    return(.amatrix_qr_or(qr$representation, "base_qr"))
  }
  if (is.list(qr) && !is.null(qr[["representation", exact = TRUE]])) {
    return(as.character(qr[["representation", exact = TRUE]]))
  }
  if (is.list(qr) && !is.null(qr[["qr", exact = TRUE]])) {
    return("base_qr")
  }
  if (is.list(qr) && (!is.null(qr[["q", exact = TRUE]]) || !is.null(qr[["q_key", exact = TRUE]])) && !is.null(qr[["r", exact = TRUE]])) {
    return("explicit_qr")
  }
  stop("unsupported QR payload", call. = FALSE)
}

.amatrix_qr_rank <- function(qr) {
  if (inherits(qr, "amQR") && !is.null(qr$rank)) {
    return(as.integer(qr$rank))
  }

  if (inherits(qr, "amDenseQR")) {
    return(as.integer(qr$rank))
  }

  payload <- .amatrix_unwrap_qr(qr)
  if (.amatrix_qr_kind(payload) == "base_qr") {
    return(as.integer(.amatrix_qr_or(payload[["rank", exact = TRUE]], payload[["qr", exact = TRUE]]$rank)))
  }

  .amatrix_explicit_qr_rank(payload)
}

.amatrix_wrap_sparse_qr <- function(qr_obj, x, method = "cpu") {
  stopifnot(inherits(qr_obj, "sparseQR"))

  pivot <- tryCatch(as.integer(qr_obj@q) + 1L, error = function(e) NULL)
  pivoted <- !is.null(pivot) && !identical(pivot, seq_along(pivot))
  source_dim <- dim(.amatrix_host_arg(x))
  rank <- as.integer(Matrix::rankMatrix(qr_obj@R)[1L])

  structure(
    list(
      qr = qr_obj,
      representation = "base_qr",
      rank = rank,
      backend = if (inherits(x, "aMatrix")) x@preferred_backend else "cpu",
      precision = if (inherits(x, "aMatrix")) x@precision else amatrix_default_precision(),
      method = method,
      source_class = class(x)[[1L]],
      source_dim = source_dim,
      thin = TRUE,
      pivoted = pivoted,
      pivot = pivot,
      state = .amatrix_qr_state(factor = qr_obj, factor_source = "native"),
      q_materialized = FALSE,
      r_materialized = FALSE
    ),
    class = "amQR"
  )
}

.amatrix_qr_source_dim <- function(qr) {
  if (inherits(qr, "amQR")) {
    return(qr$source_dim)
  }
  NULL
}

.amatrix_qr_thin <- function(qr) {
  if (inherits(qr, "amQR")) {
    return(isTRUE(qr$thin))
  }
  TRUE
}

.amatrix_qr_pivoted <- function(qr) {
  if (inherits(qr, "amQR")) {
    return(isTRUE(qr$pivoted))
  }
  FALSE
}

.amatrix_qr_pivot <- function(qr) {
  if (inherits(qr, "amQR")) {
    return(qr$pivot)
  }
  NULL
}

.amatrix_qr_backend_ops <- function(qr) {
  payload <- .amatrix_unwrap_qr(qr)
  if (is.list(payload) && !is.null(payload[["backend_ops", exact = TRUE]])) {
    return(as.character(payload[["backend_ops", exact = TRUE]]))
  }
  NULL
}

.amatrix_qr_helper_path <- function(qr) {
  if (identical(.amatrix_qr_kind(qr), "mlx_compact_qr")) {
    return("compact_mlx_factor")
  }
  backend_ops <- .amatrix_qr_backend_ops(qr)
  if (identical(.amatrix_qr_kind(qr), "explicit_qr") && identical(backend_ops, "mlx") && requireNamespace("amatrix.mlx", quietly = TRUE)) {
    helper_mode <- get(".amatrix_mlx_qr_helper_mode", envir = asNamespace("amatrix.mlx"), inherits = FALSE)()
    if (identical(helper_mode, "native")) {
      if (!is.null(.amatrix_qr_q_key(qr))) {
        return("native_resident_backend")
      }
      return("native_backend")
    }
  }
  if (identical(.amatrix_qr_kind(qr), "explicit_qr") &&
      identical(backend_ops, "opencl") &&
      requireNamespace("amatrix.opencl", quietly = TRUE) &&
      !is.null(.amatrix_qr_q_key(qr))) {
    return("native_resident_backend")
  }
  "compact_factor"
}

.amatrix_qr_compact_available <- function(qr) {
  if (!inherits(qr, "amQR")) {
    return(FALSE)
  }
  if (!is.null(qr$state$factor) || is.function(qr$state$factor_builder)) {
    return(TRUE)
  }
  identical(.amatrix_qr_kind(qr), "explicit_qr")
}

.amatrix_qr_compact_materialized <- function(qr) {
  inherits(qr, "amQR") && !is.null(qr$state$factor)
}

.amatrix_qr_q_materialized <- function(qr) {
  payload <- .amatrix_unwrap_qr(qr)
  inherits(qr, "amQR") && (!is.null(qr$state$q) || !is.null(payload[["q", exact = TRUE]]))
}

.amatrix_qr_r_materialized <- function(qr) {
  payload <- .amatrix_unwrap_qr(qr)
  inherits(qr, "amQR") && !is.null(payload[["r", exact = TRUE]])
}

.amatrix_qr_q_key <- function(qr) {
  if (!inherits(qr, "amQR")) {
    return(NULL)
  }
  .amatrix_qr_or(qr$state$q_key, .amatrix_unwrap_qr(qr)[["q_key", exact = TRUE]])
}

.amatrix_qr_reconstruct_source <- function(qr) {
  if (!identical(.amatrix_qr_kind(qr), "explicit_qr")) {
    stop("source reconstruction is only supported for explicit_qr payloads", call. = FALSE)
  }
  .amatrix_explicit_qr_q(qr) %*% .amatrix_explicit_qr_r(qr)
}

.amatrix_qr_factor <- function(qr) {
  if (!inherits(qr, "amQR")) {
    stop("qr must inherit from amQR", call. = FALSE)
  }

  factor <- qr$state$factor
  if (!is.null(factor)) {
    return(factor)
  }

  factor_builder <- qr$state$factor_builder
  if (is.function(factor_builder)) {
    factor <- factor_builder()
    qr$state$factor <- factor
    qr$state$factor_builder <- NULL
    return(factor)
  }

  if (!identical(.amatrix_qr_kind(qr), "explicit_qr")) {
    stop("compact QR factor is unavailable", call. = FALSE)
  }

  source <- .amatrix_qr_reconstruct_source(qr)
  factor <- base::qr(source)
  qr$state$factor <- factor
  qr$state$factor_source <- "reconstructed"
  factor
}

#' Inspect an amQR factorization object
#'
#' Returns a named list of metadata fields describing an \code{amQR}
#' factor produced by \code{am_qr()}.
#'
#' @param qr An \code{amQR} object.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{rank}{Integer effective rank of the factored matrix.}
#'     \item{dim}{Integer vector of length 2: \code{c(nrow, ncol)} of
#'       the source matrix.}
#'     \item{thin}{Logical; \code{TRUE} when a thin (economy) QR was
#'       computed.}
#'     \item{pivoted}{Logical; \code{TRUE} when column pivoting was
#'       used.}
#'     \item{pivot}{Integer permutation vector, or \code{NULL} when
#'       unpivoted.}
#'     \item{representation}{Character string describing the internal
#'       storage format.}
#'     \item{backend_ops}{Character string naming the backend that owns
#'       any resident buffers, or \code{NULL}.}
#'     \item{backend}{Character string: the preferred backend.}
#'     \item{precision}{Character string: \code{"strict"} or
#'       \code{"fast"}.}
#'     \item{method}{Character string: QR algorithm used.}
#'     \item{compact_factor_available}{Logical.}
#'     \item{compact_factor_source}{Character string or \code{NULL}.}
#'     \item{compact_factor_materialized}{Logical.}
#'     \item{q_materialized}{Logical.}
#'     \item{r_materialized}{Logical.}
#'   }
#'
#' @examples
#' X <- adgeMatrix(matrix(rnorm(30), nrow = 6))
#' qf <- am_qr(X)
#' info <- qr_info(qf)
#' info$rank
#'
#' @seealso \code{\link{am_qr}}, \code{\link{qr_downdate}}
#' @export
qr_info <- function(qr) {
  if (!inherits(qr, "amQR")) {
    stop("qr must inherit from amQR", call. = FALSE)
  }

  list(
    rank = .amatrix_qr_rank(qr),
    dim = .amatrix_qr_source_dim(qr),
    thin = .amatrix_qr_thin(qr),
    pivoted = .amatrix_qr_pivoted(qr),
    pivot = .amatrix_qr_pivot(qr),
    representation = .amatrix_qr_kind(qr),
    helper_path = .amatrix_qr_helper_path(qr),
    backend_ops = .amatrix_qr_backend_ops(qr),
    backend = qr$backend,
    precision = qr$precision,
    method = qr$method,
    compact_factor_available = .amatrix_qr_compact_available(qr),
    compact_factor_source = qr$state$factor_source,
    compact_factor_materialized = .amatrix_qr_compact_materialized(qr),
    q_materialized = .amatrix_qr_q_materialized(qr),
    r_materialized = .amatrix_qr_r_materialized(qr)
  )
}

.amatrix_explicit_qr_native_mlx <- function(fun, ...) {
  get(fun, envir = asNamespace("amatrix.mlx"), inherits = FALSE)(...)
}

.amatrix_explicit_qr_native_opencl <- function(fun, ...) {
  get(fun, envir = asNamespace("amatrix.opencl"), inherits = FALSE)(...)
}

.amatrix_explicit_qr_resident_materialize <- function(qr, key) {
  backend_name <- .amatrix_qr_backend_ops(qr)
  if (is.null(backend_name) || !nzchar(backend_name)) {
    return(NULL)
  }

  backend <- tryCatch(.amatrix_get_backend(backend_name), error = function(e) NULL)
  if (!is.null(backend) && is.function(backend$resident_materialize)) {
    value <- tryCatch(backend$resident_materialize(key), error = function(e) NULL)
    if (!is.null(value)) {
      return(value)
    }
  }

  fallback <- switch(
    backend_name,
    mlx = c("amatrix.mlx", "amatrix_mlx_resident_materialize"),
    opencl = c("amatrix.opencl", "amatrix_opencl_resident_materialize"),
    arrayfire = c("amatrix.arrayfire", "amatrix_arrayfire_resident_materialize"),
    NULL
  )
  if (is.null(fallback) || !requireNamespace(fallback[[1L]], quietly = TRUE)) {
    return(NULL)
  }

  fun <- get0(fallback[[2L]], envir = asNamespace(fallback[[1L]]), inherits = FALSE)
  if (!is.function(fun)) {
    return(NULL)
  }

  tryCatch(fun(key), error = function(e) NULL)
}

.amatrix_explicit_qr_rank <- function(qr_obj) {
  .amatrix_explicit_qr_rank_tol(qr_obj, tol = NULL)
}

.amatrix_explicit_qr_rank_tol <- function(qr_obj, tol = NULL) {
  if (!is.null(qr_obj[["rank", exact = TRUE]])) {
    return(as.integer(qr_obj[["rank", exact = TRUE]]))
  }

  r_mat <- as.matrix(qr_obj[["r", exact = TRUE]])
  diag_len <- min(dim(r_mat))
  if (diag_len == 0L) {
    return(0L)
  }

  diag_abs <- abs(diag(r_mat[seq_len(diag_len), seq_len(diag_len), drop = FALSE]))
  scale <- if (length(diag_abs) == 0L) 0 else max(diag_abs)
  tol_eff <- if (is.null(tol)) max(dim(r_mat)) * .Machine$double.eps * scale else as.double(tol)
  as.integer(sum(diag_abs > tol_eff))
}

.amatrix_explicit_qr_pivot <- function(qr_obj) {
  pivot <- qr_obj[["pivot", exact = TRUE]]
  if (is.null(pivot)) {
    return(NULL)
  }
  as.integer(pivot)
}

.amatrix_explicit_qr_unpivot <- function(coef, qr_obj, fill = NA_real_) {
  pivot <- .amatrix_explicit_qr_pivot(qr_obj)
  if (is.null(pivot) || identical(pivot, seq_len(length(pivot)))) {
    return(coef)
  }

  out <- matrix(
    fill,
    nrow = length(pivot),
    ncol = ncol(coef),
    dimnames = dimnames(coef)
  )
  out[pivot, ] <- coef
  out
}

.amatrix_explicit_qr_q <- function(qr) {
  payload <- .amatrix_unwrap_qr(qr)
  q_payload <- payload[["q", exact = TRUE]]
  if (!is.null(q_payload)) {
    return(as.matrix(q_payload))
  }
  if (inherits(qr, "amQR") && !is.null(qr$state$q)) {
    return(as.matrix(qr$state$q))
  }
  q_key <- .amatrix_qr_q_key(qr)
  if (!is.null(q_key)) {
    q_mat <- if (identical(.amatrix_qr_backend_ops(qr), "mlx")) {
      .amatrix_explicit_qr_native_mlx("amatrix_mlx_resident_materialize", q_key)
    } else {
      .amatrix_explicit_qr_resident_materialize(qr, q_key)
    }
    if (is.null(q_mat)) {
      stop("explicit QR resident q could not be materialized", call. = FALSE)
    }
    if (inherits(qr, "amQR")) {
      qr$state$q <- q_mat
    }
    return(as.matrix(q_mat))
  }
  stop("explicit QR payload does not contain q", call. = FALSE)
}

.amatrix_explicit_qr_r <- function(qr) {
  payload <- .amatrix_unwrap_qr(qr)
  as.matrix(payload[["r", exact = TRUE]])
}

.amatrix_base_qr_q <- function(qr, complete = FALSE) {
  payload <- .amatrix_unwrap_qr(qr)
  qr.Q(payload, complete = complete)
}

.amatrix_base_qr_r <- function(qr, complete = FALSE) {
  payload <- .amatrix_unwrap_qr(qr)
  qr.R(payload, complete = complete)
}

.amatrix_qr_q <- function(qr, complete = FALSE) {
  if (identical(.amatrix_qr_kind(qr), "mlx_compact_qr")) {
    q_key <- .amatrix_qr_q_key(qr)
    if (!isTRUE(complete) && !is.null(q_key) && identical(.amatrix_qr_backend_ops(qr), "mlx")) {
      return(.amatrix_explicit_qr_native_mlx("amatrix_mlx_resident_materialize", q_key))
    }
    factor_source <- .amatrix_qr_or(qr$state$factor_source, NULL)
    if (identical(factor_source, "tsqr_blocked")) {
      return(.amatrix_explicit_qr_native_mlx("amatrix_mlx_tsqr_q", .amatrix_qr_factor(qr), complete = complete))
    }
    return(qr.Q(.amatrix_qr_factor(qr), complete = complete))
  }
  if (identical(.amatrix_qr_kind(qr), "explicit_qr")) {
    return(.amatrix_explicit_qr_q(qr))
  }
  .amatrix_base_qr_q(qr, complete = complete)
}

.amatrix_qr_r <- function(qr, complete = FALSE) {
  if (identical(.amatrix_qr_kind(qr), "mlx_compact_qr")) {
    payload <- .amatrix_unwrap_qr(qr)
    if (!is.null(payload[["r", exact = TRUE]])) {
      return(as.matrix(payload[["r", exact = TRUE]]))
    }
    if (!is.null(payload[["r_key", exact = TRUE]]) && identical(.amatrix_qr_backend_ops(qr), "mlx")) {
      return(.amatrix_explicit_qr_native_mlx("amatrix_mlx_resident_materialize", payload[["r_key", exact = TRUE]]))
    }
    return(qr.R(.amatrix_qr_factor(qr), complete = complete))
  }
  if (identical(.amatrix_qr_kind(qr), "explicit_qr")) {
    return(.amatrix_explicit_qr_r(qr))
  }
  .amatrix_base_qr_r(qr, complete = complete)
}

.amatrix_mlx_compact_qr_qty <- function(qr, y) {
  if (identical(.amatrix_qr_or(qr$state$factor_source, NULL), "tsqr_blocked")) {
    return(.amatrix_explicit_qr_native_mlx("amatrix_mlx_tsqr_qty", .amatrix_qr_factor(qr), y))
  }
  qr.qty(.amatrix_qr_factor(qr), y)
}

.amatrix_mlx_compact_qr_qy <- function(qr, y) {
  if (identical(.amatrix_qr_or(qr$state$factor_source, NULL), "tsqr_blocked")) {
    return(.amatrix_explicit_qr_native_mlx("amatrix_mlx_tsqr_qy", .amatrix_qr_factor(qr), y))
  }
  qr.qy(.amatrix_qr_factor(qr), y)
}

.amatrix_mlx_compact_qr_coef <- function(qr, y) {
  if (identical(.amatrix_qr_or(qr$state$factor_source, NULL), "tsqr_blocked")) {
    return(.amatrix_explicit_qr_native_mlx("amatrix_mlx_tsqr_coef", .amatrix_qr_factor(qr), y))
  }
  qr.coef(.amatrix_qr_factor(qr), y)
}

.amatrix_mlx_compact_qr_solve <- function(qr, b = NULL, tol = 1e-07) {
  if (identical(.amatrix_qr_or(qr$state$factor_source, NULL), "tsqr_blocked")) {
    return(.amatrix_explicit_qr_native_mlx("amatrix_mlx_tsqr_solve", .amatrix_qr_factor(qr), b, tol = tol))
  }
  if (is.null(b)) {
    return(qr.solve(.amatrix_qr_factor(qr), tol = tol))
  }
  qr.solve(.amatrix_qr_factor(qr), b = b, tol = tol)
}

.amatrix_mlx_compact_qr_fitted <- function(qr, y, k = NULL) {
  if (identical(.amatrix_qr_or(qr$state$factor_source, NULL), "tsqr_blocked")) {
    return(.amatrix_explicit_qr_native_mlx("amatrix_mlx_tsqr_fitted", .amatrix_qr_factor(qr), y, k = k))
  }
  qr.fitted(.amatrix_qr_factor(qr), y, k = k)
}

.amatrix_mlx_compact_qr_resid <- function(qr, y) {
  if (identical(.amatrix_qr_or(qr$state$factor_source, NULL), "tsqr_blocked")) {
    return(.amatrix_explicit_qr_native_mlx("amatrix_mlx_tsqr_resid", .amatrix_qr_factor(qr), y))
  }
  qr.resid(.amatrix_qr_factor(qr), y)
}

.amatrix_explicit_qr_qty <- function(qr, y) {
  helper_path <- .amatrix_qr_helper_path(qr)
  if (helper_path %in% c("native_backend", "native_resident_backend") && identical(.amatrix_qr_backend_ops(qr), "mlx")) {
    q_key <- .amatrix_qr_q_key(qr)
    if (!is.null(q_key)) {
      return(.amatrix_explicit_qr_native_mlx("amatrix_mlx_qr_qty_key", q_key, y))
    }
    return(.amatrix_explicit_qr_native_mlx("amatrix_mlx_qr_qty", .amatrix_explicit_qr_q(qr), y))
  }
  if (helper_path %in% c("native_backend", "native_resident_backend") && identical(.amatrix_qr_backend_ops(qr), "opencl")) {
    q_key <- .amatrix_qr_q_key(qr)
    if (!is.null(q_key)) {
      return(.amatrix_explicit_qr_native_opencl("amatrix_opencl_qr_qty_key", q_key, y))
    }
  }
  crossprod(.amatrix_explicit_qr_q(qr), y)
}

.amatrix_explicit_qr_qy <- function(qr, y) {
  helper_path <- .amatrix_qr_helper_path(qr)
  if (helper_path %in% c("native_backend", "native_resident_backend") && identical(.amatrix_qr_backend_ops(qr), "mlx")) {
    q_key <- .amatrix_qr_q_key(qr)
    if (!is.null(q_key)) {
      return(.amatrix_explicit_qr_native_mlx("amatrix_mlx_qr_qy_key", q_key, y))
    }
    return(.amatrix_explicit_qr_native_mlx("amatrix_mlx_qr_qy", .amatrix_explicit_qr_q(qr), y))
  }
  if (helper_path %in% c("native_backend", "native_resident_backend") && identical(.amatrix_qr_backend_ops(qr), "opencl")) {
    q_key <- .amatrix_qr_q_key(qr)
    if (!is.null(q_key)) {
      return(.amatrix_explicit_qr_native_opencl("amatrix_opencl_qr_qy_key", q_key, y))
    }
  }
  .amatrix_explicit_qr_q(qr) %*% y
}

.amatrix_explicit_qr_coef <- function(qr, y) {
  helper_path <- .amatrix_qr_helper_path(qr)
  if (helper_path %in% c("native_backend", "native_resident_backend") && identical(.amatrix_qr_backend_ops(qr), "mlx")) {
    q_key <- .amatrix_qr_q_key(qr)
    r_mat <- .amatrix_explicit_qr_r(qr)
    rank <- .amatrix_qr_rank(qr)
    if (!is.null(q_key) && identical(rank, ncol(r_mat))) {
      return(.amatrix_explicit_qr_native_mlx("amatrix_mlx_qr_coef_key", q_key, r_mat, y))
    }
    return(.amatrix_explicit_qr_native_mlx(
      "amatrix_mlx_qr_coef",
      .amatrix_explicit_qr_q(qr),
      r_mat,
      y,
      rank = rank
    ))
  }

  r_mat <- .amatrix_explicit_qr_r(qr)
  p <- ncol(r_mat)
  rank <- .amatrix_qr_rank(qr)
  qty <- .amatrix_explicit_qr_qty(qr, y)
  coef <- matrix(NA_real_, nrow = p, ncol = ncol(y))

  if (rank > 0L) {
    r_top <- r_mat[seq_len(rank), seq_len(rank), drop = FALSE]
    qty_top <- qty[seq_len(rank), , drop = FALSE]
    coef[seq_len(rank), ] <- backsolve(r_top, qty_top)
  }

  .amatrix_explicit_qr_unpivot(coef, .amatrix_unwrap_qr(qr))
}

.amatrix_explicit_qr_solve <- function(qr, b = NULL, tol = 1e-07) {
  source_dim <- .amatrix_qr_source_dim(qr)
  r_mat <- .amatrix_explicit_qr_r(qr)
  p <- ncol(r_mat)
  rank <- .amatrix_explicit_qr_rank_tol(.amatrix_unwrap_qr(qr), tol = tol)

  if (is.null(b)) {
    if (is.null(source_dim) || source_dim[[1]] != source_dim[[2]]) {
      stop("only square matrices can be inverted", call. = FALSE)
    }
    b <- diag(p)
  } else {
    b <- as.matrix(b)
  }

  if (rank < p) {
    stop("singular matrix 'a' in solve", call. = FALSE)
  }

  .amatrix_explicit_qr_coef(qr, b)
}

.amatrix_explicit_qr_fitted <- function(qr, y, k = NULL) {
  helper_path <- .amatrix_qr_helper_path(qr)
  if (helper_path %in% c("native_backend", "native_resident_backend") && identical(.amatrix_qr_backend_ops(qr), "mlx")) {
    q_key <- .amatrix_qr_q_key(qr)
    r_mat <- .amatrix_explicit_qr_r(qr)
    if (!is.null(q_key) && !is.null(k) && identical(k, ncol(r_mat))) {
      qty <- .amatrix_explicit_qr_native_mlx("amatrix_mlx_qr_qty_key", q_key, y)
      return(.amatrix_explicit_qr_native_mlx("amatrix_mlx_qr_qy_key", q_key, qty))
    }
    return(.amatrix_explicit_qr_native_mlx(
      "amatrix_mlx_qr_fitted",
      .amatrix_explicit_qr_q(qr),
      y,
      rank = if (is.null(k)) .amatrix_qr_rank(qr) else k
    ))
  }
  if (helper_path %in% c("native_backend", "native_resident_backend") && identical(.amatrix_qr_backend_ops(qr), "opencl")) {
    q_key <- .amatrix_qr_q_key(qr)
    r_mat <- .amatrix_explicit_qr_r(qr)
    if (!is.null(q_key) && !is.null(k) && identical(k, ncol(r_mat))) {
      return(.amatrix_explicit_qr_native_opencl("amatrix_opencl_qr_fitted_key", q_key, y))
    }
  }

  q_mat <- .amatrix_explicit_qr_q(qr)
  qty <- .amatrix_explicit_qr_qty(qr, y)
  if (is.null(k)) {
    k <- .amatrix_qr_rank(qr)
  }

  if (k <= 0L) {
    return(matrix(0, nrow = nrow(q_mat), ncol = ncol(y)))
  }

  q_mat[, seq_len(k), drop = FALSE] %*% qty[seq_len(k), , drop = FALSE]
}

.amatrix_qr_solve_value <- function(a, b = NULL, tol = 1e-07) {
  if (identical(.amatrix_qr_kind(a), "mlx_compact_qr")) {
    return(.amatrix_mlx_compact_qr_solve(a, b = b, tol = tol))
  }
  if (identical(.amatrix_qr_kind(a), "explicit_qr")) {
    if (is.null(b)) {
      if (identical(.amatrix_qr_helper_path(a), "compact_factor") && .amatrix_qr_compact_available(a)) {
        return(qr.solve(.amatrix_qr_factor(a), tol = tol))
      }
      return(.amatrix_explicit_qr_solve(a, tol = tol))
    }

    if (identical(.amatrix_qr_helper_path(a), "compact_factor") && .amatrix_qr_compact_available(a)) {
      return(qr.solve(.amatrix_qr_factor(a), b = b, tol = tol))
    }
    return(.amatrix_explicit_qr_solve(a, b = b, tol = tol))
  }

  if (is.null(b)) {
    return(qr.solve(.amatrix_qr_factor(a), tol = tol))
  }

  qr.solve(.amatrix_qr_factor(a), b = b, tol = tol)
}

.amatrix_qr_coef_value <- function(qr, y) {
  if (identical(.amatrix_qr_kind(qr), "mlx_compact_qr")) {
    return(.amatrix_mlx_compact_qr_coef(qr, y))
  }
  if (identical(.amatrix_qr_kind(qr), "explicit_qr")) {
    if (identical(.amatrix_qr_helper_path(qr), "compact_factor") && .amatrix_qr_compact_available(qr)) {
      return(qr.coef(.amatrix_qr_factor(qr), y))
    }
    return(.amatrix_explicit_qr_coef(qr, y))
  }

  qr.coef(.amatrix_qr_factor(qr), y)
}

.amatrix_qr_fitted_value <- function(qr, y, k = NULL) {
  if (is.null(k)) {
    k <- .amatrix_qr_rank(qr)
  }

  if (identical(.amatrix_qr_kind(qr), "mlx_compact_qr")) {
    return(.amatrix_mlx_compact_qr_fitted(qr, y, k = k))
  }

  if (identical(.amatrix_qr_kind(qr), "explicit_qr")) {
    if (identical(.amatrix_qr_helper_path(qr), "compact_factor") && .amatrix_qr_compact_available(qr)) {
      return(qr.fitted(.amatrix_qr_factor(qr), y, k = k))
    }
    return(.amatrix_explicit_qr_fitted(qr, y, k = k))
  }

  qr.fitted(.amatrix_qr_factor(qr), y, k = k)
}

.amatrix_qr_resid_value <- function(qr, y) {
  if (identical(.amatrix_qr_kind(qr), "mlx_compact_qr")) {
    return(.amatrix_mlx_compact_qr_resid(qr, y))
  }
  if (identical(.amatrix_qr_kind(qr), "explicit_qr")) {
    if (identical(.amatrix_qr_helper_path(qr), "compact_factor") && .amatrix_qr_compact_available(qr)) {
      return(qr.resid(.amatrix_qr_factor(qr), y))
    }
    return(y - .amatrix_qr_fitted_value(qr, y, k = .amatrix_qr_rank(qr)))
  }

  qr.resid(.amatrix_qr_factor(qr), y)
}

.amatrix_rewrap_qr_result <- function(qr, value) {
  if ((inherits(value, "Matrix") || is.matrix(value)) && is.numeric(value)) {
    return(adgeMatrix(value, preferred_backend = qr$backend, precision = qr$precision))
  }
  value
}

#' @export
print.amDenseQR <- function(x, ...) {
  qr_obj <- .amatrix_unwrap_qr(x)
  cat(sprintf(
    "amQR [backend=%s|precision=%s|rank=%d|repr=%s|dim=%s|thin=%s|pivoted=%s|source=%s]\n",
    x$backend,
    x$precision,
    .amatrix_qr_rank(x),
    .amatrix_qr_kind(x),
    paste(.amatrix_qr_source_dim(x), collapse = "x"),
    .amatrix_qr_thin(x),
    .amatrix_qr_pivoted(x),
    x$source_class
  ))
  if (.amatrix_qr_compact_available(x)) {
    factor_tag <- .amatrix_qr_or(x$state$factor_source, "available")
    if (!.amatrix_qr_compact_materialized(x)) {
      factor_tag <- sprintf("%s [lazy]", factor_tag)
    }
    cat(sprintf("  compact_factor: %s\n", factor_tag))
  }
  if (!is.null(.amatrix_qr_q_key(x))) {
    cat(sprintf("  resident_q: %s\n", .amatrix_qr_q_key(x)))
  }
  print(qr_obj, ...)
  invisible(x)
}

#' @export
dim.amDenseQR <- function(x) {
  .amatrix_qr_source_dim(x)
}

#' @noRd
setMethod("qr.Q", signature(qr = "amQR"), function(qr, complete = FALSE) {
  .amatrix_rewrap_qr_result(qr, .amatrix_qr_q(qr, complete = complete))
})

#' @noRd
setMethod("qr.R", signature(qr = "amQR"), function(qr, complete = FALSE) {
  .amatrix_rewrap_qr_result(qr, .amatrix_qr_r(qr, complete = complete))
})

#' @noRd
setMethod("qr.solve", signature(a = "amQR", b = "missing"), function(a, b, tol = 1e-07) {
  .amatrix_rewrap_qr_result(a, .amatrix_qr_solve_value(a, tol = tol))
})

#' @noRd
setMethod("qr.solve", signature(a = "amQR", b = "ANY"), function(a, b, tol = 1e-07) {
  b_mat <- as.matrix(.amatrix_host_arg(b))
  .amatrix_rewrap_qr_result(a, .amatrix_qr_solve_value(a, b = b_mat, tol = tol))
})

#' @noRd
setMethod("qr.coef", signature(qr = "amQR", y = "ANY"), function(qr, y) {
  y_mat <- as.matrix(.amatrix_host_arg(y))
  .amatrix_rewrap_qr_result(qr, .amatrix_qr_coef_value(qr, y_mat))
})

#' @noRd
setMethod("qr.qty", signature(qr = "amQR", y = "ANY"), function(qr, y) {
  y_mat <- as.matrix(.amatrix_host_arg(y))
  if (identical(.amatrix_qr_kind(qr), "mlx_compact_qr")) {
    return(.amatrix_rewrap_qr_result(qr, .amatrix_mlx_compact_qr_qty(qr, y_mat)))
  }
  if (identical(.amatrix_qr_kind(qr), "explicit_qr")) {
    if (identical(.amatrix_qr_helper_path(qr), "compact_factor") && .amatrix_qr_compact_available(qr)) {
      return(.amatrix_rewrap_qr_result(qr, qr.qty(.amatrix_qr_factor(qr), y_mat)))
    }
    return(.amatrix_rewrap_qr_result(qr, .amatrix_explicit_qr_qty(qr, y_mat)))
  }
  .amatrix_rewrap_qr_result(qr, qr.qty(.amatrix_qr_factor(qr), y_mat))
})

#' @noRd
setMethod("qr.qy", signature(qr = "amQR", y = "ANY"), function(qr, y) {
  y_mat <- as.matrix(.amatrix_host_arg(y))
  if (identical(.amatrix_qr_kind(qr), "mlx_compact_qr")) {
    return(.amatrix_rewrap_qr_result(qr, .amatrix_mlx_compact_qr_qy(qr, y_mat)))
  }
  if (identical(.amatrix_qr_kind(qr), "explicit_qr")) {
    if (identical(.amatrix_qr_helper_path(qr), "compact_factor") && .amatrix_qr_compact_available(qr)) {
      return(.amatrix_rewrap_qr_result(qr, qr.qy(.amatrix_qr_factor(qr), y_mat)))
    }
    return(.amatrix_rewrap_qr_result(qr, .amatrix_explicit_qr_qy(qr, y_mat)))
  }
  .amatrix_rewrap_qr_result(qr, qr.qy(.amatrix_qr_factor(qr), y_mat))
})

#' @noRd
setMethod("qr.fitted", signature(qr = "amQR", y = "ANY"), function(qr, y, k = NULL) {
  y_mat <- as.matrix(.amatrix_host_arg(y))
  .amatrix_rewrap_qr_result(qr, .amatrix_qr_fitted_value(qr, y_mat, k = k))
})

#' @noRd
setMethod("qr.resid", signature(qr = "amQR", y = "ANY"), function(qr, y) {
  y_mat <- as.matrix(.amatrix_host_arg(y))
  .amatrix_rewrap_qr_result(qr, .amatrix_qr_resid_value(qr, y_mat))
})
