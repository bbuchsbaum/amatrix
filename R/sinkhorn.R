#' Doubly-stochastic scaling via Sinkhorn-Knopp iterations
#'
#' Alternates row and column normalization until the matrix is approximately
#' doubly stochastic. When the chosen backend supports resident broadcast and
#' reduction kernels, the hot loop stays on device via \code{resident_handle}
#' and returns a deferred \code{adgeMatrix} bound to the resident result.
#'
#' @param A A dense numeric matrix or \code{adgeMatrix}. Sparse inputs are not
#'   yet supported in this surface.
#' @param max_iter Maximum number of Sinkhorn iterations.
#' @param tol Convergence tolerance on the maximum row/column sum error.
#' @param check_every Check convergence every \code{check_every} iterations.
#' @param eps Floor applied to row/column sums before division.
#' @param mode Execution mode used when coercing a plain matrix. Default
#'   \code{"fast"} allows accelerated backends to use lower precision.
#' @param backend Backend name used when coercing a plain matrix. Ignored when
#'   \code{A} is already an \code{adgeMatrix}.
#' @param return_info When \code{TRUE}, return convergence metadata alongside
#'   the scaled matrix.
#'
#' @return By default, an \code{adgeMatrix}. With \code{return_info = TRUE}, a
#'   list containing \code{result}, \code{iterations}, \code{converged},
#'   \code{row_error}, \code{col_error}, \code{backend}, and \code{method}.
#'
#' @examples
#' A <- abs(matrix(rnorm(16), nrow = 4)) + 0.1
#' S <- sinkhorn(A, max_iter = 50L)
#' # Row sums should be close to 1
#' rowSums(as.matrix(S))
#'
#' @seealso \code{\link{dist_matrix}}
#' @export
sinkhorn <- function(
  A,
  max_iter = 200L,
  tol = 1e-8,
  check_every = 5L,
  eps = 1e-15,
  mode = "fast",
  backend = NULL,
  return_info = FALSE
) {
  max_iter <- .amatrix_sinkhorn_validate_count(max_iter, "max_iter")
  check_every <- .amatrix_sinkhorn_validate_count(check_every, "check_every")
  tol <- .amatrix_sinkhorn_validate_positive(tol, "tol", allow_zero = TRUE)
  eps <- .amatrix_sinkhorn_validate_positive(eps, "eps", allow_zero = FALSE)

  if (!is.logical(return_info) || length(return_info) != 1L || is.na(return_info)) {
    stop("return_info must be TRUE or FALSE", call. = FALSE)
  }

  prepared <- .amatrix_sinkhorn_prepare_input(A, mode = mode, backend = backend)

  resident_backend <- .amatrix_sinkhorn_resident_backend(prepared$arg)
  fit <- if (!is.null(resident_backend)) {
    .amatrix_sinkhorn_resident(
      prepared$arg,
      max_iter = max_iter,
      tol = tol,
      check_every = check_every,
      eps = eps
    )
  } else {
    .amatrix_sinkhorn_host(
      prepared$host,
      meta = prepared$meta,
      max_iter = max_iter,
      tol = tol,
      check_every = check_every,
      eps = eps
    )
  }

  if (isTRUE(return_info)) {
    fit
  } else {
    fit$result
  }
}

.amatrix_sinkhorn_validate_count <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("%s must be a single positive integer", name), call. = FALSE)
  }
  value <- as.integer(value)
  if (value <= 0L) {
    stop(sprintf("%s must be a single positive integer", name), call. = FALSE)
  }
  value
}

.amatrix_sinkhorn_validate_positive <- function(value, name, allow_zero = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || !is.finite(value)) {
    stop(sprintf("%s must be a single finite numeric value", name), call. = FALSE)
  }
  if (allow_zero) {
    if (value < 0) {
      stop(sprintf("%s must be non-negative", name), call. = FALSE)
    }
  } else if (value <= 0) {
    stop(sprintf("%s must be positive", name), call. = FALSE)
  }
  as.double(value)
}

.amatrix_sinkhorn_prepare_input <- function(A, mode, backend) {
  if (inherits(A, "adgCMatrix") || inherits(A, "sparseMatrix")) {
    stop("sinkhorn() currently requires a dense matrix or adgeMatrix", call. = FALSE)
  }

  if (inherits(A, "adgeMatrix")) {
    host <- as.matrix(amatrix_materialize_host(A))
    .amatrix_sinkhorn_validate_host(host)
    return(list(
      host = host,
      arg = A,
      meta = list(
        preferred_backend = A@preferred_backend,
        policy = A@policy,
        precision = A@precision
      )
    ))
  }

  if (!(is.matrix(A) || inherits(A, "denseMatrix"))) {
    stop("A must be a dense numeric matrix or adgeMatrix", call. = FALSE)
  }

  host <- as.matrix(A)
  .amatrix_sinkhorn_validate_host(host)
  params <- .amatrix_resolve_mode(
    mode = mode,
    backend = backend,
    preferred_backend = NULL,
    policy = NULL,
    precision = NULL
  )

  arg <- NULL
  if (!is.null(backend) && !identical(params$preferred_backend, "cpu")) {
    arg <- new_adgeMatrix(
      host,
      preferred_backend = params$preferred_backend,
      policy = params$policy,
      precision = params$precision
    )
  }

  list(
    host = host,
    arg = arg,
    meta = params
  )
}

.amatrix_sinkhorn_validate_host <- function(host) {
  if (!is.numeric(host)) {
    stop("A must be numeric", call. = FALSE)
  }
  if (length(dim(host)) != 2L) {
    stop("A must be two-dimensional", call. = FALSE)
  }
  if (!identical(nrow(host), ncol(host))) {
    stop("sinkhorn() currently requires a square matrix", call. = FALSE)
  }
  if (anyNA(host) || any(!is.finite(host))) {
    stop("A must contain only finite non-missing values", call. = FALSE)
  }
  if (any(host < 0)) {
    stop("A must be elementwise non-negative", call. = FALSE)
  }
  if (any(base::rowSums(host) <= 0) || any(base::colSums(host) <= 0)) {
    stop("A must have strictly positive row sums and column sums", call. = FALSE)
  }
  invisible(host)
}

.amatrix_sinkhorn_resident_backend <- function(x) {
  if (!inherits(x, "adgeMatrix")) {
    return(NULL)
  }

  backend_name <- .amatrix_live_resident_backend(x)
  if (is.null(backend_name)) {
    backend_name <- x@preferred_backend
  }
  if (!is.character(backend_name) || length(backend_name) != 1L || !nzchar(backend_name)) {
    return(NULL)
  }

  backend <- tryCatch(.amatrix_get_backend(backend_name), error = function(e) NULL)
  if (is.null(backend) || !.amatrix_backend_available_safe(backend) || !.amatrix_backend_residency_capable(backend)) {
    return(NULL)
  }

  required_ops <- c("broadcast_ewise", "rowSums", "colSums")
  supports <- vapply(
    required_ops,
    function(op) .amatrix_backend_supports_resident_op(backend, op),
    logical(1)
  )
  if (!all(supports)) {
    return(NULL)
  }

  backend_name
}

.amatrix_sinkhorn_finish <- function(result, iterations, converged, row_error, col_error, backend, method) {
  list(
    result = result,
    iterations = as.integer(iterations),
    converged = isTRUE(converged),
    row_error = as.double(row_error),
    col_error = as.double(col_error),
    backend = as.character(backend),
    method = as.character(method)
  )
}

.amatrix_sinkhorn_host <- function(host, meta, max_iter, tol, check_every, eps) {
  work <- host
  converged <- FALSE
  row_error <- Inf
  col_error <- Inf

  for (iter in seq_len(max_iter)) {
    row_scale <- pmax(base::rowSums(work), eps)
    work <- work / row_scale

    col_scale <- pmax(base::colSums(work), eps)
    work <- t(t(work) / col_scale)

    if ((iter %% check_every) == 0L || iter == max_iter) {
      row_error <- max(abs(base::rowSums(work) - 1))
      col_error <- max(abs(base::colSums(work) - 1))
      if (max(row_error, col_error) < tol) {
        converged <- TRUE
        break
      }
    }
  }

  result <- new_adgeMatrix(
    work,
    preferred_backend = meta$preferred_backend,
    policy = meta$policy,
    precision = meta$precision
  )

  .amatrix_sinkhorn_finish(
    result = result,
    iterations = iter,
    converged = converged,
    row_error = row_error,
    col_error = col_error,
    backend = "cpu",
    method = "host"
  )
}

.amatrix_sinkhorn_resident <- function(A, max_iter, tol, check_every, eps) {
  backend_name <- .amatrix_sinkhorn_resident_backend(A)
  if (is.null(backend_name)) {
    stop("resident sinkhorn path requires resident row/col reductions and broadcast kernels", call. = FALSE)
  }

  h <- resident_handle(A, backend = backend_name)
  # Register cleanup immediately after acquiring the handle. Only drop a key the
  # handle actually owns: when resident_handle reused an already-resident
  # matrix's buffer (owns_key = FALSE) that buffer belongs to the caller's
  # object and must not be freed here.
  on.exit({
    if (exists("h", inherits = FALSE) && isTRUE(h$active) &&
        isTRUE(h$owns_key) && !is.null(h$resident_key)) {
      backend <- tryCatch(.amatrix_get_backend(h$backend_name), error = function(e) NULL)
      if (!is.null(backend) && is.function(backend$resident_drop) &&
          isTRUE(backend$resident_has(h$resident_key))) {
        try(backend$resident_drop(h$resident_key), silent = TRUE)
      }
      h$active <- FALSE
      h$resident_key <- NULL
    }
  }, add = TRUE)
  backend <- .rh_backend(h)
  has_vector_chain <- is.function(backend$rowSums_resident_key) &&
    is.function(backend$colSums_resident_key) &&
    is.function(backend$broadcast_ewise_resident_key)

  converged <- FALSE
  row_error <- Inf
  col_error <- Inf

  for (iter in seq_len(max_iter)) {
    row_done <- FALSE
    col_done <- FALSE

    if (isTRUE(has_vector_chain)) {
      row_key <- .rh_axis_sums_key(h, 1L)
      if (is.null(row_key)) {
        has_vector_chain <- FALSE
      } else {
        .rh_sweep_inplace_key(h, 1L, row_key, "/", drop_stats = TRUE)
        row_done <- TRUE

        col_key <- .rh_axis_sums_key(h, 2L)
        if (is.null(col_key)) {
          has_vector_chain <- FALSE
        } else {
          .rh_sweep_inplace_key(h, 2L, col_key, "/", drop_stats = TRUE)
          col_done <- TRUE
        }
      }
    }

    if (!isTRUE(row_done)) {
      row_scale <- pmax(rh_rowSums(h), eps)
      am_sweep_inplace(h, 1L, row_scale, "/")
    }
    if (!isTRUE(col_done)) {
      col_scale <- pmax(rh_colSums(h), eps)
      am_sweep_inplace(h, 2L, col_scale, "/")
    }

    if (tol > 0 && ((iter %% check_every) == 0L || iter == max_iter)) {
      row_error <- max(abs(rh_rowSums(h) - 1))
      col_error <- max(abs(rh_colSums(h) - 1))
      if (max(row_error, col_error) < tol) {
        converged <- TRUE
        break
      }
    }
  }

  result <- as_adgeMatrix.resident_handle(h, defer_host = TRUE)

  .amatrix_sinkhorn_finish(
    result = result,
    iterations = iter,
    converged = converged,
    row_error = row_error,
    col_error = col_error,
    backend = backend_name,
    method = "resident"
  )
}
