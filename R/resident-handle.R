# ── Mutable Resident Handle ───────────────────────────────────────────────────
#
# A lightweight mutable wrapper around a GPU-resident array.  Unlike adgeMatrix
# (immutable S4, ~5ms per object), a handle is a plain environment (~0.01ms to
# update).  Operations mutate the handle's resident key in place, so iterative
# algorithms avoid the per-step S4 allocation overhead.
#
# Usage:
#   h <- resident_handle(A_gpu)
#   for (i in 1:100) {
#     rs <- rowSums(h)
#     am_sweep_inplace(h, 1L, 1/rs, "*")
#   }
#   result <- as.matrix(h)   # single download at the end
#
# The handle owns its resident key.  When converted back to adgeMatrix via
# as_adgeMatrix(), ownership transfers and the handle becomes inert.

#' Create a mutable GPU-resident handle
#'
#' Wraps an \code{adgeMatrix} or plain matrix in a lightweight mutable
#' environment that holds a GPU-resident buffer key. Unlike
#' \code{adgeMatrix}, the handle can be updated in place, making it
#' suitable for iterative algorithms that would otherwise incur per-step
#' S4 object allocation overhead. The handle owns its resident key and
#' releases the device buffer when garbage collected.
#'
#' @param x An \code{adgeMatrix} or plain \code{matrix}. If \code{x}
#'   is already GPU-resident on \code{backend}, the existing device
#'   buffer is reused without re-uploading.
#' @param backend Character string. Name of the backend to use.
#'   Defaults to \code{x@@preferred_backend} for \code{adgeMatrix}
#'   inputs and \code{"cpu"} for plain matrices. The backend must
#'   support GPU residency.
#'
#' @return A \code{resident_handle} environment with fields
#'   \code{backend_name}, \code{resident_key}, \code{dim},
#'   \code{dimnames}, \code{policy}, \code{precision}, and
#'   \code{active}.
#'
#' @examples
#' \donttest{
#' m <- adgeMatrix(matrix(runif(12), 3, 4), preferred_backend = "cpu")
#' # resident_handle requires a backend with residency support (e.g. MLX, OpenCL)
#' }
#'
#' @seealso \code{\link{am_sweep_inplace}}, \code{\link{rh_rowSums}},
#'   \code{\link{rh_colSums}}
#' @export
resident_handle <- function(x, backend = NULL) {
  reused_existing_key <- FALSE
  if (inherits(x, "adgeMatrix")) {
    bk_name <- if (!is.null(backend)) backend else x@preferred_backend
    policy <- x@policy
    precision <- x@precision
    handle_dimnames <- x@Dimnames
  } else if (is.matrix(x)) {
    bk_name <- if (!is.null(backend)) backend else "cpu"
    policy <- amatrix_default_policy()
    precision <- amatrix_default_precision()
    handle_dimnames <- base::dimnames(x)
  } else {
    stop("x must be an adgeMatrix or matrix")
  }

  backend_obj <- .amatrix_get_backend(bk_name)
  if (!.amatrix_backend_residency_capable(backend_obj))
    stop("backend '", bk_name, "' does not support residency")

  # Try to reuse existing resident key; otherwise upload
  key <- NULL
  if (inherits(x, "adgeMatrix")) {
    entry <- .amatrix_resident_entry(x)
    if (!is.null(entry) && identical(entry$backend, bk_name) &&
        isTRUE(backend_obj$resident_has(entry$resident_key))) {
      key <- entry$resident_key
      reused_existing_key <- TRUE
    }
  }
  if (is.null(key)) {
    key <- .amatrix_next_resident_key(bk_name)
    mat <- if (is.matrix(x)) x else as.matrix(x)
    if (!is.double(mat)) storage.mode(mat) <- "double"
    backend_obj$resident_store(key, mat)
  }

  h <- new.env(parent = emptyenv())
  h$backend_name <- bk_name
  h$resident_key <- key
  h$dim <- if (inherits(x, "adgeMatrix")) x@Dim else as.integer(dim(x))
  h$dimnames <- if (is.null(handle_dimnames)) vector("list", 2L) else handle_dimnames
  h$policy <- policy
  h$precision <- precision
  h$active <- TRUE
  h$owns_key <- !reused_existing_key
  class(h) <- "resident_handle"

  reg.finalizer(h, function(env) {
    if (isTRUE(env$active) && isTRUE(env$owns_key) && !is.null(env$resident_key)) {
      bk <- tryCatch(.amatrix_get_backend(env$backend_name), error = function(e) NULL)
      if (!is.null(bk) && is.function(bk$resident_drop) &&
          isTRUE(bk$resident_has(env$resident_key))) {
        try(bk$resident_drop(env$resident_key), silent = TRUE)
      }
    }
  }, onexit = TRUE)

  h
}

.rh_check <- function(h) {
  if (!isTRUE(h$active)) stop("resident_handle is no longer active")
  invisible(h)
}

.rh_backend <- function(h) {
  .amatrix_get_backend(h$backend_name)
}

.rh_drop_key <- function(backend, key) {
  if (!is.null(key) && is.function(backend$resident_has) &&
      isTRUE(backend$resident_has(key)) && is.function(backend$resident_drop)) {
    backend$resident_drop(key)
  }
  invisible(NULL)
}

.rh_axis_sums_key <- function(h, margin, na.rm = FALSE, dims = 1L) {
  .rh_check(h)
  backend <- .rh_backend(h)
  fn <- if (identical(as.integer(margin), 1L)) backend$rowSums_resident_key else backend$colSums_resident_key
  if (!is.function(fn)) {
    return(NULL)
  }

  out_key <- .amatrix_next_resident_key(h$backend_name)
  ok <- FALSE
  tryCatch({
    fn(h$resident_key, out_key, na.rm = na.rm, dims = dims)
    ok <- TRUE
  }, error = function(e) {
    .rh_drop_key(backend, out_key)
  })

  if (!ok) {
    return(NULL)
  }

  out_key
}

.rh_sweep_inplace_key <- function(h, MARGIN, stats_key, FUN = "+", drop_stats = FALSE) {
  .rh_check(h)
  backend <- .rh_backend(h)
  if (is.function(backend$broadcast_ewise_resident_inplace_key)) {
    ok <- FALSE
    tryCatch({
      backend$broadcast_ewise_resident_inplace_key(
        h$resident_key,
        as.character(stats_key),
        as.integer(MARGIN),
        as.character(FUN)
      )
      ok <- TRUE
    }, error = function(e) NULL)

    if (isTRUE(drop_stats)) {
      .rh_drop_key(backend, stats_key)
    }
    if (ok) {
      return(invisible(h))
    }
    stop("backend failed to apply resident vector sweep", call. = FALSE)
  }

  if (!is.function(backend$broadcast_ewise_resident_key)) {
    stop("backend does not support broadcast_ewise_resident_key")
  }

  old_key <- h$resident_key
  new_key <- .amatrix_next_resident_key(h$backend_name)
  ok <- FALSE

  tryCatch({
    backend$broadcast_ewise_resident_key(
      old_key,
      as.character(stats_key),
      as.integer(MARGIN),
      as.character(FUN),
      new_key,
      defer = TRUE
    )
    ok <- TRUE
  }, error = function(e) {
    .rh_drop_key(backend, new_key)
  })

  if (!ok) {
    if (isTRUE(drop_stats)) {
      .rh_drop_key(backend, stats_key)
    }
    stop("backend failed to apply resident vector sweep", call. = FALSE)
  }

  .rh_drop_key(backend, old_key)
  if (isTRUE(drop_stats)) {
    .rh_drop_key(backend, stats_key)
  }
  h$resident_key <- new_key
  invisible(h)
}

# ── In-place operations ──────────────────────────────────────────────────────

#' In-place broadcast sweep on a resident handle
#'
#' Applies a row-wise or column-wise arithmetic operation between the
#' resident matrix and a statistics vector, mutating the handle in
#' place. Equivalent to \code{sweep(as.matrix(h), MARGIN, STATS, FUN)}
#' but avoids downloading the matrix to host.
#'
#' @param h A \code{resident_handle}.
#' @param MARGIN Integer. \code{1L} to sweep across rows (one value per
#'   row), \code{2L} to sweep across columns.
#' @param STATS Numeric vector of length equal to the number of rows or
#'   columns selected by \code{MARGIN}.
#' @param FUN Character string. Arithmetic operator to apply:
#'   \code{"+"}, \code{"-"}, \code{"*"}, or \code{"/"}.
#'   Default \code{"+"}.
#'
#' @return \code{h}, invisibly. The handle is modified in place; the
#'   underlying device buffer is replaced with the sweep result.
#'
#' @examples
#' \donttest{
#' # requires a backend with residency support (e.g. MLX, OpenCL)
#' }
#'
#' @seealso \code{\link{resident_handle}}, \code{\link{am_ewise_inplace}}
#' @export
am_sweep_inplace <- function(h, MARGIN, STATS, FUN = "+") {
  .rh_check(h)
  backend <- .rh_backend(h)
  if (is.function(backend$broadcast_ewise_resident_inplace)) {
    backend$broadcast_ewise_resident_inplace(
      h$resident_key,
      as.double(STATS),
      as.integer(MARGIN),
      as.character(FUN)
    )
    h$owns_key <- identical(.amatrix_update_resident_aliases(h$backend_name, h$resident_key), 0L)
    return(invisible(h))
  }
  if (!is.function(backend$broadcast_ewise_resident))
    stop("backend does not support broadcast_ewise_resident")

  old_key <- h$resident_key
  new_key <- .amatrix_next_resident_key(h$backend_name)
  err <- NULL
  tryCatch(
    backend$broadcast_ewise_resident(
      old_key,
      as.double(STATS),
      as.integer(MARGIN),
      as.character(FUN),
      new_key,
      defer = TRUE
    ),
    error = function(e) {
      err <<- e
      .rh_drop_key(backend, new_key)
    }
  )
  if (!is.null(err)) {
    stop(err)
  }

  alias_count <- .amatrix_update_resident_aliases(h$backend_name, old_key, new_key)
  # Drop the old key and update handle
  if (isTRUE(backend$resident_has(old_key)))
    backend$resident_drop(old_key)
  h$resident_key <- new_key
  h$owns_key <- identical(alias_count, 0L)
  invisible(h)
}

#' In-place elementwise operation on a resident handle
#'
#' Applies an elementwise arithmetic operation between the handle's
#' resident matrix and either a scalar or another resident handle,
#' replacing the handle's device buffer with the result.
#'
#' @param h A \code{resident_handle}.
#' @param rhs A length-1 numeric scalar, or a \code{resident_handle}
#'   with identical dimensions to \code{h}.
#' @param op Character string. Arithmetic operator: \code{"+"}, \code{"-"},
#'   \code{"*"}, or \code{"/"}.
#'
#' @return \code{h}, invisibly. The handle is modified in place.
#'
#' @examples
#' \donttest{
#' # requires a backend with residency support (e.g. MLX, OpenCL)
#' }
#'
#' @seealso \code{\link{am_sweep_inplace}}, \code{\link{resident_handle}}
#' @export
am_ewise_inplace <- function(h, rhs, op) {
  .rh_check(h)
  backend <- .rh_backend(h)
  if (!is.function(backend$ewise_resident))
    stop("backend does not support ewise_resident")

  rhs_payload <- if (inherits(rhs, "resident_handle")) rhs$resident_key
                 else if (is.numeric(rhs) && length(rhs) == 1L) as.double(rhs)
                 else stop("rhs must be a scalar or resident_handle")

  old_key <- h$resident_key
  new_key <- .amatrix_next_resident_key(h$backend_name)
  err <- NULL
  tryCatch(
    backend$ewise_resident(old_key, rhs_payload, op, new_key, defer = TRUE),
    error = function(e) {
      err <<- e
      .rh_drop_key(backend, new_key)
    }
  )
  if (!is.null(err)) {
    stop(err)
  }

  alias_count <- .amatrix_update_resident_aliases(h$backend_name, old_key, new_key)
  if (isTRUE(backend$resident_has(old_key)))
    backend$resident_drop(old_key)
  h$resident_key <- new_key
  h$owns_key <- identical(alias_count, 0L)
  invisible(h)
}

# ── Reductions (return plain vectors) ────────────────────────────────────────
# Note: rowSums/colSums are S4 generics in amatrix, so S3 dispatch does not
# work for non-S4 classes.  Use rh_rowSums/rh_colSums for resident handles.

#' Row sums of a GPU-resident handle
#'
#' Computes row sums of the matrix stored in a \code{resident_handle},
#' using a GPU-resident reduction when the backend supports it to avoid
#' a round-trip download. Falls back to \code{base::rowSums} on the
#' materialized matrix when no resident reduction is available.
#'
#' @param h A \code{resident_handle}.
#'
#' @return Numeric vector of length \code{nrow(h)}.
#'
#' @examples
#' \donttest{
#' # requires a backend with residency support (e.g. MLX, OpenCL)
#' }
#'
#' @seealso \code{\link{rh_colSums}}, \code{\link{am_sweep_inplace}}
#' @export
rh_rowSums <- function(h) {
  .rh_check(h)
  backend <- .rh_backend(h)
  if (is.function(backend$rowSums_resident))
    return(backend$rowSums_resident(h$resident_key))
  base::rowSums(as.matrix(h))
}

#' Column sums of a GPU-resident handle
#'
#' Computes column sums of the matrix stored in a
#' \code{resident_handle}, using a GPU-resident reduction when the
#' backend supports it. Falls back to \code{base::colSums} on the
#' materialized matrix when no resident reduction is available.
#'
#' @param h A \code{resident_handle}.
#'
#' @return Numeric vector of length \code{ncol(h)}.
#'
#' @examples
#' \donttest{
#' # requires a backend with residency support (e.g. MLX, OpenCL)
#' }
#'
#' @seealso \code{\link{rh_rowSums}}, \code{\link{am_sweep_inplace}}
#' @export
rh_colSums <- function(h) {
  .rh_check(h)
  backend <- .rh_backend(h)
  if (is.function(backend$colSums_resident))
    return(backend$colSums_resident(h$resident_key))
  base::colSums(as.matrix(h))
}

# ── Materialization ──────────────────────────────────────────────────────────

#' @rdname amatrix-s3-methods
#' @export
as.matrix.resident_handle <- function(x, ...) {
  .rh_check(x)
  backend <- .rh_backend(x)
  mat <- backend$resident_materialize(x$resident_key)
  if (is.matrix(mat) && any(!vapply(x$dimnames, is.null, logical(1)))) {
    base::dimnames(mat) <- x$dimnames
  }
  mat
}

#' @rdname amatrix-s3-methods
#' @export
dim.resident_handle <- function(x) x$dim

#' @rdname amatrix-s3-methods
#' @export
nrow.resident_handle <- function(x) x$dim[1L]

#' @rdname amatrix-s3-methods
#' @export
ncol.resident_handle <- function(x) x$dim[2L]

#' Convert a resident handle back to an adgeMatrix
#'
#' Materialises the GPU data and creates an adgeMatrix with the resident key
#' still bound.  The handle becomes inert after this call.
#'
#' @param h A \code{resident_handle}.
#' @param ... Reserved for future use.
#' @param defer_host When \code{TRUE}, return a deferred-host
#'   \code{adgeMatrix} that materializes lazily.
#' @return An \code{adgeMatrix}.
as_adgeMatrix.resident_handle <- function(h, ..., defer_host = FALSE) {
  .rh_check(h)
  if (isTRUE(defer_host)) {
    obj <- new_adgeMatrix_deferred(
      dim = h$dim,
      dimnames = h$dimnames,
      preferred_backend = h$backend_name,
      policy = h$policy,
      precision = h$precision
    )
  } else {
    mat <- as.matrix(h)
    obj <- new_adgeMatrix(
      mat,
      preferred_backend = h$backend_name,
      policy = h$policy,
      precision = h$precision
    )
  }
  .amatrix_bind_resident(obj, h$backend_name, h$resident_key)
  # Transfer ownership: handle no longer drops the key on GC
  h$active <- FALSE
  h$owns_key <- FALSE
  h$resident_key <- NULL
  obj
}

#' @export
print.resident_handle <- function(x, ...) {
  state <- if (isTRUE(x$active)) "active" else "inert"
  cat(sprintf("resident_handle [%dx%d | %s | %s]\n",
              x$dim[1L], x$dim[2L], x$backend_name, state))
  invisible(x)
}
