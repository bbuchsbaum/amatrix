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
#' @param x An \code{adgeMatrix} or plain matrix.
#' @param backend Backend name (default: inferred from \code{x}).
#' @return A \code{resident_handle} environment.
#' @export
resident_handle <- function(x, backend = NULL) {
  if (inherits(x, "adgeMatrix")) {
    bk_name <- if (!is.null(backend)) backend else x@preferred_backend
  } else if (is.matrix(x)) {
    bk_name <- if (!is.null(backend)) backend else "cpu"
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
  h$active <- TRUE
  class(h) <- "resident_handle"

  reg.finalizer(h, function(env) {
    if (isTRUE(env$active) && !is.null(env$resident_key)) {
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

# ── In-place operations ──────────────────────────────────────────────────────

#' In-place broadcast sweep on a resident handle
#'
#' Replaces the handle's resident key with the result of the sweep.
#' Equivalent to \code{am_sweep(x, MARGIN, STATS, FUN)} but mutates in place.
#'
#' @param h A \code{resident_handle}.
#' @param MARGIN 1L (rows) or 2L (columns).
#' @param STATS Numeric vector of statistics.
#' @param FUN Operation: \code{"+"}, \code{"-"}, \code{"*"}, or \code{"/"}.
#' @return \code{h}, invisibly (modified in place).
#' @export
am_sweep_inplace <- function(h, MARGIN, STATS, FUN = "+") {
  .rh_check(h)
  backend <- .rh_backend(h)
  if (!is.function(backend$broadcast_ewise_resident))
    stop("backend does not support broadcast_ewise_resident")

  old_key <- h$resident_key
  new_key <- .amatrix_next_resident_key(h$backend_name)
  backend$broadcast_ewise_resident(old_key, as.double(STATS),
                                    as.integer(MARGIN), as.character(FUN),
                                    new_key, defer = TRUE)
  # Drop the old key and update handle
  if (isTRUE(backend$resident_has(old_key)))
    backend$resident_drop(old_key)
  h$resident_key <- new_key
  invisible(h)
}

#' In-place elementwise operation on a resident handle
#'
#' @param h A \code{resident_handle}.
#' @param rhs A scalar or another \code{resident_handle}.
#' @param op Operation string: \code{"+"}, \code{"-"}, \code{"*"}, \code{"/"}.
#' @return \code{h}, invisibly (modified in place).
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
  backend$ewise_resident(old_key, rhs_payload, op, new_key, defer = TRUE)
  if (isTRUE(backend$resident_has(old_key)))
    backend$resident_drop(old_key)
  h$resident_key <- new_key
  invisible(h)
}

# ── Reductions (return plain vectors) ────────────────────────────────────────
# Note: rowSums/colSums are S4 generics in amatrix, so S3 dispatch does not
# work for non-S4 classes.  Use rh_rowSums/rh_colSums for resident handles.

#' Row sums of a resident handle (GPU-resident reduction)
#' @param h A \code{resident_handle}.
#' @return Numeric vector of row sums.
#' @export
rh_rowSums <- function(h) {
  .rh_check(h)
  backend <- .rh_backend(h)
  if (is.function(backend$rowSums_resident))
    return(backend$rowSums_resident(h$resident_key))
  base::rowSums(as.matrix(h))
}

#' Column sums of a resident handle (GPU-resident reduction)
#' @param h A \code{resident_handle}.
#' @return Numeric vector of column sums.
#' @export
rh_colSums <- function(h) {
  .rh_check(h)
  backend <- .rh_backend(h)
  if (is.function(backend$colSums_resident))
    return(backend$colSums_resident(h$resident_key))
  base::colSums(as.matrix(h))
}

# ── Materialization ──────────────────────────────────────────────────────────

#' @export
as.matrix.resident_handle <- function(x, ...) {
  .rh_check(x)
  backend <- .rh_backend(x)
  backend$resident_materialize(x$resident_key)
}

#' @export
dim.resident_handle <- function(x) x$dim

#' @export
nrow.resident_handle <- function(x) x$dim[1L]

#' @export
ncol.resident_handle <- function(x) x$dim[2L]

#' Convert a resident handle back to an adgeMatrix
#'
#' Materialises the GPU data and creates an adgeMatrix with the resident key
#' still bound.  The handle becomes inert after this call.
#'
#' @param h A \code{resident_handle}.
#' @return An \code{adgeMatrix}.
as_adgeMatrix.resident_handle <- function(h, ...) {
  .rh_check(h)
  mat <- as.matrix(h)
  obj <- new_adgeMatrix(mat, preferred_backend = h$backend_name)
  .amatrix_bind_resident(obj, h$backend_name, h$resident_key)
  # Transfer ownership: handle no longer drops the key on GC
  h$active <- FALSE
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
