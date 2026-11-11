# TRUE when `value` is a logical matrix. Deliberately inspects only the class:
# calling as.matrix() here would allocate a dense n x p copy of every sparse
# result just to look at its type.
.amatrix_value_is_logical <- function(value) {
  if (inherits(value, "Matrix")) {
    return(inherits(value, c("lMatrix", "nMatrix")))
  }
  is.matrix(value) && is.logical(value)
}

.amatrix_rewrap_like <- function(template, value) {
  if (inherits(template, "adlgeMatrix")) {
    if (inherits(value, "sparseMatrix") && inherits(value, "lsparseMatrix")) {
      return(new_adlgCMatrix(
        as(value, "lgCMatrix"),
        preferred_backend = template@preferred_backend,
        policy = template@policy,
        precision = template@precision
      ))
    }
    if (.amatrix_value_is_logical(value)) {
      return(new_adlgeMatrix(
        value,
        preferred_backend = template@preferred_backend,
        policy = template@policy,
        precision = template@precision
      ))
    }
    return(new_adgeMatrix(
      as.matrix(value),
      preferred_backend = template@preferred_backend,
      policy = template@policy,
      precision = template@precision
    ))
  }

  if (inherits(template, "adlgCMatrix")) {
    if (inherits(value, "sparseMatrix") && inherits(value, "lsparseMatrix")) {
      return(new_adlgCMatrix(
        as(value, "lgCMatrix"),
        preferred_backend = template@preferred_backend,
        policy = template@policy,
        precision = template@precision
      ))
    }
    if (.amatrix_value_is_logical(value)) {
      return(new_adlgeMatrix(
        value,
        preferred_backend = template@preferred_backend,
        policy = template@policy,
        precision = template@precision
      ))
    }
    if (inherits(value, "sparseMatrix")) {
      return(new_adgCMatrix(
        as(value, "dgCMatrix"),
        preferred_backend = template@preferred_backend,
        policy = template@policy,
        precision = template@precision
      ))
    }
    return(new_adgeMatrix(
      as.matrix(value),
      preferred_backend = template@preferred_backend,
      policy = template@policy,
      precision = template@precision
    ))
  }

  if (inherits(template, "adgeMatrix")) {
    # Matrix may hand back a sparse or diagonal result for a dense operand
    # (e.g. A * <mostly-zero matrix>). The template is dense, so densify.
    if (inherits(value, "Matrix") && !inherits(value, "denseMatrix")) {
      value <- as.matrix(value)
    }
    if (.amatrix_value_is_logical(value)) {
      return(new_adlgeMatrix(
        value,
        preferred_backend = template@preferred_backend,
        policy = template@policy,
        precision = template@precision
      ))
    }
    return(new_adgeMatrix(
      value,
      preferred_backend = template@preferred_backend,
      policy = template@policy,
      precision = template@precision
    ))
  }

  if (inherits(template, "adgCMatrix")) {
    if (inherits(value, "sparseMatrix") && inherits(value, "lsparseMatrix")) {
      return(new_adlgCMatrix(
        as(value, "lgCMatrix"),
        preferred_backend = template@preferred_backend,
        policy = template@policy,
        precision = template@precision
      ))
    }
    if (.amatrix_value_is_logical(value)) {
      return(new_adlgeMatrix(
        value,
        preferred_backend = template@preferred_backend,
        policy = template@policy,
        precision = template@precision
      ))
    }
    if (inherits(value, "sparseMatrix")) {
      return(new_adgCMatrix(
        as(value, "dgCMatrix"),
        preferred_backend = template@preferred_backend,
        policy = template@policy,
        precision = template@precision
      ))
    }
    return(new_adgeMatrix(
      as.matrix(value),
      preferred_backend = template@preferred_backend,
      policy = template@policy,
      precision = template@precision
    ))
  }

  value
}

# Promote a plain base matrix / Matrix operand to the matching amatrix class so
# the exported wrappers accept the inputs their documentation promises.
# amatrix objects and everything else (vectors, NULL, ...) pass through.
.amatrix_promote_arg <- function(x) {
  if (inherits(x, c("aMatrix", "aTransposeView"))) {
    return(x)
  }
  if (is.matrix(x) && (is.numeric(x) || is.logical(x))) {
    storage.mode(x) <- "double"
    return(as_adgeMatrix(x))
  }
  if (inherits(x, "sparseMatrix")) {
    return(as_adgCMatrix(methods::as(methods::as(methods::as(x, "dMatrix"), "generalMatrix"), "CsparseMatrix")))
  }
  if (inherits(x, "denseMatrix")) {
    return(as_adgeMatrix(methods::as(methods::as(methods::as(x, "dMatrix"), "generalMatrix"), "unpackedMatrix")))
  }
  x
}

.amatrix_host_arg <- function(value) {
  if (inherits(value, "aMatrix")) {
    return(amatrix_materialize_host(value))
  }
  value
}

.amatrix_subscript_arg <- function(value) {
  if (inherits(value, "aMatrix")) {
    return(as.matrix(amatrix_materialize_host(value)))
  }
  if (inherits(value, "Matrix")) {
    return(as.matrix(value))
  }
  value
}

.amatrix_is_logical_matrix_subscript <- function(value, target) {
  is.matrix(value) && is.logical(value) && identical(dim(value), dim(target))
}

.amatrix_logical_host_arg <- function(value) {
  if (inherits(value, "adlgCMatrix")) {
    return(as(amatrix_materialize_host(value), "lgCMatrix"))
  }

  if (inherits(value, "adlgeMatrix")) {
    return(as(amatrix_materialize_host(value), "lgeMatrix"))
  }

  if (inherits(value, "aMatrix")) {
    return(as.matrix(amatrix_materialize_host(value)))
  }

  value
}

.amatrix_is_numeric_matrix_value <- function(value) {
  if (inherits(value, "Matrix")) {
    if ("x" %in% slotNames(value)) {
      return(is.numeric(value@x))
    }
    return(FALSE)
  }

  if (is.matrix(value)) {
    return(is.numeric(value))
  }

  FALSE
}

.amatrix_is_logical_matrix_value <- function(value) {
  if (inherits(value, "Matrix")) {
    if ("x" %in% slotNames(value)) {
      return(is.logical(value@x))
    }
    return(FALSE)
  }

  if (is.matrix(value)) {
    return(is.logical(value))
  }

  FALSE
}

.amatrix_rewrap_value <- function(template, value) {
  if ((inherits(value, "Matrix") || is.matrix(value)) &&
      (.amatrix_is_numeric_matrix_value(value) || .amatrix_is_logical_matrix_value(value))) {
    return(.amatrix_rewrap_like(template, value))
  }
  value
}

.amatrix_template <- function(e1, e2 = NULL) {
  if (inherits(e1, "aMatrix")) {
    return(e1)
  }

  if (inherits(e2, "aMatrix")) {
    return(e2)
  }

  NULL
}

# Host operand for cbind2/rbind2. Logical amatrix classes must be stripped to
# their plain Matrix parents explicitly: amatrix_materialize_host() returns
# them unchanged, which made cbind2(aMatrix, aMatrix) re-dispatch to itself
# until the C stack overflowed.
.amatrix_bind_host_arg <- function(value) {
  if (inherits(value, c("adlgeMatrix", "adlgCMatrix"))) {
    return(.amatrix_logical_host_arg(value))
  }
  .amatrix_host_arg(value)
}

.amatrix_bind2 <- function(kind, x, y) {
  template <- .amatrix_template(x, y)
  value <- switch(
    kind,
    cbind2 = methods::cbind2(.amatrix_bind_host_arg(x), .amatrix_bind_host_arg(y)),
    rbind2 = methods::rbind2(.amatrix_bind_host_arg(x), .amatrix_bind_host_arg(y)),
    stop(sprintf("unsupported bind kind '%s'", kind), call. = FALSE)
  )

  if (is.null(template)) {
    return(value)
  }

  .amatrix_rewrap_value(template, value)
}

.amatrix_is_dense_matrix_like <- function(value) {
  inherits(value, "adgeMatrix") || inherits(value, "dgeMatrix") || is.matrix(value)
}

.amatrix_sparse_rhs_lowerable <- function(x, y) {
  inherits(y, "adgCMatrix") && (
    inherits(x, "adgeMatrix") ||
      inherits(x, "aTransposeView") ||
      inherits(x, "denseMatrix") ||
      is.matrix(x) ||
      (is.numeric(x) && is.null(dim(x)))
  )
}

.amatrix_retarget_for_sparse_rhs <- function(x, y, vector_as = c("row", "col")) {
  stopifnot(inherits(y, "adgCMatrix"))
  vector_as <- match.arg(vector_as)

  if (is.numeric(x) && is.null(dim(x))) {
    x <- if (identical(vector_as, "row")) {
      matrix(x, nrow = 1L)
    } else {
      matrix(x, ncol = 1L)
    }
  }

  if (inherits(x, "dgCMatrix")) {
    return(new_adgCMatrix(
      x,
      preferred_backend = y@preferred_backend,
      policy = y@policy,
      precision = y@precision
    ))
  }

  if (inherits(x, "aMatrix") || inherits(x, "denseMatrix")) {
    x <- as.matrix(amatrix_materialize_host(x))
  }

  new_adgeMatrix(
    as.matrix(x),
    preferred_backend = y@preferred_backend,
    policy = y@policy,
    precision = y@precision
  )
}

.amatrix_transpose_dimnames <- function(dn) {
  if (is.null(dn) || length(dn) != 2L) {
    return(list(NULL, NULL))
  }
  rev(dn)
}

.amatrix_transpose_dense_result <- function(value, template) {
  if (inherits(value, "adgeMatrix")) {
    backend_name <- .amatrix_live_resident_backend(value)
    if (!is.null(backend_name)) {
      backend <- .amatrix_get_backend(backend_name)
      if (is.function(backend$transpose_resident) &&
          .amatrix_backend_supports_resident_op(backend, "transpose", x = value)) {
        source <- .amatrix_prepare_resident_arg(value, backend_name, promote_amatrix = FALSE)
        if (!is.null(source)) {
          out_key <- .amatrix_next_resident_key(backend_name)
          ok <- tryCatch(
            {
              backend$transpose_resident(source$key, out_key)
              TRUE
            },
            error = function(e) {
              try(backend$resident_drop(out_key), silent = TRUE)
              FALSE
            }
          )
          if (ok) {
            obj <- new_adgeMatrix_deferred(
              dim = rev(dim(value)),
              dimnames = .amatrix_transpose_dimnames(value@Dimnames),
              preferred_backend = template@preferred_backend,
              policy = template@policy,
              precision = template@precision
            )
            return(.amatrix_bind_resident(obj, backend_name, out_key))
          }
        }
      }
    }
  }

  .amatrix_rewrap_like(template, t(as.matrix(.amatrix_host_arg(value))))
}

.amatrix_lower_sparse_rhs_matmul <- function(x, y) {
  lhs <- .amatrix_retarget_for_sparse_rhs(x, y, vector_as = "row")
  if (inherits(y, "adgCMatrix")) {
    choice <- .amatrix_backend_for(y, "matmul", y = lhs)
    resident <- .amatrix_try_dense_sparse_resident_matmul(lhs, y, choice$name)
    if (!is.null(resident)) {
      return(.amatrix_resident_wrap(lhs, resident, out_dim = c(nrow(lhs), ncol(y))))
    }
  }
  lowered <- am_crossprod(y, t(lhs))
  .amatrix_transpose_dense_result(lowered, lhs)
}

.amatrix_lower_sparse_rhs_crossprod <- function(x, y, ...) {
  lhs <- .amatrix_retarget_for_sparse_rhs(x, y, vector_as = "col")
  lowered <- am_crossprod(y, lhs, ...)
  .amatrix_transpose_dense_result(lowered, lhs)
}

.amatrix_lower_sparse_rhs_tcrossprod <- function(x, y, ...) {
  lhs <- .amatrix_retarget_for_sparse_rhs(x, y, vector_as = "row")
  lowered <- am_tcrossprod(y, lhs, ...)
  .amatrix_transpose_dense_result(lowered, lhs)
}

.amatrix_prepare_resident_arg <- function(value, backend_name, promote_amatrix = TRUE) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_residency_capable(backend)) {
    return(NULL)
  }

  if (inherits(value, "adgCMatrix") && is.function(backend$sparse_resident_store)) {
    current_backend <- .amatrix_live_resident_backend(value)
    resident_key <- .amatrix_resident_key(value, backend = backend_name)
    if (!is.null(resident_key) && isTRUE(backend$sparse_resident_has(resident_key))) {
      return(list(key = resident_key, temporary = FALSE, tracked = TRUE, sparse = TRUE))
    }
    if (!isTRUE(promote_amatrix)) {
      return(NULL)
    }
    resident_key <- .amatrix_next_resident_key(backend_name)
    host <- amatrix_materialize_host(value)  # returns dgCMatrix
    backend$sparse_resident_store(resident_key, host)
    if (!is.null(current_backend) && !identical(current_backend, backend_name)) {
      return(list(key = resident_key, temporary = TRUE, tracked = FALSE, sparse = TRUE))
    }
    .amatrix_bind_resident(value, backend_name, resident_key, sparse = TRUE)
    return(list(key = resident_key, temporary = FALSE, tracked = TRUE, sparse = TRUE))
  }

  if (inherits(value, "aTransposeView")) {
    if (is.function(backend$transpose_resident)) {
      source <- .amatrix_prepare_resident_arg(value@source, backend_name, promote_amatrix = promote_amatrix)
      if (is.null(source) || isTRUE(source$sparse)) {
        .amatrix_cleanup_temp_resident(list(source), backend_name)
        return(NULL)
      }

      resident_key <- .amatrix_next_resident_key(backend_name)
      ok <- tryCatch(
        {
          backend$transpose_resident(source$key, resident_key)
          TRUE
        },
        error = function(e) {
          try(backend$resident_drop(resident_key), silent = TRUE)
          FALSE
        }
      )
      .amatrix_cleanup_temp_resident(list(source), backend_name)
      if (ok) {
        return(list(key = resident_key, temporary = TRUE, tracked = FALSE))
      }
    }

    if (!isTRUE(promote_amatrix)) {
      return(NULL)
    }

    resident_key <- .amatrix_next_resident_key(backend_name)
    backend$resident_store(resident_key, t(as.matrix(amatrix_materialize_dense(value@source))))
    return(list(key = resident_key, temporary = TRUE, tracked = FALSE))
  }

  if (inherits(value, "adgeMatrix")) {
    current_backend <- .amatrix_live_resident_backend(value)
    resident_key <- .amatrix_resident_key(value, backend = backend_name)
    if (!is.null(resident_key) && isTRUE(backend$resident_has(resident_key))) {
      return(list(key = resident_key, temporary = FALSE, tracked = TRUE))
    }
    if (!isTRUE(promote_amatrix)) {
      return(NULL)
    }

    resident_key <- .amatrix_next_resident_key(backend_name)
    backend$resident_store(resident_key, amatrix_materialize_host(value))
    if (!is.null(current_backend) && !identical(current_backend, backend_name)) {
      return(list(key = resident_key, temporary = TRUE, tracked = FALSE))
    }
    .amatrix_bind_resident(value, backend_name, resident_key)
    return(list(key = resident_key, temporary = FALSE, tracked = TRUE))
  }

  if (inherits(value, "aMatrix")) {
    return(NULL)
  }

  if (.amatrix_is_dense_matrix_like(value)) {
    resident_key <- .amatrix_next_resident_key(backend_name)
    backend$resident_store(resident_key, .amatrix_host_arg(value))
    return(list(key = resident_key, temporary = TRUE, tracked = FALSE))
  }

  NULL
}

.amatrix_cleanup_temp_resident <- function(args, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_residency_capable(backend)) {
    return(invisible(NULL))
  }

  for (arg in args) {
    if (!is.null(arg) && isTRUE(arg$temporary)) {
      if (isTRUE(arg$sparse) && is.function(backend$sparse_resident_has) &&
          isTRUE(backend$sparse_resident_has(arg$key))) {
        backend$sparse_resident_drop(arg$key)
      } else if (isTRUE(backend$resident_has(arg$key))) {
        backend$resident_drop(arg$key)
      }
    }
  }

  invisible(NULL)
}

.amatrix_release_resident_key <- function(backend_name, key, sparse = FALSE) {
  if (is.null(key) || !nzchar(key)) {
    return(invisible(NULL))
  }

  backend <- tryCatch(.amatrix_get_backend(backend_name), error = function(e) NULL)
  if (is.null(backend) || !.amatrix_backend_residency_capable(backend)) {
    return(invisible(NULL))
  }

  if (isTRUE(sparse)) {
    if (is.function(backend$sparse_resident_has) &&
        is.function(backend$sparse_resident_drop) &&
        isTRUE(backend$sparse_resident_has(key))) {
      try(backend$sparse_resident_drop(key), silent = TRUE)
    }
    return(invisible(NULL))
  }

  if (is.function(backend$resident_has) &&
      is.function(backend$resident_drop) &&
      isTRUE(backend$resident_has(key))) {
    try(backend$resident_drop(key), silent = TRUE)
  }

  invisible(NULL)
}

.amatrix_cleanup_temp_resident_safe <- function(args, backend_name) {
  try(.amatrix_cleanup_temp_resident(args, backend_name), silent = TRUE)
  invisible(NULL)
}

# ── Deferred host materialization helpers ──────────────────────────────────
# When active, resident ops skip the GPU→CPU download.  The result is wrapped
# in a deferred adgeMatrix whose @x holds NaN until the first host access
# triggers a transparent download via amatrix_materialize_dense().

# Returns TRUE when the caller should skip host materialization.
# Default remains FALSE for the public adgeMatrix path. Backends can opt into
# the deferred resident result path by setting the global option during
# benchmarks or fused pipelines.
.amatrix_defer_host_active <- function() {
  isTRUE(getOption("amatrix.defer_host", FALSE))
}

# Wrap a resident op result as either a deferred or eager adgeMatrix.
# `template` is the source adgeMatrix (provides backend/policy/precision).
# `resident` is the list(value, backend, resident_key) from _try_resident_*.
# `out_dim` is required when resident$value is NULL (deferred mode).
.amatrix_resident_wrap <- function(template, resident, out_dim = NULL) {
  release_output <- !is.null(resident$backend) &&
    !is.null(resident$resident_key) &&
    nzchar(resident$resident_key)
  on.exit({
    if (release_output) {
      .amatrix_release_resident_key(resident$backend, resident$resident_key)
    }
  }, add = TRUE)

  if (is.null(resident$value) && !is.null(out_dim)) {
    # Deferred path
    dn <- if (inherits(template, "aMatrix")) template@Dimnames else list(NULL, NULL)
    # Dims may have changed (e.g., crossprod changes shape) — only reuse
    # dimnames if they match the output shape.
    if (!identical(as.integer(out_dim), template@Dim)) dn <- list(NULL, NULL)
    obj <- new_adgeMatrix_deferred(
      dim      = out_dim,
      dimnames = dn,
      preferred_backend = template@preferred_backend,
      policy    = template@policy,
      precision = template@precision
    )
    bound <- .amatrix_bind_resident(obj, resident$backend, resident$resident_key)
    release_output <- FALSE
    return(bound)
  }
  # Eager path
  value <- .amatrix_rewrap_like(template, resident$value)
  bound <- .amatrix_bind_resident(value, resident$backend, resident$resident_key)
  release_output <- FALSE
  bound
}

.amatrix_try_resident_matmul <- function(x, y, backend_name) {
  backend <- .amatrix_get_backend(backend_name)

  # Sparse resident path: LHS is adgCMatrix with spmm_resident support
  if (inherits(x, "adgCMatrix")) {
    if (!.amatrix_backend_supports_resident_op(backend, "matmul", x = x, y = y)) {
      return(NULL)
    }
    lhs <- .amatrix_prepare_resident_arg(x, backend_name, promote_amatrix = TRUE)
    if (is.null(lhs) || !isTRUE(lhs$sparse)) return(NULL)

    if (is.function(backend$spmm_resident_key)) {
      rhs <- .amatrix_prepare_resident_arg(y, backend_name, promote_amatrix = TRUE)
      if (!is.null(rhs) && !isTRUE(rhs$sparse)) {
        defer <- .amatrix_defer_host_active()
        out_key <- .amatrix_next_resident_key(backend_name)
        ok <- FALSE
        value <- tryCatch(
          {
            ok <- TRUE
            backend$spmm_resident_key(lhs$key, rhs$key, out_key, trans_lhs = FALSE, defer = defer)
          },
          error = function(e) {
            try(backend$resident_drop(out_key), silent = TRUE)
            NULL
          }
        )
        .amatrix_cleanup_temp_resident(list(rhs), backend_name)
        if (isTRUE(ok)) {
          return(list(value = value, backend = backend_name, resident_key = out_key))
        }
      } else {
        .amatrix_cleanup_temp_resident(list(rhs), backend_name)
      }
    }

    if (is.function(backend$spmm_resident)) {
      rhs_mat <- if (is.matrix(y)) y else as.matrix(.amatrix_host_arg(y))
      if (!is.double(rhs_mat)) storage.mode(rhs_mat) <- "double"
      value <- backend$spmm_resident(lhs$key, rhs_mat, trans_lhs = FALSE)
      return(list(value = value, backend = NULL, resident_key = NULL, host_only = TRUE))
    }
  }

  if (!.amatrix_backend_supports_resident_op(backend, "matmul", x = x, y = y)) {
    return(NULL)
  }

  lhs <- .amatrix_prepare_resident_arg(x, backend_name, promote_amatrix = TRUE)
  rhs <- .amatrix_prepare_resident_arg(y, backend_name, promote_amatrix = TRUE)
  if (is.null(lhs) || is.null(rhs)) {
    .amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name)
    return(NULL)
  }

  defer <- .amatrix_defer_host_active()
  out_key <- .amatrix_next_resident_key(backend_name)
  release_output <- TRUE
  on.exit({
    if (release_output) {
      .amatrix_release_resident_key(backend_name, out_key)
    }
    .amatrix_cleanup_temp_resident_safe(list(lhs, rhs), backend_name)
  }, add = TRUE)
  value <- tryCatch(
    backend$matmul_resident(lhs$key, rhs$key, out_key, defer = defer),
    error = function(e) NULL
  )
  if (is.null(value)) return(NULL)

  release_output <- FALSE
  list(value = value, backend = backend_name, resident_key = out_key)
}

.amatrix_try_dense_sparse_resident_matmul <- function(x, y, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!inherits(x, "adgeMatrix") || !inherits(y, "adgCMatrix") ||
      !is.function(backend$dense_sparse_matmul_resident_key)) {
    return(NULL)
  }

  lhs <- .amatrix_prepare_resident_arg(x, backend_name, promote_amatrix = TRUE)
  rhs <- .amatrix_prepare_resident_arg(y, backend_name, promote_amatrix = TRUE)
  if (is.null(lhs) || is.null(rhs) || isTRUE(rhs$sparse) != TRUE) {
    .amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name)
    return(NULL)
  }

  defer <- .amatrix_defer_host_active()
  out_key <- .amatrix_next_resident_key(backend_name)
  release_output <- TRUE
  on.exit({
    if (release_output) {
      .amatrix_release_resident_key(backend_name, out_key)
    }
    .amatrix_cleanup_temp_resident_safe(list(lhs, rhs), backend_name)
  }, add = TRUE)
  value <- tryCatch(
    backend$dense_sparse_matmul_resident_key(lhs$key, rhs$key, out_key, defer = defer),
    error = function(e) NULL
  )
  if (is.null(value)) {
    return(NULL)
  }

  release_output <- FALSE
  list(value = value, backend = backend_name, resident_key = out_key)
}

.amatrix_try_resident_crossprod <- function(x, y, backend_name) {
  backend <- .amatrix_get_backend(backend_name)

  # Sparse resident path: crossprod(X, Y) = t(X) %*% Y → spmm_resident with trans_lhs=TRUE
  if (inherits(x, "adgCMatrix") && !is.null(y)) {
    if (!.amatrix_backend_supports_resident_op(backend, "crossprod", x = x, y = y)) {
      return(NULL)
    }
    lhs <- .amatrix_prepare_resident_arg(x, backend_name, promote_amatrix = TRUE)
    if (is.null(lhs) || !isTRUE(lhs$sparse)) return(NULL)

    if (is.function(backend$spmm_resident_key)) {
      rhs <- .amatrix_prepare_resident_arg(y, backend_name, promote_amatrix = TRUE)
      if (!is.null(rhs) && !isTRUE(rhs$sparse)) {
        defer <- .amatrix_defer_host_active()
        out_key <- .amatrix_next_resident_key(backend_name)
        release_output <- TRUE
        on.exit({
          if (release_output) {
            .amatrix_release_resident_key(backend_name, out_key)
          }
          .amatrix_cleanup_temp_resident_safe(list(lhs, rhs), backend_name)
        }, add = TRUE)
        value <- tryCatch(
          backend$spmm_resident_key(lhs$key, rhs$key, out_key, trans_lhs = TRUE, defer = defer),
          error = function(e) NULL
        )
        if (!is.null(value)) {
          release_output <- FALSE
          return(list(value = value, backend = backend_name, resident_key = out_key))
        }
      } else {
        .amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name)
      }
    }

    if (is.function(backend$spmm_resident)) {
      rhs_mat <- if (is.matrix(y)) y else as.matrix(.amatrix_host_arg(y))
      if (!is.double(rhs_mat)) storage.mode(rhs_mat) <- "double"
      value <- backend$spmm_resident(lhs$key, rhs_mat, trans_lhs = TRUE)
      return(list(value = value, backend = NULL, resident_key = NULL, host_only = TRUE))
    }

    .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  }

  if (!.amatrix_backend_supports_resident_op(backend, "crossprod", x = x, y = y)) {
    return(NULL)
  }

  lhs <- .amatrix_prepare_resident_arg(x, backend_name, promote_amatrix = TRUE)
  rhs <- if (is.null(y)) NULL else .amatrix_prepare_resident_arg(y, backend_name, promote_amatrix = TRUE)
  if (is.null(lhs) || (!is.null(y) && is.null(rhs))) {
    .amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name)
    return(NULL)
  }

  defer <- .amatrix_defer_host_active()
  out_key <- .amatrix_next_resident_key(backend_name)
  rhs_key <- if (is.null(rhs)) NULL else rhs$key
  release_output <- TRUE
  on.exit({
    if (release_output) {
      .amatrix_release_resident_key(backend_name, out_key)
    }
    .amatrix_cleanup_temp_resident_safe(list(lhs, rhs), backend_name)
  }, add = TRUE)
  value <- tryCatch(
    backend$crossprod_resident(lhs$key, rhs_key, out_key, defer = defer),
    error = function(e) NULL
  )
  if (is.null(value)) return(NULL)

  release_output <- FALSE
  list(value = value, backend = backend_name, resident_key = out_key)
}

.amatrix_try_resident_tcrossprod <- function(x, y, backend_name) {
  backend <- .amatrix_get_backend(backend_name)

  # Sparse resident path: tcrossprod(X, Y) = X %*% t(Y) → spmm_resident with
  # the dense RHS transposed on the host side.
  if (inherits(x, "adgCMatrix") && !is.null(y)) {
    if (!.amatrix_backend_supports_resident_op(backend, "tcrossprod", x = x, y = y)) {
      return(NULL)
    }
    lhs <- .amatrix_prepare_resident_arg(x, backend_name, promote_amatrix = TRUE)
    if (is.null(lhs) || !isTRUE(lhs$sparse)) return(NULL)

    if (is.function(backend$spmm_resident_key)) {
      rhs_t <- .amatrix_prepare_resident_arg(am_transpose(y), backend_name, promote_amatrix = TRUE)
      if (!is.null(rhs_t) && !isTRUE(rhs_t$sparse)) {
        defer <- .amatrix_defer_host_active()
        out_key <- .amatrix_next_resident_key(backend_name)
        release_output <- TRUE
        on.exit({
          if (release_output) {
            .amatrix_release_resident_key(backend_name, out_key)
          }
          .amatrix_cleanup_temp_resident_safe(list(lhs, rhs_t), backend_name)
        }, add = TRUE)
        value <- tryCatch(
          backend$spmm_resident_key(lhs$key, rhs_t$key, out_key, trans_lhs = FALSE, defer = defer),
          error = function(e) NULL
        )
        if (!is.null(value)) {
          release_output <- FALSE
          return(list(value = value, backend = backend_name, resident_key = out_key))
        }
      } else {
        .amatrix_cleanup_temp_resident(list(lhs, rhs_t), backend_name)
      }
    }

    if (is.function(backend$spmm_resident)) {
      rhs_mat <- if (is.matrix(y)) y else as.matrix(.amatrix_host_arg(y))
      if (!is.double(rhs_mat)) storage.mode(rhs_mat) <- "double"
      value <- backend$spmm_resident(lhs$key, t(rhs_mat), trans_lhs = FALSE)
      return(list(value = value, backend = NULL, resident_key = NULL, host_only = TRUE))
    }

    .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  }

  if (!.amatrix_backend_supports_resident_op(backend, "tcrossprod", x = x, y = y)) {
    return(NULL)
  }

  lhs <- .amatrix_prepare_resident_arg(x, backend_name, promote_amatrix = TRUE)
  rhs <- if (is.null(y)) NULL else .amatrix_prepare_resident_arg(y, backend_name, promote_amatrix = TRUE)
  if (is.null(lhs) || (!is.null(y) && is.null(rhs))) {
    .amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name)
    return(NULL)
  }

  defer <- .amatrix_defer_host_active()
  out_key <- .amatrix_next_resident_key(backend_name)
  rhs_key <- if (is.null(rhs)) NULL else rhs$key
  release_output <- TRUE
  on.exit({
    if (release_output) {
      .amatrix_release_resident_key(backend_name, out_key)
    }
    .amatrix_cleanup_temp_resident_safe(list(lhs, rhs), backend_name)
  }, add = TRUE)
  value <- tryCatch(
    backend$tcrossprod_resident(lhs$key, rhs_key, out_key, defer = defer),
    error = function(e) NULL
  )
  if (is.null(value)) return(NULL)

  release_output <- FALSE
  list(value = value, backend = backend_name, resident_key = out_key)
}

.amatrix_try_resident_ewise <- function(op, e1, e2, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  template <- .amatrix_template(e1, e2)
  if (is.null(template) || !inherits(template, "adgeMatrix")) {
    return(NULL)
  }
  rhs_arg <- if (inherits(e1, "adgeMatrix")) e2 else e1
  if (!.amatrix_backend_supports_resident_op(backend, "ewise", x = template, y = rhs_arg)) {
    return(NULL)
  }

  lhs <- if (inherits(e1, "adgeMatrix")) {
    .amatrix_prepare_resident_arg(e1, backend_name, promote_amatrix = TRUE)
  } else {
    .amatrix_prepare_resident_arg(e2, backend_name, promote_amatrix = TRUE)
  }
  if (is.null(lhs)) {
    return(NULL)
  }

  rhs <- NULL
  rhs_payload <- rhs_arg

  if (inherits(rhs_arg, "adgeMatrix") || .amatrix_is_dense_matrix_like(rhs_arg)) {
    rhs <- .amatrix_prepare_resident_arg(rhs_arg, backend_name, promote_amatrix = TRUE)
    if (is.null(rhs)) {
      .amatrix_cleanup_temp_resident(list(lhs), backend_name)
      return(NULL)
    }
    rhs_payload <- rhs$key
  } else if (!is.null(rhs_arg) && is.numeric(rhs_arg) && length(rhs_arg) == 1L) {
    rhs_payload <- as.double(rhs_arg)
  } else if (!is.null(rhs_arg)) {
    .amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name)
    return(NULL)
  }

  defer <- .amatrix_defer_host_active()
  out_key <- .amatrix_next_resident_key(backend_name)
  release_output <- TRUE
  on.exit({
    if (release_output) {
      .amatrix_release_resident_key(backend_name, out_key)
    }
    .amatrix_cleanup_temp_resident_safe(list(lhs, rhs), backend_name)
  }, add = TRUE)
  value <- tryCatch(
    backend$ewise_resident(lhs$key, rhs_payload, op, out_key, defer = defer),
    error = function(e) NULL
  )

  if (is.null(value)) {
    return(NULL)
  }

  release_output <- FALSE
  list(value = value, backend = backend_name, resident_key = out_key)
}

# rowSums/colSums: output is a vector — no resident binding on result.
.amatrix_try_resident_rowSums <- function(x, na.rm, dims, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_supports_resident_op(backend, "rowSums", x = x)) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(x, backend_name, promote_amatrix = FALSE)
  if (is.null(lhs)) return(NULL)
  result <- backend$rowSums_resident(lhs$key, na.rm, dims)
  .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  result
}

.amatrix_try_resident_colSums <- function(x, na.rm, dims, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_supports_resident_op(backend, "colSums", x = x)) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(x, backend_name, promote_amatrix = FALSE)
  if (is.null(lhs)) return(NULL)
  result <- backend$colSums_resident(lhs$key, na.rm, dims)
  .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  result
}

# am_solve: output is a matrix; store at out_key and bind resident.
.amatrix_try_resident_solve <- function(a, b, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_supports_resident_op(backend, "solve", x = a, y = b)) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(a, backend_name, promote_amatrix = TRUE)
  if (is.null(lhs)) return(NULL)
  out_key <- .amatrix_next_resident_key(backend_name)
  release_output <- TRUE
  if (is.null(b)) {
    on.exit({
      if (release_output) {
        .amatrix_release_resident_key(backend_name, out_key)
      }
      .amatrix_cleanup_temp_resident_safe(list(lhs), backend_name)
    }, add = TRUE)
    value <- tryCatch(
      backend$solve_resident(lhs$key, NULL, out_key),
      error = function(e) NULL
    )
  } else {
    b_arg <- if (is.vector(b)) matrix(b, ncol = 1L) else b
    rhs <- .amatrix_prepare_resident_arg(b_arg, backend_name, promote_amatrix = TRUE)
    if (is.null(rhs)) {
      .amatrix_cleanup_temp_resident(list(lhs), backend_name)
      return(NULL)
    }
    on.exit({
      if (release_output) {
        .amatrix_release_resident_key(backend_name, out_key)
      }
      .amatrix_cleanup_temp_resident_safe(list(lhs, rhs), backend_name)
    }, add = TRUE)
    value <- tryCatch(
      backend$solve_resident(lhs$key, rhs$key, out_key),
      error = function(e) NULL
    )
  }
  if (is.null(value)) return(NULL)
  release_output <- FALSE
  list(value = value, backend = backend_name, resident_key = out_key)
}

# chol: output is a matrix; store at out_key and bind resident.
.amatrix_try_resident_chol <- function(x, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_supports_resident_op(backend, "chol", x = x)) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(x, backend_name, promote_amatrix = TRUE)
  if (is.null(lhs)) return(NULL)
  out_key <- .amatrix_next_resident_key(backend_name)
  release_output <- TRUE
  on.exit({
    if (release_output) {
      .amatrix_release_resident_key(backend_name, out_key)
    }
    .amatrix_cleanup_temp_resident_safe(list(lhs), backend_name)
  }, add = TRUE)
  value <- tryCatch(
    backend$chol_resident(lhs$key, out_key),
    error = function(e) NULL
  )
  if (is.null(value)) return(NULL)
  release_output <- FALSE
  list(value = value, backend = backend_name, resident_key = out_key)
}

#' Matrix multiplication
#'
#' Multiplies two matrices, routing to an accelerated backend when available.
#' Plain numeric vectors supplied as \code{y} are promoted to a column matrix
#' and the result is dropped back to a vector.
#'
#' @param x A matrix or \code{aMatrix} object.
#' @param y A matrix, \code{aMatrix} object, or numeric vector.
#'
#' @return A matrix (or numeric vector when \code{y} was a vector) of
#'   dimensions \code{nrow(x)} by \code{ncol(y)}.
#'
#' @examples
#' A <- adgeMatrix(matrix(1:6, 2, 3))
#' B <- adgeMatrix(matrix(1:6, 3, 2))
#' am_matmul(A, B)
#'
#' @export
am_matmul <- function(x, y) {
  x <- .amatrix_promote_arg(x)
  if (.amatrix_sparse_rhs_lowerable(x, y)) {
    return(.amatrix_lower_sparse_rhs_matmul(x, y))
  }

  # irlba's hot path passes a plain numeric vector for v (A %*% v).
  # _prepare_resident_arg only accepts matrices, so without promotion the
  # resident path silently fails: _try_resident_matmul returns NULL,
  # amatrix_dispatch_op drops the resident binding, and A is re-uploaded on
  # every Lanczos step.  Promoting to a column matrix fixes that.
  y_vec <- is.numeric(y) && is.null(dim(y))
  y_eff <- if (y_vec) matrix(y, ncol = 1L) else y

  choice <- .amatrix_backend_for(x, "matmul", y = y_eff)
  resident <- .amatrix_try_resident_matmul(x, y_eff, choice$name)
  if (!is.null(resident)) {
    if (isTRUE(resident$host_only)) {
      return(.amatrix_rewrap_like(x, resident$value))
    }
    return(.amatrix_resident_wrap(x, resident,
                                   out_dim = c(nrow(x), ncol(y_eff))))
  }

  if (identical(choice$name, "cpu")) {
    cpu_value <- choice$backend$matmul(x, y_eff)
    return(.amatrix_rewrap_like(x, cpu_value))
  }

  .amatrix_rewrap_like(
    x,
    amatrix_dispatch_op(
      x = x,
      op = "matmul",
      method = "matmul",
      y = y,
      args = list(y = .amatrix_host_arg(y)),
      fallback = function() amatrix_materialize_host(x) %*% .amatrix_host_arg(y)
    )
  )
}

am_crossprod <- function(x, y = NULL, ...) {
  if (.amatrix_sparse_rhs_lowerable(x, y)) {
    return(.amatrix_lower_sparse_rhs_crossprod(x, y, ...))
  }

  if (!inherits(x, "aMatrix")) {
    if (is.null(y)) return(base::crossprod(x, ...))
    return(base::crossprod(x, y = y, ...))
  }
  choice <- .amatrix_backend_for(x, "crossprod", y = y)
  resident <- .amatrix_try_resident_crossprod(x, y, choice$name)
  if (!is.null(resident)) {
    if (isTRUE(resident$host_only)) {
      return(.amatrix_rewrap_like(x, resident$value))
    }
    nc <- if (is.null(y)) ncol(x) else ncol(y)
    return(.amatrix_resident_wrap(x, resident, out_dim = c(ncol(x), nc)))
  }

  if (identical(choice$name, "cpu")) {
    return(.amatrix_rewrap_like(x, choice$backend$crossprod(x, y = y, ...)))
  }

  .amatrix_rewrap_like(
    x,
    amatrix_dispatch_op(
      x = x,
      op = "crossprod",
      method = "crossprod",
      y = y,
      args = list(y = .amatrix_host_arg(y), ...),
      fallback = function() {
        if (is.null(y)) {
          return(base::crossprod(as.matrix(amatrix_materialize_host(x)), ...))
        }
        base::crossprod(as.matrix(amatrix_materialize_host(x)), y = .amatrix_host_arg(y), ...)
      }
    )
  )
}

am_tcrossprod <- function(x, y = NULL, ...) {
  if (.amatrix_sparse_rhs_lowerable(x, y)) {
    return(.amatrix_lower_sparse_rhs_tcrossprod(x, y, ...))
  }

  if (!inherits(x, "aMatrix")) {
    if (is.null(y)) return(base::tcrossprod(x, ...))
    return(base::tcrossprod(x, y = y, ...))
  }
  choice <- .amatrix_backend_for(x, "tcrossprod", y = y)
  resident <- .amatrix_try_resident_tcrossprod(x, y, choice$name)
  if (!is.null(resident)) {
    if (isTRUE(resident$host_only)) {
      return(.amatrix_rewrap_like(x, resident$value))
    }
    nr <- if (is.null(y)) nrow(x) else nrow(y)
    return(.amatrix_resident_wrap(x, resident, out_dim = c(nrow(x), nr)))
  }

  if (identical(choice$name, "cpu")) {
    return(.amatrix_rewrap_like(x, choice$backend$tcrossprod(x, y = y, ...)))
  }

  .amatrix_rewrap_like(
    x,
    amatrix_dispatch_op(
      x = x,
      op = "tcrossprod",
      method = "tcrossprod",
      y = y,
      args = list(y = .amatrix_host_arg(y), ...),
      fallback = function() {
        if (is.null(y)) {
          return(base::tcrossprod(as.matrix(amatrix_materialize_host(x)), ...))
        }
        base::tcrossprod(as.matrix(amatrix_materialize_host(x)), y = .amatrix_host_arg(y), ...)
      }
    )
  )
}

#' Generalised matrix multiply (BLAS DGEMM interface)
#'
#' Computes \code{alpha * op(A) \%*\% op(B) + beta * C}, where
#' \code{op(X) = t(X)} when the corresponding \code{trans} flag is
#' \code{TRUE}. Routes internally to the most efficient resident
#' operation for the chosen transpose combination.
#'
#' @param A A matrix or \code{aMatrix}.
#' @param B A matrix or \code{aMatrix}.
#' @param C Optional matrix or \code{aMatrix} to add after scaling;
#'   \code{NULL} omits the addition term.
#' @param alpha Numeric scalar multiplier for \code{op(A) \%*\% op(B)}.
#'   Default \code{1.0}.
#' @param beta Numeric scalar multiplier for \code{C}. Default \code{1.0}.
#' @param transA Logical; transpose \code{A} before multiplying.
#'   Default \code{FALSE}.
#' @param transB Logical; transpose \code{B} before multiplying.
#'   Default \code{FALSE}.
#'
#' @return A matrix of dimensions \code{nrow(op(A))} by
#'   \code{ncol(op(B))}.
#'
#' @examples
#' A <- adgeMatrix(matrix(1:6, 2, 3))
#' B <- adgeMatrix(matrix(1:6, 2, 3))
#' am_gemm(A, B, transA = TRUE)          # t(A) %*% B
#' am_gemm(A, B, transB = TRUE)          # A %*% t(B)
#'
#' @export
am_gemm <- function(A, B, C = NULL, alpha = 1.0, beta = 1.0,
                    transA = FALSE, transB = FALSE) {
  A <- .amatrix_promote_arg(A)
  if (!is.null(C)) C <- .amatrix_promote_arg(C)
  AB <- if (transA && transB) {
    # t(A) %*% t(B) = t(B %*% A): use identity to avoid materialising t(B) to host.
    am_transpose(am_matmul(B, A))
  } else if (transA) {
    am_crossprod(A, B)
  } else if (transB) {
    am_tcrossprod(A, B)
  } else {
    am_matmul(A, B)
  }

  if (alpha != 1.0) AB <- am_ewise("*", AB, alpha)

  if (!is.null(C)) {
    C_scaled <- if (beta != 1.0) am_ewise("*", C, beta) else C
    am_ewise("+", AB, C_scaled)
  } else {
    AB
  }
}

#' Row and column sums
#'
#' Compute row or column sums of a matrix or \code{aMatrix}, dispatching
#' to an accelerated backend when one is available.
#'
#' @param x A matrix or \code{aMatrix} object.
#' @param na.rm Logical; if \code{TRUE}, missing values are removed before
#'   summing. Default \code{FALSE}.
#' @param dims Integer; the number of dimensions to regard as rows
#'   (for \code{am_rowsums}) or columns (for \code{am_colsums}). Default \code{1L}.
#'
#' @return A numeric vector of length \code{nrow(x)} (\code{am_rowsums}) or
#'   \code{ncol(x)} (\code{am_colsums}).
#'
#' @examples
#' m <- adgeMatrix(matrix(1:12, 3, 4))
#' am_rowsums(m)
#' am_colsums(m)
#'
#' @export
am_rowsums <- function(x, na.rm = FALSE, dims = 1L) {
  x <- .amatrix_promote_arg(x)
  choice <- .amatrix_backend_for(x, "rowSums")
  resident <- .amatrix_try_resident_rowSums(x, na.rm, dims, choice$name)
  if (!is.null(resident)) return(resident)
  amatrix_dispatch_op(
    x = x,
    op = "rowSums",
    method = "rowSums",
    args = list(na.rm = na.rm, dims = dims),
    fallback = function() Matrix::rowSums(amatrix_materialize_host(x), na.rm = na.rm, dims = dims)
  )
}

#' @rdname am_rowsums
#' @export
am_colsums <- function(x, na.rm = FALSE, dims = 1L) {
  x <- .amatrix_promote_arg(x)
  choice <- .amatrix_backend_for(x, "colSums")
  resident <- .amatrix_try_resident_colSums(x, na.rm, dims, choice$name)
  if (!is.null(resident)) return(resident)
  amatrix_dispatch_op(
    x = x,
    op = "colSums",
    method = "colSums",
    args = list(na.rm = na.rm, dims = dims),
    fallback = function() Matrix::colSums(amatrix_materialize_host(x), na.rm = na.rm, dims = dims)
  )
}

am_transpose <- function(x) {
  if (inherits(x, "adgeMatrix")) {
    return(.new_aTransposeView(x))
  }
  if (inherits(x, "adgCMatrix")) {
    host <- amatrix_materialize_host(x)   # dgCMatrix
    return(new_adgCMatrix(t(host), preferred_backend = x@preferred_backend,
                          precision = x@precision, policy = x@policy))
  }
  .amatrix_rewrap_like(x, t(as.matrix(amatrix_materialize_host(x))))
}

# `.linear` is TRUE for the one-subscript form `x[i]` (no comma), which the
# calling S4 method detects with nargs(). It selects elements in column-major
# order, exactly as for a base matrix, rather than rows.
am_subset <- function(x, i, j, ..., drop = TRUE, .linear = FALSE) {
  host_x <- amatrix_materialize_host(x)
  i_arg <- if (!missing(i)) .amatrix_subscript_arg(i) else NULL
  j_arg <- if (!missing(j)) .amatrix_subscript_arg(j) else NULL
  value <- if (missing(i) && missing(j)) {
    host_x[..., drop = drop]
  } else if (missing(j)) {
    if (isTRUE(.linear) || .amatrix_is_logical_matrix_subscript(i_arg, host_x)) {
      host_x[i_arg]
    } else {
      host_x[i_arg, , ..., drop = drop]
    }
  } else if (missing(i)) {
    host_x[, j_arg, ..., drop = drop]
  } else {
    host_x[i_arg, j_arg, ..., drop = drop]
  }
  .amatrix_rewrap_value(x, value)
}

am_subassign <- function(x, i, j, ..., value, .linear = FALSE) {
  host_x <- amatrix_materialize_host(x)
  host_value <- .amatrix_host_arg(value)
  i_arg <- if (!missing(i)) .amatrix_subscript_arg(i) else NULL
  j_arg <- if (!missing(j)) .amatrix_subscript_arg(j) else NULL
  if (missing(i) && missing(j)) {
    host_x[...] <- host_value
  } else if (missing(j)) {
    if (isTRUE(.linear) || .amatrix_is_logical_matrix_subscript(i_arg, host_x)) {
      host_x[i_arg] <- host_value
    } else {
      host_x[i_arg, , ...] <- host_value
    }
  } else if (missing(i)) {
    host_x[, j_arg, ...] <- host_value
  } else {
    host_x[i_arg, j_arg, ...] <- host_value
  }
  .amatrix_cache_invalidate_object(x@object_id)
  .amatrix_rewrap_value(x, host_x)
}

am_solve <- function(a, b = NULL, ...) {
  b_arg <- if (missing(b)) NULL else b
  if (!inherits(a, "aMatrix")) {
    if (is.null(b_arg)) return(base::solve(a, ...))
    return(base::solve(a, b = b_arg, ...))
  }
  if (inherits(a, "adgCMatrix")) {
    host <- amatrix_materialize_host(a)
    if (is.null(b_arg)) {
      result <- Matrix::solve(host, ...)
      if (inherits(result, "sparseMatrix")) return(new_adgCMatrix(result, preferred_backend = a@preferred_backend, precision = a@precision, policy = a@policy))
      return(as_adgeMatrix(as.matrix(result), preferred_backend = a@preferred_backend, precision = a@precision, policy = a@policy))
    }
    b_mat <- if (inherits(b_arg, "adgCMatrix")) amatrix_materialize_host(b_arg) else b_arg
    result <- Matrix::solve(host, b_mat, ...)
    if (inherits(result, "sparseMatrix")) return(new_adgCMatrix(result, preferred_backend = a@preferred_backend, precision = a@precision, policy = a@policy))
    return(as_adgeMatrix(as.matrix(result), preferred_backend = a@preferred_backend, precision = a@precision, policy = a@policy))
  }
  b_was_vector <- !is.null(b_arg) && is.numeric(b_arg) && is.null(dim(b_arg))
  choice <- .amatrix_backend_for(a, "solve", y = b_arg)
  resident <- .amatrix_try_resident_solve(a, b_arg, choice$name)
  if (!is.null(resident)) {
    if (b_was_vector) return(as.vector(resident$value))
    value <- .amatrix_rewrap_value(a, resident$value)
    return(.amatrix_bind_resident(value, resident$backend, resident$resident_key))
  }

  if (is.null(b_arg)) {
    return(.amatrix_rewrap_value(
      a,
      amatrix_dispatch_op(
        x = a,
        op = "solve",
        method = "solve",
        args = list(...),
        fallback = function() base::solve(as.matrix(amatrix_materialize_host(a)), ...)
      )
    ))
  }

  result <- amatrix_dispatch_op(
    x = a,
    op = "solve",
    method = "solve",
    y = b_arg,
    args = list(b = .amatrix_host_arg(b_arg), ...),
    fallback = function() base::solve(as.matrix(amatrix_materialize_host(a)), as.matrix(.amatrix_host_arg(b_arg)), ...)
  )
  # When b was a plain vector, base::solve() returns a named numeric vector.
  # Preserve that contract: don't wrap back to adgeMatrix.
  if (b_was_vector) {
    out <- as.vector(result)
    names(out) <- colnames(a)
    return(out)
  }
  .amatrix_rewrap_value(a, result)
}

am_chol <- function(x, ...) {
  if (inherits(x, "adgCMatrix")) {
    host <- amatrix_materialize_host(x)
    result <- Matrix::chol(host, ...)
    if (inherits(result, "sparseMatrix")) return(new_adgCMatrix(result, preferred_backend = x@preferred_backend, precision = x@precision, policy = x@policy))
    return(result)
  }
  choice <- .amatrix_backend_for(x, "chol")
  resident <- .amatrix_try_resident_chol(x, choice$name)
  if (!is.null(resident)) {
    value <- .amatrix_rewrap_value(x, resident$value)
    return(.amatrix_bind_resident(value, resident$backend, resident$resident_key))
  }
  .amatrix_rewrap_value(
    x,
    amatrix_dispatch_op(
      x = x,
      op = "chol",
      method = "chol",
      args = list(...),
      fallback = function() base::chol(as.matrix(amatrix_materialize_host(x)), ...)
    )
  )
}

.amatrix_qr_arg <- function(x) {
  if (inherits(x, "adgCMatrix") || inherits(x, "adgeMatrix") || inherits(x, "amQR")) {
    return(x)
  }

  if (inherits(x, "aMatrix")) {
    return(x)
  }

  if (inherits(x, "dgCMatrix") || (inherits(x, "sparseMatrix") && !inherits(x, "denseMatrix"))) {
    return(new_adgCMatrix(
      x,
      preferred_backend = "cpu",
      policy = amatrix_default_policy(),
      precision = amatrix_default_precision()
    ))
  }

  if (inherits(x, "Matrix") || is.matrix(x)) {
    return(as_adgeMatrix(
      as.matrix(x),
      preferred_backend = "cpu",
      policy = amatrix_default_policy(),
      precision = amatrix_default_precision()
    ))
  }

  stop("x must be a matrix-like object", call. = FALSE)
}

#' QR decomposition of an amatrix object
#'
#' Computes the QR decomposition of a matrix or \code{aMatrix}, routing to
#' a backend-specific implementation when available.
#'
#' @param x A matrix or \code{aMatrix} object.
#' @param ... Additional arguments passed to the underlying QR routine.
#'
#' @return An object of class \code{amDenseQR} (or a wrapped sparse QR
#'   for \code{adgCMatrix} input) containing the factorisation components.
#'
#' @examples
#' m <- adgeMatrix(matrix(rnorm(12), 4, 3))
#' qr_obj <- am_qr(m)
#' qr.R(qr_obj)
#'
#' @export
am_qr <- function(x, ...) {
  x <- .amatrix_qr_arg(x)

  if (inherits(x, "adgCMatrix")) {
    host <- amatrix_materialize_host(x)
    return(.amatrix_wrap_sparse_qr(Matrix::qr(host, ...), x, method = "cpu"))
  }
  qr_value <- amatrix_dispatch_op(
    x = x,
    op = "qr",
    method = "qr",
    args = list(...),
    fallback = function() base::qr(as.matrix(amatrix_materialize_host(x)), ...)
  )

  .amatrix_wrap_qr(qr_value, x)
}

# ── Resident QR helper ────────────────────────────────────────────────────────
#
# Given an adgeMatrix z_am that is already (or can be) resident on the backend,
# run QR and return Q as a new resident adgeMatrix — no R-memory round-trip for
# the QR step itself.
#
# Falls back to qr.Q(qr(as.matrix(z_am))) + re-upload when the backend does
# not support qr_Q_resident (e.g. the cpu backend).
.amatrix_try_resident_qr_Q <- function(z_am) {
  if (!inherits(z_am, "adgeMatrix")) {
    z_am <- adgeMatrix(as.matrix(z_am))
  }
  choice       <- .amatrix_backend_for(z_am, "qr")
  backend_name <- choice$name
  backend      <- choice$backend

  if (is.function(backend$qr_Q_resident) &&
      .amatrix_backend_residency_capable(backend)) {
    z_info <- .amatrix_prepare_resident_arg(z_am, backend_name)
    if (!is.null(z_info)) {
      q_key <- .amatrix_next_resident_key(backend_name)
      backend$qr_Q_resident(z_info$key, q_key)
      .amatrix_cleanup_temp_resident(list(z_info), backend_name)
      # Materialize Q for the host slot (dims/fallback paths).
      # Q is p×k_over — small relative to the data matrix.
      q_mat <- backend$resident_materialize(q_key)
      q_am  <- new_adgeMatrix(q_mat,
                  preferred_backend = backend_name,
                  policy            = z_am@policy,
                  precision         = z_am@precision)
      .amatrix_bind_resident(q_am, backend_name, q_key)
      return(q_am)
    }
  }
  {
    q_mat <- qr.Q(qr(as.matrix(amatrix_materialize_host(z_am))))
    adgeMatrix(q_mat,
      preferred_backend = backend_name,
      policy            = z_am@policy,
      precision         = z_am@precision)
  }
}

am_svd <- function(x, nu = min(dim(x)), nv = min(dim(x)), LINPACK = FALSE, ...) {
  amatrix_dispatch_op(
    x = x,
    op = "svd",
    method = "svd",
    args = list(nu = nu, nv = nv, LINPACK = LINPACK, ...),
    fallback = function() {
      if (isTRUE(LINPACK)) {
        stop("LINPACK is not supported", call. = FALSE)
      }
      base::svd(as.matrix(amatrix_materialize_host(x)), nu = nu, nv = nv, ...)
    }
  )
}

am_eigen <- function(x, symmetric = NULL, only.values = FALSE, EISPACK = FALSE) {
  # Mirror base::eigen behaviour: if symmetric is not supplied, auto-detect
  # from the host matrix so callers don't have to know the structure.
  if (is.null(symmetric)) {
    x_host <- as.matrix(amatrix_materialize_host(x))
    symmetric <- isSymmetric(x_host)
  }

  # Some GPU backends currently only expose a dense symmetric eigen surface.
  # Keep the nonsymmetric path honest and fall back to the host implementation
  # rather than pretending there is native support for the full general problem.
  if (!isTRUE(symmetric) && inherits(x, "aMatrix")) {
    choice <- .amatrix_backend_for(x, "eigen")
    backend_features <- tryCatch(choice$backend$features(), error = function(e) character())
    if (!identical(choice$name, "cpu") && "eigen_sym" %in% backend_features) {
      return(base::eigen(as.matrix(amatrix_materialize_host(x)), symmetric = FALSE, only.values = only.values, EISPACK = EISPACK))
    }
  }

  amatrix_dispatch_op(
    x = x,
    op = "eigen",
    method = "eigen",
    args = list(symmetric = symmetric, only.values = only.values, EISPACK = EISPACK),
    fallback = function() base::eigen(as.matrix(amatrix_materialize_host(x)), symmetric = symmetric, only.values = only.values, EISPACK = EISPACK)
  )
}

#' Symmetric eigendecomposition
#'
#' Computes eigenvalues and eigenvectors of a real symmetric matrix by
#' dispatching to the active backend via \code{\link[base]{eigen}} with
#' \code{symmetric = TRUE}.
#'
#' @param x A real symmetric numeric matrix, \code{adgeMatrix}, or other
#'   object accepted by \code{\link[base]{eigen}}.
#'
#' @return A list with components \code{values} (numeric vector of eigenvalues
#'   in decreasing order, as for \code{\link[base]{eigen}}) and \code{vectors} (numeric matrix whose columns are
#'   the corresponding eigenvectors).
#'
#' @examples
#' S <- crossprod(matrix(rnorm(25), nrow = 5))
#' ev <- am_eigh(adgeMatrix(S))
#' length(ev$values)
#'
#' @seealso \code{\link{am_rsvd}}
#' @export
am_eigh <- function(x) {
  am_eigen(.amatrix_promote_arg(x), symmetric = TRUE)
}

# ── Weighted am_crossprod helpers ─────────────────────────────────────────────

#' Weighted cross-product X'WX
#'
#' Computes \eqn{X^T \mathrm{diag}(w) X}, a \code{p x p} weighted
#' cross-product. A GPU-resident fast path is used when available.
#'
#' @param X Numeric matrix or \code{adgeMatrix} of shape \code{[n, p]}.
#' @param w Positive numeric vector of length \code{n}; observation
#'   weights.
#'
#' @return An \code{adgeMatrix} of shape \code{[p, p]}.
#'
#' @examples
#' X <- matrix(rnorm(20), nrow = 5)
#' w <- runif(5)
#' am_crossprod_weighted(X, w)
#'
#' @seealso \code{\link{am_tcrossprod_weighted}}, \code{\link{am_xty_weighted}}
#' @export
# X' diag(w) X  (p x p)
am_crossprod_weighted <- function(X, w) {
  X_arg <- .amatrix_model_dense_arg(X)
  w <- as.double(w)
  if (length(w) != nrow(X_arg)) {
    stop("length(w) must equal nrow(X)", call. = FALSE)
  }
  sqrt_w <- sqrt(w)

  # GPU-resident path: scale rows on GPU, then crossprod on GPU — no host round-trip
  backend_name <- .amatrix_live_resident_backend(X_arg)
  if (!is.null(backend_name)) {
    backend <- .amatrix_get_backend(backend_name)
    if (.amatrix_backend_supports_resident_op(backend, "broadcast_ewise", x = X_arg) &&
        .amatrix_backend_supports_resident_op(backend, "crossprod", x = X_arg)) {
      lhs <- .amatrix_prepare_resident_arg(X_arg, backend_name)
      if (!is.null(lhs)) {
        scaled_key <- .amatrix_next_resident_key(backend_name)
        out_key    <- .amatrix_next_resident_key(backend_name)
        release_scaled <- TRUE
        release_output <- TRUE
        on.exit({
          if (release_output) {
            .amatrix_release_resident_key(backend_name, out_key)
          }
          if (release_scaled) {
            .amatrix_release_resident_key(backend_name, scaled_key)
          }
          .amatrix_cleanup_temp_resident_safe(list(lhs), backend_name)
        }, add = TRUE)
        result <- tryCatch({
          backend$broadcast_ewise_resident(lhs$key, sqrt_w, 1L, "*", scaled_key)
          val <- backend$crossprod_resident(scaled_key, NULL, out_key)
          val
        }, error = function(e) NULL)
        if (!is.null(result)) {
          value <- .amatrix_rewrap_value(X_arg, result)
          bound <- .amatrix_bind_resident(value, backend_name, out_key)
          release_output <- FALSE
          return(bound)
        }
      }
    }
  }

  # CPU fallback
  x_host <- as.matrix(amatrix_materialize_host(X_arg))
  x_scaled <- x_host * sqrt_w
  am_crossprod(.amatrix_rewrap_like(X_arg, x_scaled))
}

#' Row-weighted outer cross-product
#'
#' Computes \eqn{W^{1/2} X X^T W^{1/2}} with
#' \eqn{W = \mathrm{diag}(w)}: the \code{n x n} outer cross-product of
#' the row-weighted matrix \eqn{W^{1/2} X}. A GPU-resident fast path is used when available.
#'
#' @param X Numeric matrix or \code{adgeMatrix} of shape \code{[n, p]}.
#' @param w Positive numeric vector of length \code{n}; observation
#'   weights.
#'
#' @return An \code{adgeMatrix} of shape \code{[n, n]}.
#'
#' @examples
#' X <- matrix(rnorm(20), nrow = 5)
#' w <- runif(5)
#' am_tcrossprod_weighted(X, w)
#'
#' @seealso \code{\link{am_crossprod_weighted}}, \code{\link{am_xty_weighted}}
#' @export
# X diag(w) X'  (n x n)
am_tcrossprod_weighted <- function(X, w) {
  X_arg <- .amatrix_model_dense_arg(X)
  w <- as.double(w)
  if (length(w) != nrow(X_arg)) {
    stop("length(w) must equal nrow(X)", call. = FALSE)
  }
  sqrt_w <- sqrt(w)

  # GPU-resident path: scale rows on GPU, then tcrossprod on GPU
  backend_name <- .amatrix_live_resident_backend(X_arg)
  if (!is.null(backend_name)) {
    backend <- .amatrix_get_backend(backend_name)
    if (.amatrix_backend_supports_resident_op(backend, "broadcast_ewise", x = X_arg) &&
        .amatrix_backend_supports_resident_op(backend, "tcrossprod", x = X_arg)) {
      lhs <- .amatrix_prepare_resident_arg(X_arg, backend_name)
      if (!is.null(lhs)) {
        scaled_key <- .amatrix_next_resident_key(backend_name)
        out_key    <- .amatrix_next_resident_key(backend_name)
        release_scaled <- TRUE
        release_output <- TRUE
        on.exit({
          if (release_output) {
            .amatrix_release_resident_key(backend_name, out_key)
          }
          if (release_scaled) {
            .amatrix_release_resident_key(backend_name, scaled_key)
          }
          .amatrix_cleanup_temp_resident_safe(list(lhs), backend_name)
        }, add = TRUE)
        result <- tryCatch({
          backend$broadcast_ewise_resident(lhs$key, sqrt_w, 1L, "*", scaled_key)
          val <- backend$tcrossprod_resident(scaled_key, NULL, out_key)
          val
        }, error = function(e) NULL)
        if (!is.null(result)) {
          value <- .amatrix_rewrap_value(X_arg, result)
          bound <- .amatrix_bind_resident(value, backend_name, out_key)
          release_output <- FALSE
          return(bound)
        }
      }
    }
  }

  # CPU fallback
  x_host <- as.matrix(amatrix_materialize_host(X_arg))
  x_scaled <- x_host * sqrt_w
  am_tcrossprod(.amatrix_rewrap_like(X_arg, x_scaled))
}

#' Weighted cross-product X'Wy
#'
#' Computes \eqn{X^T \mathrm{diag}(w) y}, a \code{p x k} weighted
#' cross-product between \code{X} and response matrix \code{y}. A
#' GPU-resident fast path is used when available.
#'
#' @param X Numeric matrix or \code{adgeMatrix} of shape \code{[n, p]}.
#' @param w Positive numeric vector of length \code{n}; observation
#'   weights.
#' @param y Numeric vector or matrix of shape \code{[n, k]};
#'   response(s).
#'
#' @return An \code{adgeMatrix} of shape \code{[p, k]}.
#'
#' @examples
#' X <- matrix(rnorm(20), nrow = 5)
#' w <- runif(5)
#' y <- rnorm(5)
#' am_xty_weighted(X, w, y)
#'
#' @seealso \code{\link{am_crossprod_weighted}},
#'   \code{\link{am_tcrossprod_weighted}}
#' @export
# X' diag(w) y  (p x k)
am_xty_weighted <- function(X, w, y) {
  X_arg <- .amatrix_model_dense_arg(X)
  w <- as.double(w)
  if (length(w) != nrow(X_arg)) {
    stop("length(w) must equal nrow(X)", call. = FALSE)
  }
  sqrt_w <- sqrt(w)
  y_mat <- if (is.vector(y)) matrix(y, ncol = 1L) else as.matrix(y)
  if (nrow(y_mat) != nrow(X_arg)) {
    stop("nrow(y) must equal nrow(X)", call. = FALSE)
  }

  # GPU-resident path: scale X and y on GPU, crossprod on GPU
  backend_name <- .amatrix_live_resident_backend(X_arg)
  if (!is.null(backend_name)) {
    backend <- .amatrix_get_backend(backend_name)
    if (.amatrix_backend_supports_resident_op(backend, "broadcast_ewise", x = X_arg) &&
        .amatrix_backend_supports_resident_op(backend, "crossprod", x = X_arg)) {
      lhs <- .amatrix_prepare_resident_arg(X_arg, backend_name)
      if (!is.null(lhs)) {
        x_scaled_key <- .amatrix_next_resident_key(backend_name)
        y_scaled_key <- .amatrix_next_resident_key(backend_name)
        out_key      <- .amatrix_next_resident_key(backend_name)
        y_key <- .amatrix_next_resident_key(backend_name)
        release_x_scaled <- TRUE
        release_y_scaled <- TRUE
        release_y <- TRUE
        release_output <- TRUE
        on.exit({
          if (release_output) {
            .amatrix_release_resident_key(backend_name, out_key)
          }
          if (release_y_scaled) {
            .amatrix_release_resident_key(backend_name, y_scaled_key)
          }
          if (release_x_scaled) {
            .amatrix_release_resident_key(backend_name, x_scaled_key)
          }
          if (release_y) {
            .amatrix_release_resident_key(backend_name, y_key)
          }
          .amatrix_cleanup_temp_resident_safe(list(lhs), backend_name)
        }, add = TRUE)
        result <- tryCatch({
          # Scale X rows by sqrt(w) on GPU
          backend$broadcast_ewise_resident(lhs$key, sqrt_w, 1L, "*", x_scaled_key)
          # Upload and scale y rows by sqrt(w) on GPU
          backend$resident_store(y_key, y_mat)
          backend$broadcast_ewise_resident(y_key, sqrt_w, 1L, "*", y_scaled_key)
          # crossprod(X_scaled, y_scaled)
          val <- backend$crossprod_resident(x_scaled_key, y_scaled_key, out_key)
          val
        }, error = function(e) NULL)
        if (!is.null(result)) {
          value <- .amatrix_rewrap_value(X_arg, result)
          bound <- .amatrix_bind_resident(value, backend_name, out_key)
          release_output <- FALSE
          return(bound)
        }
      }
    }
  }

  # CPU fallback
  x_host <- as.matrix(amatrix_materialize_host(X_arg))
  x_scaled <- x_host * sqrt_w
  y_scaled <- y_mat * sqrt_w
  am_crossprod(
    .amatrix_rewrap_like(X_arg, x_scaled),
    .amatrix_rewrap_like(X_arg, y_scaled)
  )
}

# ── Diagonal scaling ──────────────────────────────────────────────────────────

#' Row and column diagonal scaling
#'
#' Scale each row or column of a matrix by a numeric vector, equivalent to
#' left- or right-multiplying by a diagonal matrix.
#' \code{am_rowscale} computes \code{diag(d) \%*\% X} (row \eqn{i} scaled by
#' \code{d[i]}); \code{am_colscale} computes \code{X \%*\% diag(d)} (column
#' \eqn{j} scaled by \code{d[j]}).
#'
#' @param X A matrix or \code{aMatrix} object.
#' @param d Numeric vector of scale factors. Length must equal \code{nrow(X)}
#'   for \code{am_rowscale} and \code{ncol(X)} for \code{am_colscale}.
#'
#' @return A matrix or \code{aMatrix} of the same dimensions as \code{X}.
#'
#' @examples
#' m <- matrix(1:6, 2, 3)
#' am_rowscale(m, c(2, 0.5))
#' am_colscale(m, c(1, 2, 3))
#'
#' @export
am_rowscale <- function(X, d) {
  if (inherits(X, c("adgCMatrix", "sparseMatrix"))) {
    if (length(d) != nrow(X)) stop("length(d) must equal nrow(X)", call. = FALSE)
    return(.amatrix_sparse_scale(X, 1, d, "*"))   # stays sparse
  }
  X_arg <- .amatrix_model_dense_arg(X)
  d <- as.double(d)
  if (length(d) != nrow(X_arg)) {
    stop("length(d) must equal nrow(X)", call. = FALSE)
  }
  x_host <- as.matrix(amatrix_materialize_host(X_arg))
  .amatrix_rewrap_like(X_arg, x_host * d)
}

#' @rdname am_rowscale
#' @export
am_colscale <- function(X, d) {
  if (inherits(X, c("adgCMatrix", "sparseMatrix"))) {
    if (length(d) != ncol(X)) stop("length(d) must equal ncol(X)", call. = FALSE)
    return(.amatrix_sparse_scale(X, 2, d, "*"))   # stays sparse
  }
  X_arg <- .amatrix_model_dense_arg(X)
  d <- as.double(d)
  if (length(d) != ncol(X_arg)) {
    stop("length(d) must equal ncol(X)", call. = FALSE)
  }
  x_host <- as.matrix(amatrix_materialize_host(X_arg))
  .amatrix_rewrap_like(X_arg, t(t(x_host) * d))
}

am_diag <- function(x, nrow, ncol, names = TRUE) {
  # Extract mode: diag(matrix_x) → numeric vector of diagonal elements
  extract_mode <- (missing(nrow) && missing(ncol)) &&
                  (is.matrix(x) || inherits(x, "aMatrix"))
  if (extract_mode) {
    # Read the diagonal straight from the host representation; as.matrix()
    # would copy a dense matrix and fully densify a sparse one.
    x_host <- amatrix_materialize_host(x)
    if (inherits(x_host, "Matrix")) {
      return(Matrix::diag(x_host, names = names))
    }
    return(base::diag(x_host, names = names))
  }
  # Create mode: diag(d) → diagonal adgeMatrix; nrow/ncol set the size
  nrow <- if (missing(nrow)) NULL else nrow
  ncol <- if (missing(ncol)) NULL else ncol
  amatrix_dispatch_op(
    x = x,
    op = "diag",
    method = "diag",
    args = list(nrow = nrow, ncol = ncol, names = names),
    fallback = function() {
      args <- list(as.matrix(amatrix_materialize_host(x)), names = names)
      if (!is.null(nrow)) args$nrow <- nrow
      if (!is.null(ncol)) args$ncol <- ncol
      do.call(base::diag, args)
    }
  )
}

# ── Fused crossprod + diagonal add ───────────────────────────────────────────

#' Cross-product plus diagonal perturbation
#'
#' Computes \eqn{X^T X + \lambda I} (scalar \code{lambda}) or
#' \eqn{X^T X + \mathrm{diag}(\lambda)} (vector \code{lambda}) in a
#' single fused call.
#'
#' @param X Numeric matrix or \code{adgeMatrix} of shape \code{[n, p]}.
#' @param lambda Scalar or numeric vector of length \code{p}; diagonal
#'   perturbation to add to the cross-product.
#'
#' @return An \code{adgeMatrix} of shape \code{[p, p]}: the perturbed
#'   cross-product.
#'
#' @examples
#' X <- matrix(rnorm(20), nrow = 5)
#' am_crossprod_add_diag(X, lambda = 0.1)
#'
#' @seealso \code{\link{am_crossprod_weighted}}
#' @export
# X'X + lambda*I  or  X'X + diag(d)
am_crossprod_add_diag <- function(X, lambda) {
  X_arg  <- .amatrix_model_dense_arg(X)
  xtx    <- am_crossprod(X_arg)
  p      <- ncol(X_arg)
  xtx_m  <- as.matrix(amatrix_materialize_host(xtx))
  if (length(lambda) == 1L) {
    diag(xtx_m) <- diag(xtx_m) + as.double(lambda)
  } else {
    if (length(lambda) != p)
      stop("lambda must be a scalar or length ncol(X)", call. = FALSE)
    diag(xtx_m) <- diag(xtx_m) + as.double(lambda)
  }
  .amatrix_rewrap_like(X_arg, xtx_m)
}

# ── Matrix functions (via symmetric eigendecomposition) ───────────────────────

.mat_fun <- function(X, f, check_positive = TRUE) {
  X_arg <- .amatrix_model_dense_arg(X)
  res   <- am_eigh(X_arg)
  lam   <- res$values
  if (isTRUE(check_positive) && any(lam <= 0))
    warning("matrix has non-positive eigenvalues; result may be complex or NaN",
            call. = FALSE)
  Q     <- as.matrix(res$vectors)
  new_lam <- f(lam)
  # Scale the columns of Q by the transformed eigenvalues. (diag(new_lam) is
  # avoided: for a 1 x 1 input diag(<scalar>) builds an identity of that size.)
  result <- (Q * rep(new_lam, each = nrow(Q))) %*% t(Q)
  dimnames(result) <- dimnames(X_arg)
  .amatrix_rewrap_like(X_arg, result)
}

#' Matrix functions via symmetric eigendecomposition
#'
#' Apply an elementwise function to the eigenvalues of a symmetric
#' positive definite matrix and reconstruct the result:
#' \eqn{f(X) = Q \, \mathrm{diag}(f(\lambda)) \, Q^T}.
#'
#' @param X Symmetric positive definite numeric matrix or
#'   \code{adgeMatrix} of shape \code{[p, p]}.
#' @param p Numeric scalar exponent (used by \code{am_powm} only).
#'
#' @return An \code{adgeMatrix} of shape \code{[p, p]}: the matrix
#'   function applied to \code{X}.
#'
#' @examples
#' S <- crossprod(matrix(rnorm(16), 4)) + diag(4)
#' am_sqrtm(S)
#' am_logm(S)
#' am_powm(S, -1)
#'
#' @name am_matfun
#' @rdname am_matfun
#' @export
am_sqrtm <- function(X) .mat_fun(X, sqrt)

#' @rdname am_matfun
#' @export
am_powm  <- function(X, p) .mat_fun(X, function(lam) lam^p)

#' @rdname am_matfun
#' @export
am_logm  <- function(X) .mat_fun(X, log)

# ── Stochastic trace estimator (Hutchinson) ───────────────────────────────────

#' Stochastic trace estimator (Hutchinson)
#'
#' Estimates \eqn{\mathrm{tr}(A)} or \eqn{\mathrm{tr}(A^{-1})} using
#' Hutchinson's method with \code{k} Rademacher probe vectors. Supply
#' \code{solve_fn} to estimate the trace of an inverse without forming it
#' explicitly.
#'
#' @param A Square matrix or \code{aMatrix}; required when \code{solve_fn}
#'   is \code{NULL}.
#' @param k Integer number of Rademacher probe vectors. Default \code{30L}.
#' @param seed Optional integer random seed for reproducibility.
#' @param solve_fn Optional function \code{function(V)} that returns
#'   \code{A^{-1} \%*\% V}; use this to estimate \eqn{\mathrm{tr}(A^{-1})}
#'   without materialising the inverse.
#' @param n Integer dimension of the matrix; required when \code{solve_fn}
#'   is supplied.
#'
#' @return A single numeric scalar estimate of the trace.
#'
#' @examples
#' A <- crossprod(matrix(rnorm(25), 5, 5)) + 5 * diag(5)
#' am_trace_estim(A, k = 50L, seed = 1L)
#'
#' @export
am_trace_estim <- function(A = NULL, k = 30L, seed = NULL,
                        solve_fn = NULL, n = NULL) {
  if (!is.null(seed)) {
    # Seed only this call: restore the caller's RNG state on exit.
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = globalenv())
      } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      }
    }, add = TRUE)
    set.seed(seed)
  }

  if (!is.null(solve_fn)) {
    if (is.null(n)) stop("n must be supplied when using solve_fn", call. = FALSE)
    probes <- matrix(sample(c(-1, 1), n * k, replace = TRUE), n, k)
    sols   <- solve_fn(probes)
    return(mean(colSums(probes * as.matrix(sols))))
  }

  A_host <- as.matrix(amatrix_materialize_host(A))
  n      <- nrow(A_host)
  probes <- matrix(sample(c(-1, 1), n * k, replace = TRUE), n, k)
  mean(colSums(probes * (A_host %*% probes)))
}

# ── Row / column means ────────────────────────────────────────────────────────

#' Row and column means
#'
#' Compute row or column means of a matrix or \code{aMatrix}, dispatching
#' to an accelerated backend when one is available.
#'
#' @param x A matrix or \code{aMatrix} object.
#' @param na.rm Logical; if \code{TRUE}, \code{NA} values are excluded before
#'   averaging. Default \code{FALSE}.
#'
#' @return A numeric vector of length \code{nrow(x)} (\code{am_rowmeans}) or
#'   \code{ncol(x)} (\code{am_colmeans}).
#'
#' @examples
#' m <- matrix(1:12, 3, 4)
#' am_rowmeans(m)
#' am_colmeans(m)
#'
#' @export
am_rowmeans <- function(x, na.rm = FALSE) {
  if (inherits(x, "adgCMatrix")) {
    return(Matrix::rowMeans(amatrix_materialize_host(x), na.rm = na.rm))
  }
  if (inherits(x, "adgeMatrix")) {
    if (isTRUE(na.rm)) {
      return(base::rowMeans(as.matrix(amatrix_materialize_host(x)), na.rm = TRUE))
    }
    return(am_rowsums(x, na.rm = na.rm) / ncol(x))
  }
  x_host <- as.matrix(amatrix_materialize_host(x))
  base::rowMeans(x_host, na.rm = na.rm)
}

#' @rdname am_rowmeans
#' @export
am_colmeans <- function(x, na.rm = FALSE) {
  if (inherits(x, "adgCMatrix")) {
    return(Matrix::colMeans(amatrix_materialize_host(x), na.rm = na.rm))
  }
  if (inherits(x, "adgeMatrix")) {
    if (isTRUE(na.rm)) {
      return(base::colMeans(as.matrix(amatrix_materialize_host(x)), na.rm = TRUE))
    }
    return(am_colsums(x, na.rm = na.rm) / nrow(x))
  }
  x_host <- as.matrix(amatrix_materialize_host(x))
  base::colMeans(x_host, na.rm = na.rm)
}

# ── Matrix trace ──────────────────────────────────────────────────────────────

#' Matrix trace
#'
#' Returns the trace (sum of diagonal elements) of a square matrix or
#' \code{aMatrix}.
#'
#' @param x A square matrix, sparse \code{sparseMatrix}, or \code{aMatrix}.
#'
#' @return A single numeric scalar equal to the sum of diagonal elements.
#'
#' @examples
#' am_trace(diag(1:4))
#'
#' @export
am_trace <- function(x) {
  host <- amatrix_materialize_host(x)
  if (inherits(host, "sparseMatrix")) return(sum(Matrix::diag(host)))
  sum(base::diag(as.matrix(host)))
}

# ── Symmetry enforcement ──────────────────────────────────────────────────────

#' Symmetrise a matrix
#'
#' Returns \code{(x + t(x)) / 2}, enforcing exact symmetry. Handles both
#' dense \code{aMatrix} and sparse \code{adgCMatrix} inputs.
#'
#' @param x A square matrix or \code{aMatrix} object.
#'
#' @return A symmetric matrix or \code{aMatrix} of the same class and
#'   dimensions as \code{x}.
#'
#' @examples
#' m <- matrix(c(1, 2, 3, 4), 2, 2)
#' am_sym(m)
#'
#' @export
am_sym <- function(x) {
  if (inherits(x, "adgCMatrix")) {
    host <- amatrix_materialize_host(x)
    result <- (host + Matrix::t(host)) / 2
    return(new_adgCMatrix(result, preferred_backend = x@preferred_backend,
                          precision = x@precision, policy = x@policy))
  }
  X_arg <- .amatrix_model_dense_arg(x)
  x_host <- as.matrix(amatrix_materialize_host(X_arg))
  .amatrix_rewrap_like(X_arg, (x_host + t(x_host)) / 2)
}

# ── Inner product ─────────────────────────────────────────────────────────────

#' Inner product of two vectors or matrices
#'
#' Computes the element-wise inner product \code{sum(x * y)}, equivalent to
#' \code{as.numeric(t(x) \%*\% y)} for vectors.
#'
#' @param x A numeric vector, matrix, or \code{aMatrix}.
#' @param y A numeric vector, matrix, or \code{aMatrix} conformable with
#'   \code{x}.
#'
#' @return A single numeric scalar.
#'
#' @examples
#' am_dot(1:4, 4:1)
#'
#' @export
am_dot <- function(x, y) {
  x_host <- .amatrix_host_arg(x)
  y_host <- .amatrix_host_arg(y)
  sum(x_host * y_host)
}

# ── segment_sum / segment_mean (amatrix-ylo) ──────────────────────────────
# First-class grouped-reduction primitives.  GPU path stores result as a
# resident adgeMatrix (no data downloaded). CPU path uses base::rowsum.

.am_segment_sum_cpu <- function(X_mat, labels, K) {
  # NA labels belong to no segment (the sparse C kernel skips them too).
  if (anyNA(labels)) {
    keep <- !is.na(labels)
    X_mat <- X_mat[keep, , drop = FALSE]
    labels <- labels[keep]
  }
  if (length(labels) == 0L) {
    return(matrix(0, K, ncol(X_mat)))
  }
  sums_raw <- rowsum(X_mat, labels, reorder = FALSE)
  out <- matrix(0, K, ncol(X_mat))
  idx <- as.integer(rownames(sums_raw))
  valid <- idx >= 1L & idx <= K
  if (any(valid)) out[idx[valid], ] <- sums_raw[valid, , drop = FALSE]
  out
}

.am_segment_mean_cpu <- function(X_mat, labels, K) {
  if (anyNA(labels)) {
    keep <- !is.na(labels)
    X_mat <- X_mat[keep, , drop = FALSE]
    labels <- labels[keep]
  }
  if (length(labels) == 0L) {
    return(matrix(NA_real_, K, ncol(X_mat)))
  }
  sums_raw <- rowsum(X_mat, labels, reorder = FALSE)
  idx    <- as.integer(rownames(sums_raw))
  counts <- tabulate(labels, nbins = K)
  out    <- matrix(NA_real_, K, ncol(X_mat))
  valid  <- idx >= 1L & idx <= K
  if (any(valid)) {
    k  <- idx[valid]
    nz <- counts[k] > 0L
    if (any(nz))
      out[k[nz], ] <- sums_raw[valid, , drop = FALSE][nz, , drop = FALSE] /
                      counts[k[nz]]
  }
  out
}

.am_try_resident_segment_op <- function(x, labels, K, backend_name, op_name) {
  backend <- .amatrix_get_backend(backend_name)
  fn_name <- paste0(op_name, "_resident")
  if (!is.function(backend[[fn_name]])) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(x, backend_name)
  if (is.null(lhs)) return(NULL)
  out_key <- .amatrix_next_resident_key(backend_name)
  value   <- backend[[fn_name]](lhs$key, labels, K, out_key)
  .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  # Bridge returns a plain R matrix; fall back to resident_materialize if not
  if (!is.matrix(value)) value <- backend$resident_materialize(out_key)
  list(value = value, key = out_key)
}

.am_segment_resident_wrap <- function(x, resident, choice_name) {
  wrapped <- new_adgeMatrix(resident$value,
                            preferred_backend = x@preferred_backend,
                            precision        = x@precision,
                            policy           = x@policy)
  .amatrix_bind_resident(wrapped, choice_name, resident$key)
}

#' Segment sum by group labels
#'
#' Sum rows of \code{x} grouped by integer \code{labels}, dispatching to
#' GPU when available.
#'
#' @param x A numeric matrix or \code{adgeMatrix}.
#' @param labels Integer vector of group labels (1-based).
#' @param K Number of groups.
#' @return A \code{K}-by-\code{ncol(x)} matrix of group sums.
#'
#' @examples
#' x <- matrix(1:12, 6, 2)
#' labels <- c(1L, 2L, 1L, 3L, 2L, 1L)
#' am_segment_sum(x, labels, K = 3L)
#'
#' @seealso \code{\link{am_segment_mean}}, \code{\link{am_scatter_mean}}
#' @export
am_segment_sum <- function(x, labels, K) {
  labels <- as.integer(labels)
  K      <- as.integer(K)
  if (inherits(x, "adgCMatrix")) {
    sums_raw <- rowsum.adgCMatrix(x, group = labels, reorder = FALSE)
    out <- matrix(0, K, ncol(x))
    idx <- as.integer(rownames(sums_raw))
    valid <- idx >= 1L & idx <= K
    if (any(valid)) out[idx[valid], ] <- sums_raw[valid, , drop = FALSE]
    return(out)
  }
  if (!inherits(x, "adgeMatrix"))
    return(.am_segment_sum_cpu(as.matrix(x), labels, K))
  # Size gate: GPU kernel launch overhead exceeds compute for small matrices.
  # Fall back to CPU when n*p < amatrix.segment_min_size (default 500k).
  seg_min <- as.numeric(getOption("amatrix.segment_min_size", 500000L))
  if (nrow(x) * ncol(x) < seg_min)
    return(.am_segment_sum_cpu(as.matrix(amatrix_materialize_host(x)), labels, K))
  choice   <- .amatrix_backend_for(x, "segment_sum")
  resident <- .am_try_resident_segment_op(x, labels, K, choice$name, "segment_sum")
  if (is.null(resident))
    return(.am_segment_sum_cpu(as.matrix(amatrix_materialize_host(x)), labels, K))
  .am_segment_resident_wrap(x, resident, choice$name)
}

#' Segment mean by group labels
#'
#' Compute the mean of rows of \code{x} grouped by integer \code{labels},
#' dispatching to GPU when available.
#'
#' @param x A numeric matrix or \code{adgeMatrix}.
#' @param labels Integer vector of group labels (1-based).
#' @param K Number of groups.
#' @return A \code{K}-by-\code{ncol(x)} matrix of group means.
#'
#' @examples
#' x <- matrix(1:12, 6, 2)
#' labels <- c(1L, 2L, 1L, 3L, 2L, 1L)
#' am_segment_mean(x, labels, K = 3L)
#'
#' @seealso \code{\link{am_segment_sum}}, \code{\link{am_scatter_mean}}
#' @export
am_segment_mean <- function(x, labels, K) {
  labels <- as.integer(labels)
  K      <- as.integer(K)
  if (inherits(x, "adgCMatrix")) {
    sums_raw <- rowsum.adgCMatrix(x, group = labels, reorder = FALSE)
    idx    <- as.integer(rownames(sums_raw))
    counts <- tabulate(labels, nbins = K)
    out    <- matrix(NA_real_, K, ncol(x))
    valid  <- idx >= 1L & idx <= K
    if (any(valid)) {
      k  <- idx[valid]
      nz <- counts[k] > 0L
      if (any(nz))
        out[k[nz], ] <- sums_raw[valid, , drop = FALSE][nz, , drop = FALSE] /
          counts[k[nz]]
    }
    return(out)
  }
  if (!inherits(x, "adgeMatrix"))
    return(.am_segment_mean_cpu(as.matrix(x), labels, K))
  # Size gate: same as segment_sum — CPU is faster for small n*p.
  seg_min <- as.numeric(getOption("amatrix.segment_min_size", 500000L))
  if (nrow(x) * ncol(x) < seg_min)
    return(.am_segment_mean_cpu(as.matrix(amatrix_materialize_host(x)), labels, K))
  choice   <- .amatrix_backend_for(x, "segment_mean")
  resident <- .am_try_resident_segment_op(x, labels, K, choice$name, "segment_mean")
  if (is.null(resident))
    return(.am_segment_mean_cpu(as.matrix(amatrix_materialize_host(x)), labels, K))
  .am_segment_resident_wrap(x, resident, choice$name)
}

# ── rowsum.adgeMatrix ──────────────────────────────────────────────────────
# Intercept base::rowsum(X, group) for GPU matrices and route to segment_sum.
# group can be integer, numeric, character, or factor — we map to 1..K labels.

#' @export
rowsum.adgeMatrix <- function(x, group, reorder = TRUE, na.rm = FALSE, ...) {
  # na.rm = TRUE and missing group values follow base::rowsum exactly
  # (including its NA row and warning), so defer to it on the host matrix.
  if (isTRUE(na.rm) || anyNA(group)) {
    host <- as.matrix(amatrix_materialize_host(x))
    dimnames(host) <- dimnames(x)
    return(base::rowsum(host, group, reorder = reorder, na.rm = na.rm, ...))
  }
  if (is.factor(group)) {
    # reorder=TRUE: honour factor's own level ordering (user-defined, preserves
    #   all levels including unoccupied ones — e.g. empty k-means clusters).
    # reorder=FALSE: first-occurrence order, matching base::rowsum semantics.
    lvls <- if (reorder) levels(group) else unique(as.character(group))
    grp  <- factor(as.character(group), levels = lvls)
  } else {
    g_chr <- as.character(group)
    if (reorder) {
      u     <- unique(g_chr)
      num   <- suppressWarnings(as.numeric(u))
      lvls  <- if (!anyNA(num)) u[order(num)] else sort(u)
    } else {
      lvls  <- unique(g_chr)
    }
    grp   <- factor(g_chr, levels = lvls)
  }
  labels <- as.integer(grp)
  K      <- nlevels(grp)
  result <- am_segment_sum(x, labels, K)
  if (is.matrix(result)) dimnames(result) <- list(lvls, colnames(x))
  result
}

# ── rowsum.adgCMatrix ──────────────────────────────────────────────────────────
# Group row sums for sparse matrices. Result is always dense (K × p).
# Avoids densifying the full matrix — uses Matrix::rowsum on the dgCMatrix host.

#' @export
rowsum.adgCMatrix <- function(x, group, reorder = TRUE, na.rm = FALSE, ...) {
  host <- amatrix_materialize_host(x)   # dgCMatrix

  if (isTRUE(na.rm)) {
    return(base::rowsum(as.matrix(host), group, reorder = reorder, na.rm = TRUE, ...))
  }

  # Match base::rowsum behavior: unique groups in first-occurrence order
  ugrp <- unique(group)
  if (anyNA(ugrp)) {
    # As base::rowsum: missing group values form their own (last) group.
    warning("missing values for 'group'", call. = FALSE)
  }
  if (isTRUE(reorder)) ugrp <- sort(ugrp, na.last = TRUE)
  K <- length(ugrp)
  # Map each row to 1-based group index
  labels <- match(group, ugrp)

  result <- .Call("am_sparse_segment_sum_c", PACKAGE = "amatrix",
                  as.double(host@x), as.integer(host@p), as.integer(host@i),
                  as.integer(host@Dim), as.integer(labels), as.integer(K))
  rownames(result) <- as.character(ugrp)
  colnames(result) <- colnames(host)

  result
}

# ── sweep() for amatrix objects ───────────────────────────────────────────────
# base::sweep() is an ordinary function, not a generic, so S3-style methods
# (the former sweep.adgeMatrix / sweep.adgCMatrix) were never called. sweep is
# declared as a standard S4 generic instead, with base::sweep as its default,
# so idiomatic sweep(x, MARGIN, STATS, FUN) reaches the backend kernels for
# dense input and the O(nnz) path for sparse input.

# Scale the rows (MARGIN = 1) or columns (MARGIN = 2) of a sparse matrix by
# multiplying or dividing its stored values: O(nnz) and sparsity preserving.
# Returns NULL when the request is not of that form.
.amatrix_sparse_scale <- function(x, MARGIN, STATS, FUN) {
  if (!is.character(FUN) || length(FUN) != 1L || !FUN %in% c("*", "/")) return(NULL)
  if (!is.numeric(MARGIN) || length(MARGIN) != 1L || !MARGIN %in% c(1, 2)) return(NULL)
  if (!is.numeric(STATS) || !is.null(dim(STATS))) return(NULL)
  host <- if (inherits(x, "adgCMatrix")) {
    amatrix_materialize_host(x)
  } else {
    methods::as(methods::as(methods::as(x, "dMatrix"), "generalMatrix"), "CsparseMatrix")
  }
  if (length(STATS) != dim(host)[MARGIN]) return(NULL)
  STATS <- as.double(STATS)
  factor <- if (MARGIN == 1) {
    STATS[host@i + 1L]                 # @i holds 0-based row indices
  } else {
    rep.int(STATS, diff(host@p))       # one entry per stored value, by column
  }
  host@x <- if (identical(FUN, "*")) host@x * factor else host@x / factor
  host@factors <- list()
  if (inherits(x, "adgCMatrix")) {
    new_adgCMatrix(host, preferred_backend = x@preferred_backend,
                   precision = x@precision, policy = x@policy)
  } else {
    new_adgCMatrix(host)
  }
}

#' @noRd
setGeneric("sweep")

#' @noRd
setMethod("sweep", "adgeMatrix", function(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...) {
  if (...length() > 0L) {
    return(.amatrix_rewrap_value(
      x, base::sweep(as.matrix(amatrix_materialize_host(x)), MARGIN, STATS, FUN, check.margin, ...)
    ))
  }
  am_sweep(x, MARGIN, STATS, FUN)
})

#' @noRd
setMethod("sweep", "adgCMatrix", function(x, MARGIN, STATS, FUN = "-", check.margin = TRUE, ...) {
  if (...length() == 0L) {
    scaled <- .amatrix_sparse_scale(x, MARGIN, STATS, FUN)
    if (!is.null(scaled)) return(scaled)
  }
  # Other operations (e.g. "+" / "-") fill in the zeros, so the result is dense.
  .amatrix_rewrap_value(
    x, base::sweep(as.matrix(amatrix_materialize_host(x)), MARGIN, STATS, FUN, check.margin, ...)
  )
})

# ── addmm (amatrix-uaj) ─────────────────────────────────────────────────
# alpha*(A%*%B) + beta*C  — BLAS-3 fused scaled matmul with optional bias.
# A: n×p adgeMatrix (resident if GPU); B: p×k R matrix; C: n×k R matrix or NULL.
# GPU path uses mlx_addmm directly; CPU path uses plain R arithmetic.

.am_addmm_cpu <- function(A_mat, B_mat, C_mat, alpha, beta) {
  result <- alpha * (A_mat %*% B_mat)
  if (!is.null(C_mat) && beta != 0) result <- result + beta * C_mat
  result
}

.am_try_addmm_gpu <- function(A, B_mat, C_mat, alpha, beta, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!is.function(backend$addmm_resident)) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(A, backend_name)
  if (is.null(lhs)) return(NULL)
  out_key <- .amatrix_next_resident_key(backend_name)
  value   <- backend$addmm_resident(lhs$key, B_mat, C_mat, alpha, beta, out_key)
  .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  if (!is.matrix(value)) value <- backend$resident_materialize(out_key)
  list(value = value, key = out_key)
}

#' Scaled matrix multiply with optional bias: alpha*(A\%*\%B) + beta*C
#'
#' @param A  n×p \code{adgeMatrix} or plain matrix.
#' @param B  p×k numeric matrix.
#' @param C  n×k numeric matrix or \code{NULL} (treated as zeros).
#' @param alpha Scalar multiplier for \code{A\%*\%B} (default 1).
#' @param beta  Scalar multiplier for \code{C} (default 1).
#' @return \code{adgeMatrix} if A is resident, otherwise plain matrix.
#'
#' @examples
#' A <- adgeMatrix(matrix(1:6, 2, 3))
#' B <- matrix(1:12, 3, 4)
#' C <- matrix(1, 2, 4)
#' # 2 * (A %*% B) + 0.5 * C
#' am_addmm(A, B, C, alpha = 2, beta = 0.5)
#'
#' @export
am_addmm <- function(A, B, C = NULL, alpha = 1.0, beta = 1.0) {
  B_mat <- as.matrix(B); storage.mode(B_mat) <- "double"
  C_mat <- if (!is.null(C)) { m <- as.matrix(C); storage.mode(m) <- "double"; m } else NULL

  if (!inherits(A, "adgeMatrix")) {
    A_mat <- as.matrix(A); storage.mode(A_mat) <- "double"
    return(.am_addmm_cpu(A_mat, B_mat, C_mat, alpha, beta))
  }

  choice   <- .amatrix_backend_for(A, "addmm")
  resident <- .am_try_addmm_gpu(A, B_mat, C_mat, alpha, beta, choice$name)
  if (!is.null(resident)) {
    wrapped <- new_adgeMatrix(resident$value,
                              preferred_backend = A@preferred_backend,
                              precision = A@precision,
                              policy = A@policy)
    return(.amatrix_bind_resident(wrapped, choice$name, resident$key))
  }

  A_mat <- as.matrix(amatrix_materialize_host(A)); storage.mode(A_mat) <- "double"
  .am_addmm_cpu(A_mat, B_mat, C_mat, alpha, beta)
}

# ── pairwise_sqdist_argmin (amatrix-zas) ───────────────────────────────────
# Fused nearest-centroid assignment via the squared-distance identity:
#   D[i,k] = ||xi||^2 - 2*(X@Ct)[i,k] + ||ck||^2
# GPU path chains resident operations (no intermediate host round-trips):
#   1. cross  = matmul_resident(X, Ct)         [n×K]
#   2. neg2   = ewise_resident(cross, -2, "*") [n×K]
#   3. d1     = broadcast_ewise(neg2, x_norms, margin=1, "+") [add row norms]
#   4. d      = broadcast_ewise(d1,   c_norms, margin=2, "+") [add col norms]
#   5. labels = rowargmin_resident(d) + 1L     [0→1-indexed]
# CPU fallback: base R distance matrix + max.col.
.pairwise_sqdist_argmin_cpu <- function(X_mat, Ct_mat, x_norms, c_norms) {
  cross <- X_mat %*% Ct_mat                          # n×K
  D <- sweep(-2 * cross + x_norms, 2L, c_norms, "+")
  max.col(-D, ties.method = "first")
}

.pairwise_sqdist_argmin_gpu <- function(X, Ct_mat, x_norms, c_norms,
                                            backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  # All required resident ops must be present
  needed <- c("matmul_resident", "ewise_resident",
               "broadcast_ewise_resident", "rowargmin_resident")
  if (!all(vapply(needed, function(f) is.function(backend[[f]]), logical(1L))))
    return(NULL)

  lhs_X  <- .amatrix_prepare_resident_arg(X, backend_name)
  if (is.null(lhs_X)) return(NULL)

  # Upload Ct (p×K plain matrix) as temporary resident
  Ct_key  <- .amatrix_next_resident_key(backend_name)
  backend$resident_store(Ct_key, Ct_mat)
  temps   <- list(lhs_X, list(key = Ct_key, is_temp = TRUE))

  cross_key  <- .amatrix_next_resident_key(backend_name)
  neg2_key   <- .amatrix_next_resident_key(backend_name)
  d1_key     <- .amatrix_next_resident_key(backend_name)
  d_key      <- .amatrix_next_resident_key(backend_name)

  tryCatch({
    backend$matmul_resident(lhs_X$key, Ct_key, cross_key)
    backend$ewise_resident(cross_key, -2.0, "*", neg2_key)
    backend$broadcast_ewise_resident(neg2_key, as.double(x_norms), 1L, "+", d1_key)
    backend$broadcast_ewise_resident(d1_key, as.double(c_norms), 2L, "+", d_key)
    labels0 <- backend$rowargmin_resident(d_key)
    labels0 + 1L   # 0-indexed → 1-indexed
  }, error = function(e) NULL,
  finally = {
    .amatrix_cleanup_temp_resident(temps, backend_name)
    for (k in c(cross_key, neg2_key, d1_key, d_key))
      tryCatch(backend$resident_drop(k), error = function(e) invisible(NULL))
  })
}

#' Nearest-centroid assignment via fused squared-distance computation
#'
#' Computes \eqn{D[i,k] = \|x_i\|^2 - 2 x_i^\top c_k + \|c_k\|^2} and
#' returns \eqn{\arg\min_k D[i,k]} for each row \eqn{i}, 1-indexed.
#' GPU path avoids host round-trips by chaining resident operations.
#'
#' @param X       n×p \code{adgeMatrix} or plain matrix (query points).
#' @param Ct      p×K numeric matrix holding the centroids \strong{transposed}:
#'   \code{nrow(Ct)} is the feature dimension \eqn{p} and each \emph{column} is a
#'   centroid. Pass \code{t(centroids)} when your centroids are stored k×p with
#'   one centroid per row. A dimension mismatch (\code{ncol(X) != nrow(Ct)})
#'   errors, but a square \code{k == p} matrix passed untransposed cannot be
#'   detected and will silently assign to the wrong (column) centroids.
#' @param x_norms Optional n-vector of precomputed \eqn{\|x_i\|^2}. Computed if \code{NULL}.
#' @param c_norms Optional K-vector of precomputed \eqn{\|c_k\|^2}. Computed if \code{NULL}.
#' @return Integer vector of length n, 1-indexed nearest centroid per row.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(40), 10, 4)
#' centroids <- X[c(1, 5, 9), ]
#' # Centroids go in COLUMNS: pass t(centroids)
#' am_pairwise_sqdist_argmin(X, t(centroids))
#'
#' @export
am_pairwise_sqdist_argmin <- function(X, Ct, x_norms = NULL, c_norms = NULL) {
  Ct_mat  <- as.matrix(Ct);  storage.mode(Ct_mat) <- "double"
  if (ncol(X) != nrow(Ct_mat))
    stop(sprintf(
      paste0("ncol(X) (%d) must equal nrow(Ct) (%d): Ct must be p x K with ",
             "centroids in COLUMNS (did you forget t(centroids)?)"),
      ncol(X), nrow(Ct_mat)), call. = FALSE)
  if (is.null(x_norms)) {
    X_mat   <- if (inherits(X, "adgeMatrix")) as.matrix(amatrix_materialize_host(X))
               else as.matrix(X)
    x_norms <- rowSums(X_mat^2)
  } else {
    x_norms <- as.double(x_norms)
    X_mat   <- NULL
  }
  if (is.null(c_norms)) c_norms <- colSums(Ct_mat^2)
  else c_norms <- as.double(c_norms)
  # A wrong-length norm vector would otherwise be recycled silently.
  if (length(x_norms) != nrow(X))
    stop(sprintf("x_norms must have length nrow(X) = %d, not %d", nrow(X), length(x_norms)), call. = FALSE)
  if (length(c_norms) != ncol(Ct_mat))
    stop(sprintf("c_norms must have length ncol(Ct) = %d, not %d", ncol(Ct_mat), length(c_norms)), call. = FALSE)

  if (!inherits(X, "adgeMatrix")) {
    X_mat <- if (is.null(X_mat)) as.matrix(X) else X_mat
    storage.mode(X_mat) <- "double"
    return(.pairwise_sqdist_argmin_cpu(X_mat, Ct_mat, x_norms, c_norms))
  }

  choice <- .amatrix_backend_for(X, "matmul")
  result <- .pairwise_sqdist_argmin_gpu(X, Ct_mat, x_norms, c_norms, choice$name)
  if (!is.null(result)) return(result)

  # Fallback: materialize and use CPU
  X_mat <- if (is.null(X_mat)) as.matrix(amatrix_materialize_host(X)) else X_mat
  storage.mode(X_mat) <- "double"
  .pairwise_sqdist_argmin_cpu(X_mat, Ct_mat, x_norms, c_norms)
}

.am_scatter_mean_cpu <- function(X_mat, labels, K) {
  # One grouped pass (rowsum) instead of a scan of `labels` per cluster, which
  # was O(n * K). Labels that are NA or outside 1..K are ignored; empty
  # clusters are NA rows.
  centroids <- matrix(NA_real_, K, ncol(X_mat))
  valid <- !is.na(labels) & labels >= 1L & labels <= K
  if (!any(valid)) {
    return(centroids)
  }
  if (!all(valid)) {
    X_mat <- X_mat[valid, , drop = FALSE]
    labels <- labels[valid]
  }
  sums <- rowsum(X_mat, labels, reorder = FALSE)
  idx <- as.integer(rownames(sums))
  centroids[idx, ] <- sums / tabulate(labels, nbins = K)[idx]
  centroids
}

.amatrix_try_resident_scatter_mean <- function(x, labels, K, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!is.function(backend$scatter_mean_resident)) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(x, backend_name)
  if (is.null(lhs)) return(NULL)
  result <- backend$scatter_mean_resident(lhs$key, labels, K)
  .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  result
}

#' Scatter mean by group labels
#'
#' Compute the mean of rows of \code{x} grouped by integer \code{labels}.
#'
#' @param x A numeric matrix or \code{adgeMatrix}.
#' @param labels Integer vector of group labels (1-based).
#' @param K Number of groups.
#' @return A \code{K}-by-\code{ncol(x)} matrix of group means.
#'
#' @examples
#' x <- matrix(1:12, 6, 2)
#' labels <- c(1L, 2L, 1L, 3L, 2L, 1L)
#' # One row of column means per cluster
#' am_scatter_mean(x, labels, K = 3L)
#'
#' @export
am_scatter_mean <- function(x, labels, K) {
  labels <- as.integer(labels)
  K      <- as.integer(K)
  if (length(labels) != nrow(x)) {
    stop("length(labels) must equal nrow(x)", call. = FALSE)
  }
  counts <- tabulate(labels, nbins = K)   # O(n), always on CPU

  if (!inherits(x, "adgeMatrix")) {
    return(.am_scatter_mean_cpu(as.matrix(x), labels, K))
  }

  choice <- .amatrix_backend_for(x, "scatter_mean")
  sums   <- .amatrix_try_resident_scatter_mean(x, labels, K, choice$name)

  if (is.null(sums)) {
    return(.am_scatter_mean_cpu(as.matrix(amatrix_materialize_host(x)), labels, K))
  }

  # sums is K×p; divide each row by cluster count
  means <- sums
  nonzero <- counts > 0L
  if (any(nonzero))
    means[nonzero, ] <- sums[nonzero, , drop = FALSE] / counts[nonzero]
  means[!nonzero, ] <- NA_real_
  means
}

.amatrix_try_resident_broadcast_ewise <- function(x, v, margin, op, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_supports_resident_op(backend, "broadcast_ewise", x = x)) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(x, backend_name)
  if (is.null(lhs)) return(NULL)
  defer <- .amatrix_defer_host_active()
  out_key <- .amatrix_next_resident_key(backend_name)
  value <- backend$broadcast_ewise_resident(lhs$key, v, margin, op, out_key,
                                             defer = defer)
  .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  list(value = value, backend = backend_name, resident_key = out_key)
}

#' Backend-dispatched sweep
#'
#' Apply a function to each row or column of a matrix, dispatching to the
#' preferred GPU backend when available.
#'
#' @param x A numeric matrix or \code{adgeMatrix}.
#' @param MARGIN 1 for rows, 2 for columns.
#' @param STATS Numeric vector of statistics to apply.
#' @param FUN Operation: \code{"+"}, \code{"-"} (the default, as for
#'   \code{\link[base]{sweep}}), \code{"*"}, or \code{"/"}.
#' @return A matrix of the same dimensions as \code{x}.
#'
#' @examples
#' x <- adgeMatrix(matrix(1:6, 2, 3))
#' # Subtract the column means
#' am_sweep(x, 2, colMeans(as.matrix(x)), FUN = "-")
#'
#' @seealso \code{\link{am_sweep_inplace}}
#' @export
am_sweep <- function(x, MARGIN, STATS, FUN = "-") {
  if (inherits(x, c("adgCMatrix", "sparseMatrix"))) {
    scaled <- .amatrix_sparse_scale(x, MARGIN, STATS, FUN)
    if (!is.null(scaled)) return(scaled)
  }
  if (!inherits(x, "adgeMatrix")) {
    return(base::sweep(as.matrix(x), MARGIN = MARGIN, STATS = STATS, FUN = FUN))
  }
  op <- if (is.character(FUN) && length(FUN) == 1L) FUN else NULL
  if (is.null(op) || !op %in% c("+", "-", "*", "/")) {
    return(base::sweep(as.matrix(amatrix_materialize_host(x)), MARGIN = MARGIN,
                       STATS = STATS, FUN = FUN))
  }
  if (!is.numeric(STATS) || is.matrix(STATS)) {
    return(base::sweep(as.matrix(amatrix_materialize_host(x)), MARGIN = MARGIN,
                       STATS = STATS, FUN = FUN))
  }
  v <- as.double(STATS)
  choice <- .amatrix_backend_for(x, "broadcast_ewise")
  resident <- .amatrix_try_resident_broadcast_ewise(x, v, MARGIN, op, choice$name)
  if (!is.null(resident)) {
    return(.amatrix_resident_wrap(x, resident, out_dim = dim(x)))
  }
  result <- amatrix_dispatch_op(
    x = x,
    op = "broadcast_ewise",
    method = "broadcast_ewise",
    args = list(lhs = as.matrix(amatrix_materialize_host(x)), v = v,
                margin = MARGIN, op = op),
    fallback = function() base::sweep(as.matrix(amatrix_materialize_host(x)),
                                      MARGIN = MARGIN, STATS = STATS, FUN = FUN)
  )
  .amatrix_rewrap_value(x, result)
}

.amatrix_try_resident_argreduce <- function(x, kind, backend_name) {
  backend  <- .amatrix_get_backend(backend_name)
  fn_name  <- paste0(kind, "_resident")
  if (!is.function(backend[[fn_name]])) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(x, backend_name)
  if (is.null(lhs)) return(NULL)
  result <- backend[[fn_name]](lhs$key)
  .amatrix_cleanup_temp_resident(list(lhs), backend_name)
  result
}

.am_argreduce_cpu <- function(x, kind) {
  mat <- as.matrix(x)
  switch(kind,
    rowargmax = max.col(mat,  ties.method = "first"),
    rowargmin = max.col(-mat, ties.method = "first"),
    colargmax = max.col(t(mat),  ties.method = "first"),
    colargmin = max.col(-t(mat), ties.method = "first")
  )
}

.am_argreduce <- function(x, kind) {
  if (!inherits(x, "adgeMatrix")) return(.am_argreduce_cpu(x, kind))
  choice <- .amatrix_backend_for(x, "argmax")
  result <- .amatrix_try_resident_argreduce(x, kind, choice$name)
  if (!is.null(result)) return(result)
  .am_argreduce_cpu(x, kind)
}

#' Row and column argmax/argmin
#'
#' Return the index of the maximum or minimum value in each row or column.
#'
#' @param x A numeric matrix or \code{adgeMatrix}.
#' @return An integer vector of indices.
#'
#' @examples
#' x <- adgeMatrix(matrix(c(3, 1, 2, 9, 7, 8), 2, 3))
#' am_rowargmax(x)
#' am_rowargmin(x)
#' am_colargmax(x)
#' am_colargmin(x)
#'
#' @name am_argreduce
#' @rdname am_argreduce
#' @export
am_rowargmax <- function(x) .am_argreduce(x, "rowargmax")

#' @rdname am_argreduce
#' @export
am_rowargmin <- function(x) .am_argreduce(x, "rowargmin")

#' @rdname am_argreduce
#' @export
am_colargmax <- function(x) .am_argreduce(x, "colargmax")

#' @rdname am_argreduce
#' @export
am_colargmin <- function(x) .am_argreduce(x, "colargmin")

#' Element-wise operations
#'
#' Apply an element-wise arithmetic operation to one or two matrices,
#' dispatching to the preferred GPU backend when available.
#'
#' @param op Character string: \code{"+"}, \code{"-"}, \code{"*"}, or \code{"/"}.
#' @param e1 A numeric matrix or \code{adgeMatrix}.
#' @param e2 A numeric matrix, \code{adgeMatrix}, or \code{NULL} for unary ops.
#' @return A matrix of the same dimensions as \code{e1}.
#'
#' @examples
#' x <- adgeMatrix(matrix(1:6, 2, 3))
#' am_ewise("+", x, 1)
#' am_ewise("*", x, x)
#' am_ewise("exp", x)
#'
#' @export
am_ewise <- function(op, e1, e2 = NULL) {
  if (!inherits(e1, "aMatrix") && !inherits(e2, "aMatrix")) {
    e1 <- .amatrix_promote_arg(e1)
  }
  template <- .amatrix_template(e1, e2)

  # Try GPU resident path FIRST — no host materialization needed
  if (!is.null(template)) {
    choice <- .amatrix_backend_for(template, "ewise", y = e2)
    resident <- .amatrix_try_resident_ewise(op, e1, e2, choice$name)
    if (!is.null(resident)) {
      return(.amatrix_resident_wrap(template, resident, out_dim = dim(template)))
    }
  }

  # Only materialize on the cold/fallback path
  host_e1 <- .amatrix_host_arg(e1)
  host_e2 <- .amatrix_host_arg(e2)

  value <- amatrix_dispatch_op(
    x = template,
    op = "ewise",
    method = "ewise",
    y = e2,
    args = list(lhs = host_e1, rhs = host_e2, op = op),
    fallback = function() {
      if (is.null(e2)) {
        return(do.call(op, list(host_e1)))
      }
      do.call(op, list(host_e1, host_e2))
    }
  )

  if (is.null(template)) {
    return(value)
  }

  .amatrix_rewrap_value(template, value)
}

am_set_dimnames <- function(x, value) {
  host_x <- amatrix_materialize_host(x)
  dimnames(host_x) <- value
  .amatrix_rewrap_like(x, host_x)
}

# ── Distance / Kernel helpers ──────────────────────────────────────────────

.am_as_double_matrix <- function(x) {
  if (!(inherits(x, c("aMatrix", "Matrix")) || is.matrix(x))) {
    stop("x must be a matrix-like object", call. = FALSE)
  }
  if (inherits(x, c("adgeMatrix", "adgCMatrix")))
    x <- amatrix_materialize_host(x)
  if (!is.matrix(x)) x <- as.matrix(x)
  if (length(dim(x)) != 2L) {
    stop("x must be a matrix-like object", call. = FALSE)
  }
  if (is.complex(x)) {
    stop("complex matrices are not supported", call. = FALSE)
  }
  if (!is.double(x)) storage.mode(x) <- "double"
  x
}

.am_metric_checked_matrix <- function(x, name) {
  mat <- .am_as_double_matrix(x)
  if (anyNA(mat) || any(!is.finite(mat))) {
    stop(sprintf("%s must contain only finite non-missing values", name), call. = FALSE)
  }
  mat
}

.am_kernel_finalize <- function(out, kernel, y_host, zero_diag) {
  if (is.null(y_host) && identical(kernel, "rbf")) {
    diag(out) <- 1
  }
  if (isTRUE(zero_diag) && is.null(y_host)) {
    diag(out) <- 0
  }
  out
}

# GPU dispatch helpers for distance/kernel computation.
# .dist_matrix_sq_gpu: returns squared Euclidean distance matrix [m×n].
# .am_kernel_gpu:  returns kernel matrix [m×n].
# Both use dedicated column-major AF bridges (all dims) or MLX.

.am_af_ok <- function() {
  tryCatch(
    requireNamespace("amatrix.arrayfire", quietly = TRUE) &&
      amatrix.arrayfire::amatrix_arrayfire_is_available(),
    error = function(e) FALSE
  )
}
.am_mlx_ok <- function() {
  tryCatch(
    requireNamespace("amatrix.mlx", quietly = TRUE) &&
      amatrix.mlx::amatrix_mlx_is_available(),
    error = function(e) FALSE
  )
}

.am_metric_backend_spec <- function(X, Y = NULL, preferred_backend = NULL) {
  if (!is.null(preferred_backend)) {
    return(list(
      name = preferred_backend,
      precision = amatrix_default_precision(),
      policy = amatrix_default_policy()
    ))
  }

  template <- NULL
  if (inherits(X, "aMatrix")) {
    template <- X
  } else if (!is.null(Y) && inherits(Y, "aMatrix")) {
    template <- Y
  }

  if (is.null(template)) {
    return(list(
      name = "auto",
      precision = amatrix_default_precision(),
      policy = amatrix_default_policy()
    ))
  }

  list(
    name = template@preferred_backend,
    precision = template@precision,
    policy = template@policy
  )
}

.am_metric_wrap_dense_operand <- function(x, spec) {
  if (inherits(x, "adgeMatrix")) {
    return(x)
  }

  adgeMatrix(
    .am_as_double_matrix(x),
    preferred_backend = spec$name,
    precision = spec$precision,
    policy = spec$policy
  )
}

.am_metric_tcrossprod_backend <- function(X, Y = NULL, spec = .am_metric_backend_spec(X, Y)) {
  x_host <- .am_as_double_matrix(X)
  y_host <- if (is.null(Y)) NULL else .am_as_double_matrix(Y)

  if (identical(spec$name, "cpu")) {
    return(if (is.null(y_host)) tcrossprod(x_host) else tcrossprod(x_host, y_host))
  }

  if (identical(spec$name, "arrayfire") && isTRUE(.am_af_ok())) {
    return(getExportedValue("amatrix.arrayfire", "amatrix_arrayfire_tcrossprod_host_bridge")(x_host, y_host))
  }

  if (identical(spec$name, "mlx") && isTRUE(.am_mlx_ok())) {
    return(getExportedValue("amatrix.mlx", "amatrix_mlx_tcrossprod_host_bridge")(x_host, y_host))
  }

  if (identical(spec$name, "opencl")) {
    x_arg <- .am_metric_wrap_dense_operand(X, spec)
    y_arg <- if (is.null(Y)) NULL else .am_metric_wrap_dense_operand(Y, spec)
    return(as.matrix(amatrix_materialize_host(am_tcrossprod(x_arg, y_arg))))
  }

  if (identical(spec$name, "auto")) {
    if (isTRUE(.am_af_ok())) {
      return(getExportedValue("amatrix.arrayfire", "amatrix_arrayfire_tcrossprod_host_bridge")(x_host, y_host))
    }
    if (isTRUE(.am_mlx_ok())) {
      return(getExportedValue("amatrix.mlx", "amatrix_mlx_tcrossprod_host_bridge")(x_host, y_host))
    }
  }

  if (is.null(y_host)) tcrossprod(x_host) else tcrossprod(x_host, y_host)
}

# MLX is preferred on Apple Silicon (39x vs 9x speedup in benchmarks).
# AF is the fallback for CUDA/other platforms where MLX is unavailable.

# AF bridge does D² entirely in C (GEMM + rowSums + broadcast + clamp) with no
# R-level allocation overhead, so it's fast regardless of GPU vs CPU backend.
# MLX path only does the GEMM on GPU; subsequent R ops on 3M-element matrices
# add ~150ms overhead for large matrices — only worthwhile when AF unavailable.
.dist_matrix_sq_gpu <- function(X, Y = NULL, spec = .am_metric_backend_spec(X, Y)) {
  x_host <- .am_as_double_matrix(X)
  y_host <- if (is.null(Y)) NULL else .am_as_double_matrix(Y)
  y_eff <- if (is.null(y_host)) x_host else y_host

  if (identical(spec$name, "arrayfire") && isTRUE(.am_af_ok())) {
    return(getExportedValue("amatrix.arrayfire", "amatrix_arrayfire_dist_sq_host_bridge")(x_host, y_host))
  }

  if (identical(spec$name, "auto") && isTRUE(.am_af_ok())) {
    return(getExportedValue("amatrix.arrayfire", "amatrix_arrayfire_dist_sq_host_bridge")(x_host, y_host))
  }

  G <- .am_metric_tcrossprod_backend(X, Y, spec = spec)
  nx <- rowSums(x_host^2)
  ny <- if (is.null(y_host)) nx else rowSums(y_eff^2)
  pmax(outer(nx, ny, "+") - 2 * G, 0)
}

.am_kernel_gpu <- function(X, Y = NULL, kernel, sigma, degree, coef, zero_diag = FALSE, spec = .am_metric_backend_spec(X, Y)) {
  x_host <- .am_as_double_matrix(X)
  y_host <- if (is.null(Y)) NULL else .am_as_double_matrix(Y)
  y_eff <- if (is.null(y_host)) x_host else y_host

  # AF bridge computes entire kernel in C — no R allocation overhead.
  if ((identical(spec$name, "arrayfire") || identical(spec$name, "auto")) && isTRUE(.am_af_ok())) {
    out <- getExportedValue("amatrix.arrayfire", "amatrix_arrayfire_kernel_host_bridge")(
      x_host, y_host, kernel,
      as.double(sigma), as.integer(degree), as.double(coef)
    )
    return(.am_kernel_finalize(out, kernel, y_host, zero_diag))
  }
  # MLX: GPU GEMM + R-level transforms (cheap for linear/poly/cosine; heavier for rbf/lap)
  if (identical(spec$name, "mlx") && isTRUE(.am_mlx_ok())) {
    G <- getExportedValue("amatrix.mlx", "amatrix_mlx_tcrossprod_host_bridge")(x_host, y_host)
    return(.am_kernel_finalize(switch(kernel,
      linear     = G,
      polynomial = (coef + G)^degree,
      cosine     = {
        nx <- sqrt(rowSums(x_host^2)); ny <- if (is.null(y_host)) nx else sqrt(rowSums(y_eff^2))
        G / pmax(outer(nx, ny), .Machine$double.eps)
      },
      rbf        = {
        nx <- rowSums(x_host^2); ny <- if (is.null(y_host)) nx else rowSums(y_eff^2)
        D_sq <- pmax(outer(nx, ny, "+") - 2 * G, 0)
        if (is.null(y_host)) diag(D_sq) <- 0
        out <- exp(-D_sq / (2 * sigma^2))
        if (isTRUE(zero_diag) && is.null(y_host)) diag(out) <- 0
        out
      },
      laplacian  = {
        nx <- rowSums(x_host^2); ny <- if (is.null(y_host)) nx else rowSums(y_eff^2)
        D_sq <- pmax(outer(nx, ny, "+") - 2 * G, 0)
        if (is.null(y_host)) diag(D_sq) <- 0
        out <- exp(-sqrt(D_sq) / sigma)
        if (isTRUE(zero_diag) && is.null(y_host)) diag(out) <- 0
        out
      }
    ), kernel, y_host, zero_diag))
  }
  # CPU fallback
  G <- .am_metric_tcrossprod_backend(X, Y, spec = spec)

  .am_kernel_finalize(switch(kernel,
    linear     = G,
    polynomial = (coef + G)^degree,
    cosine     = {
      nx <- sqrt(rowSums(x_host^2))
      ny <- if (is.null(y_host)) nx else sqrt(rowSums(y_eff^2))
      G / pmax(outer(nx, ny), .Machine$double.eps)
    },
    rbf        = {
      nx   <- rowSums(x_host^2); ny <- if (is.null(y_host)) nx else rowSums(y_eff^2)
      D_sq <- pmax(outer(nx, ny, "+") - 2 * G, 0)
      if (is.null(y_host)) diag(D_sq) <- 0
      out <- exp(-D_sq / (2 * sigma^2))
      if (isTRUE(zero_diag) && is.null(y_host)) diag(out) <- 0
      out
    },
    laplacian  = {
      nx   <- rowSums(x_host^2); ny <- if (is.null(y_host)) nx else rowSums(y_eff^2)
      D_sq <- pmax(outer(nx, ny, "+") - 2 * G, 0)
      if (is.null(y_host)) diag(D_sq) <- 0
      out <- exp(-sqrt(D_sq) / sigma)
      if (isTRUE(zero_diag) && is.null(y_host)) diag(out) <- 0
      out
    }
  ), kernel, y_host, zero_diag)
}

# Tiled pairwise distance for large n (avoids GPU OOM on n > 50k).
# Processes row-blocks of X (and Y) independently, assembling the host result
# block by block.  Exploits symmetry when Y = NULL to halve the GEMM count.
.dist_matrix_tiled <- function(X, Y, method, tile_size, spec = .am_metric_backend_spec(X, Y)) {
  m <- nrow(X)
  symmetric <- is.null(Y)
  Y_eff <- if (symmetric) X else Y
  n <- nrow(Y_eff)

  result <- matrix(0.0, nrow = m, ncol = n)

  i_breaks <- c(seq(1L, m, by = tile_size), m + 1L)
  j_breaks <- c(seq(1L, n, by = tile_size), n + 1L)
  ni <- length(i_breaks) - 1L
  nj <- length(j_breaks) - 1L

  for (ii in seq_len(ni)) {
    i0 <- i_breaks[ii]; i1 <- i_breaks[ii + 1L] - 1L
    Xi <- X[i0:i1, , drop = FALSE]

    # Symmetric: compute only lower-triangular blocks, mirror to upper.
    j_end <- if (symmetric) ii else nj
    for (jj in seq_len(j_end)) {
      j0 <- j_breaks[jj]; j1 <- j_breaks[jj + 1L] - 1L
      Xj <- Y_eff[j0:j1, , drop = FALSE]

      D_sq <- .dist_matrix_sq_gpu(Xi, Xj, spec = spec)
      D_block <- if (method == "euclidean") sqrt(D_sq) else D_sq

      result[i0:i1, j0:j1] <- D_block
      if (symmetric && ii != jj) result[j0:j1, i0:i1] <- t(D_block)
    }
  }

  if (symmetric) diag(result) <- 0
  result
}

#' GPU-accelerated pairwise distance matrix
#'
#' Computes the pairwise distance matrix between rows of \code{X} and \code{Y}.
#' The dominant cost (row inner-products via am_tcrossprod) is dispatched to the
#' active GPU backend (ArrayFire or MLX); norm computation and final transforms
#' run on CPU where they are O(mp + np) — negligible versus the O(mnp) GEMM.
#'
#' @param X Numeric matrix or \code{adgeMatrix}, shape [m, p].
#' @param Y Numeric matrix or \code{adgeMatrix}, shape [n, p], or \code{NULL}
#'   to compute pairwise distances within \code{X} (returns [m, m] matrix).
#' @param method One of \code{"euclidean"} (default), \code{"sqeuclidean"},
#'   or \code{"cosine"} (cosine distance,
#'   \eqn{1 - x^\top y / (\|x\| \|y\|)}; for the similarity itself use
#'   \code{\link{am_kernel_matrix}(kernel = "cosine")}).
#' @param tile_size Integer row-block size for tiled computation, or \code{NULL}
#'   (default) to auto-tile when \code{nrow(X) > 50000} (self-distance only).
#'   Set explicitly to process any size in row-blocks; useful when GPU memory
#'   is limited.  Not supported for \code{method = "cosine"}.
#' @return Numeric matrix [m, n] of pairwise distances.
#'
#' @examples
#' X <- matrix(rnorm(30), nrow = 6)
#' D <- am_dist_matrix(X)
#' dim(D)
#'
#' @seealso \code{\link{am_kernel_matrix}}
#' @export
am_dist_matrix <- function(X, Y = NULL,
                    method = c("euclidean", "sqeuclidean", "cosine"),
                    tile_size = NULL) {
  method <- match.arg(method)
  spec <- .am_metric_backend_spec(X, Y)
  X_mat <- .am_metric_checked_matrix(X, "X")
  Y_mat <- if (!is.null(Y)) .am_metric_checked_matrix(Y, "Y") else NULL
  if (any(dim(X_mat) == 0L) || (!is.null(Y_mat) && any(dim(Y_mat) == 0L))) {
    spec$name <- "cpu"
  }

  # Auto-tile for large self-distance to prevent GPU OOM
  if (is.null(tile_size) && is.null(Y_mat) && nrow(X_mat) > 50000L)
    tile_size <- 10000L

  if (!is.null(tile_size) && method != "cosine") {
    return(.dist_matrix_tiled(X_mat, Y_mat, method, as.integer(tile_size), spec = spec))
  }

  if (method == "cosine") {
    # Cosine distance = 1 - cosine similarity (the similarity itself is
    # am_kernel_matrix(kernel = "cosine")). Clamp rounding noise into [0, 2] and
    # make self-distances exactly zero.
    sim <- .am_kernel_gpu(X, Y, "cosine", 1.0, 2L, 0.0, spec = spec)
    D <- 1 - as.matrix(sim)
    D[] <- pmin(pmax(D, 0), 2)
    if (is.null(Y_mat) && nrow(D) == ncol(D)) diag(D) <- 0
    dimnames(D) <- dimnames(sim)
    return(if (inherits(sim, "aMatrix")) .amatrix_rewrap_like(sim, D) else D)
  }

  D_sq <- .dist_matrix_sq_gpu(X, Y, spec = spec)
  if (!is.null(Y_mat) && identical(X_mat, Y_mat)) diag(D_sq) <- 0
  if (is.null(Y_mat)) diag(D_sq) <- 0   # fix float32/float64 diagonal mismatch
  if (method == "sqeuclidean") return(D_sq)
  sqrt(D_sq)
}

#' GPU-accelerated pairwise kernel matrix
#'
#' Computes the pairwise kernel matrix between rows of \code{X} and \code{Y}.
#' The expensive am_tcrossprod is GPU-dispatched; element-wise transforms (exp,
#' sqrt, pow) run on CPU.
#'
#' Kernels:
#' \describe{
#'   \item{linear}{\eqn{k(x,y) = x^\top y}}
#'   \item{rbf}{\eqn{k(x,y) = \exp(-\|x - y\|^2 / (2 \sigma^2))}}
#'   \item{polynomial}{\eqn{k(x,y) = (\mathrm{coef} + x^\top y)^{\mathrm{degree}}}}
#'   \item{cosine}{\eqn{k(x,y) = (x^\top y) / (\|x\| \|y\|)}}
#'   \item{laplacian}{\eqn{k(x,y) = \exp(-\|x - y\| / \sigma)}}
#' }
#'
#' @param X Numeric matrix or \code{adgeMatrix}, shape [m, p].
#' @param Y Numeric matrix or \code{adgeMatrix}, shape [n, p], or \code{NULL}.
#' @param kernel Kernel type string (see Details).
#' @param sigma Bandwidth for \code{"rbf"} and \code{"laplacian"}.
#' @param degree Polynomial degree for \code{"polynomial"}.
#' @param coef Constant term for \code{"polynomial"}: \eqn{(\mathrm{coef} + x^\top y)^{\mathrm{degree}}}.
#' @param preferred_backend Optional backend name to override the default
#'   dispatch (e.g., \code{"mlx"}, \code{"opencl"}).
#' @param zero_diag When \code{TRUE} and \code{Y} is \code{NULL}, set the
#'   diagonal of the kernel matrix to zero.
#' @return Numeric matrix [m, n] of kernel values.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(20), 5, 4)
#' K <- am_kernel_matrix(X, kernel = "rbf", sigma = 1.5)
#' dim(K)
#' # Cross-kernel between two sets of rows
#' am_kernel_matrix(X, X[1:2, ], kernel = "linear")
#'
#' @seealso \code{\link{am_dist_matrix}}
#' @export
am_kernel_matrix <- function(X, Y = NULL,
                      kernel = c("linear", "rbf", "polynomial",
                                 "cosine", "laplacian"),
                      sigma = 1.0, degree = 2L, coef = 0.0,
                      preferred_backend = NULL, zero_diag = FALSE) {
  kernel <- match.arg(kernel)
  if (kernel %in% c("rbf", "laplacian")) {
    if (!is.numeric(sigma) || length(sigma) != 1L || is.na(sigma) || !is.finite(sigma) || sigma <= 0) {
      stop("sigma must be a single positive finite number", call. = FALSE)
    }
  }
  spec <- .am_metric_backend_spec(X, Y, preferred_backend = preferred_backend)
  X_mat  <- .am_metric_checked_matrix(X, "X")
  Y_mat  <- if (!is.null(Y)) .am_metric_checked_matrix(Y, "Y") else NULL
  if (any(dim(X_mat) == 0L) || (!is.null(Y_mat) && any(dim(Y_mat) == 0L))) {
    spec$name <- "cpu"
  }

  # Resident path: compute on GPU and store directly, skipping CPU round-trip.
  # Returns an adgeMatrix with a live resident key bound, so the next GPU op
  # uses the in-device array without re-uploading.
  if (!is.null(preferred_backend)) {
    backend <- .amatrix_get_backend(preferred_backend)
    if (!is.null(backend) && is.function(backend[["kernel_resident"]])) {
      out_key    <- .amatrix_next_resident_key(preferred_backend)
      result_mat <- backend$kernel_resident(
        out_key, X_mat, Y_mat, kernel,
        as.double(sigma), as.integer(degree), as.double(coef),
        isTRUE(zero_diag) && is.null(Y_mat)
      )
      obj <- new_adgeMatrix(result_mat, preferred_backend = preferred_backend)
      return(.amatrix_bind_resident(obj, preferred_backend, out_key))
    }
  }

  .am_kernel_gpu(
    X,
    Y,
    kernel,
    sigma = sigma,
    degree = degree,
    coef = coef,
    zero_diag = zero_diag,
    spec = spec
  )
}

# ---------------------------------------------------------------------------
# Kronecker product
# ---------------------------------------------------------------------------


# norm: matrix and vector norms.
#
# Matrix types follow base::norm exactly, including its default ("O"):
#   "O"/"1" - max absolute column sum (one norm)
#   "I"     - max absolute row sum (infinity norm)
#   "F"/"E" - Frobenius norm
#   "M"     - max absolute entry
#   "2"     - spectral norm (largest singular value), computed exactly
#
# For numeric vector inputs (which base::norm rejects), "2"/"F"/"E" return the
# Euclidean norm (the default), "O"/"1" return sum(|x|), "I"/"M" max(|x|).
#
# norm() is an S4 generic that dispatches on (x, type), and S4 ignores default
# values in method definitions, so each class needs an explicit
# type = "missing" method. No method is registered for plain base matrices:
# doing so would change norm() for every matrix in the session.

.norm_type <- function(type) {
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    stop("'type' must be a single character string", call. = FALSE)
  }
  key <- toupper(substr(type, 1L, 1L))
  if (!nzchar(key) || !key %in% c("O", "1", "I", "F", "E", "M", "2")) {
    stop(sprintf("invalid 'type' for norm(): \"%s\"", type), call. = FALSE)
  }
  switch(key, "1" = "O", "E" = "F", key)
}

.amatrix_dense_norm <- function(x, type) {
  base::norm(.am_as_double_matrix(x), type = .norm_type(type))
}

.amatrix_vector_norm <- function(x, type) {
  switch(.norm_type(type),
    "2" = ,
    "F" = sqrt(sum(x^2)),
    "O" = sum(abs(x)),
    "I" = ,
    "M" = max(abs(x))
  )
}

#' @noRd
setMethod("norm", signature(x = "adgeMatrix", type = "missing"), function(x, type, ...) {
  .amatrix_dense_norm(x, "O")
})

#' @noRd
setMethod("norm", signature(x = "adgeMatrix", type = "character"), function(x, type, ...) {
  .amatrix_dense_norm(x, type)
})

# Vector norms (base::norm errors on a bare numeric vector).
#' @noRd
setMethod("norm", signature(x = "numeric", type = "missing"), function(x, type, ...) {
  .amatrix_vector_norm(x, "2")
})

#' @noRd
setMethod("norm", signature(x = "numeric", type = "character"), function(x, type, ...) {
  .amatrix_vector_norm(x, type)
})

# ---------------------------------------------------------------------------
# Batch factorization helpers
# ---------------------------------------------------------------------------

# Normalise a batch argument to a list of plain R matrices.
# Accepts: list of matrices, or a 3-D array dim c(n, n, B).
.am_batch_to_list <- function(A, arg_name = "A") {
  if (is.list(A)) {
    lapply(A, function(a) {
      m <- as.matrix(a)
      if (!is.double(m)) storage.mode(m) <- "double"
      m
    })
  } else if (is.array(A) && length(dim(A)) == 3L) {
    B <- dim(A)[[3L]]
    lapply(seq_len(B), function(b) {
      m <- A[, , b, drop = FALSE]
      dim(m) <- dim(m)[1:2]
      if (!is.double(m)) storage.mode(m) <- "double"
      m
    })
  } else {
    stop(arg_name, " must be a list of matrices or a 3-D array [n, n, B]",
         call. = FALSE)
  }
}

#' Batch Cholesky factorization
#'
#' Factorize B symmetric positive-definite matrices in parallel.  Each matrix
#' is dispatched through the same backend as \code{\link{am_chol_factor}}, so MLX
#' GPU acceleration applies to every element when available.
#'
#' @param A A list of square numeric matrices, or a 3-D array \code{[n, n, B]}.
#' @return A list of \code{amChol} objects, one per input matrix.
#'
#' @examples
#' mats <- lapply(1:3, function(i) crossprod(matrix(rnorm(12), 4, 3)) + diag(3))
#' facs <- am_batch_chol(mats)
#' length(facs)
#'
#' @seealso \code{\link{am_chol_factor}}, \code{\link{am_batch_solve}}
#' @export
am_batch_chol <- function(A) {
  mats <- .am_batch_to_list(A, "A")
  lapply(mats, function(m) {
    X <- adgeMatrix(m)
    am_chol_factor(X)
  })
}

#' Batch triangular solve
#'
#' Solve B linear systems \code{A_b x_b = B_b} where each \code{A_b} is
#' represented by its Cholesky factor from \code{\link{am_batch_chol}}.
#'
#' @param Ls A list of \code{amChol} objects (output of \code{am_batch_chol}).
#' @param B A list of right-hand-side matrices/vectors, or a 3-D array
#'   \code{[n, k, B]}.  Length / third dimension must match \code{Ls}.
#' @return A list of solution matrices (or vectors when each rhs is a vector).
#'
#' @examples
#' mats <- lapply(1:2, function(i) crossprod(matrix(rnorm(12), 4, 3)) + diag(3))
#' facs <- am_batch_chol(mats)
#' rhs <- lapply(1:2, function(i) rnorm(3))
#' sols <- am_batch_solve(facs, rhs)
#' all.equal(sols[[1]], solve(mats[[1]], rhs[[1]]))
#'
#' @seealso \code{\link{am_batch_chol}}, \code{\link{am_chol_solve}}
#' @export
am_batch_solve <- function(Ls, B) {
  if (!is.list(Ls) || !all(vapply(Ls, inherits, logical(1L), "amChol"))) {
    stop("Ls must be a list of amChol objects from am_batch_chol()", call. = FALSE)
  }

  rhs_list <- if (is.list(B)) {
    B
  } else if (is.array(B) && length(dim(B)) == 3L) {
    nb <- dim(B)[[3L]]
    lapply(seq_len(nb), function(b) B[, , b, drop = TRUE])
  } else {
    stop("B must be a list or a 3-D array [n, k, B]", call. = FALSE)
  }

  if (length(Ls) != length(rhs_list)) {
    stop("Ls and B must have the same batch size", call. = FALSE)
  }

  Map(am_chol_solve, Ls, rhs_list)
}

#' Eager Kronecker product
#'
#' Computes the Kronecker product \eqn{A \otimes B} and returns the result as an \code{adgeMatrix}.
#' Accepts plain matrices or any \code{aMatrix} subclass.
#' For a lazy variant that avoids forming the full product see
#' \code{\link{kron_matrix}}.
#'
#' @param A,B Matrices or \code{aMatrix} objects.
#' @return An \code{adgeMatrix} of dimension \code{(nrow(A)*nrow(B)) x (ncol(A)*ncol(B))}.
#'
#' @examples
#' A <- matrix(1:4, 2, 2)
#' B <- diag(2)
#' am_kron(A, B)
#'
#' @seealso \code{\link{kron_matrix}}
#' @export
am_kron <- function(A, B) {
  A_mat <- .am_as_double_matrix(A)
  B_mat <- .am_as_double_matrix(B)
  adgeMatrix(base::kronecker(A_mat, B_mat))
}

#' Batch crossproduct
#'
#' Compute \code{t(A_b) \%*\% A_b} for each matrix in a batch.
#'
#' @param A A list of numeric matrices, or a 3-D array \code{[n, p, B]}.
#' @return A list of \code{p x p} crossproduct matrices.
#'
#' @examples
#' mats <- lapply(1:3, function(i) matrix(rnorm(12), 4, 3))
#' xtx <- am_batch_crossprod(mats)
#' dim(xtx[[1]])
#'
#' @seealso \code{\link{am_batch_chol}}
#' @export
am_batch_crossprod <- function(A) {
  mats <- .am_batch_to_list(A, "A")
  lapply(mats, function(m) crossprod(m))
}
