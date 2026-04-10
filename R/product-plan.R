.amatrix_product_plan_meta <- function(plan) {
  attr(plan, "amatrix_plan_meta", exact = TRUE)
}

.amatrix_product_plan_lhs <- function(plan) {
  if (!inherits(plan, "am_product_plan")) {
    return(NULL)
  }

  env <- environment(plan)
  if (!is.environment(env) || !exists("lhs_bound", envir = env, inherits = FALSE)) {
    return(NULL)
  }

  get("lhs_bound", envir = env, inherits = FALSE)
}

.amatrix_release_product_plan <- function(plan) {
  meta <- .amatrix_product_plan_meta(plan)
  if (!inherits(plan, "am_product_plan") || !isTRUE(meta$owned_resident)) {
    return(invisible(FALSE))
  }

  lhs <- .amatrix_product_plan_lhs(plan)
  if (!inherits(lhs, "aMatrix")) {
    return(invisible(FALSE))
  }

  .amatrix_release_resident(lhs)
}

.amatrix_product_plan_drop_output <- function(resident) {
  if (is.null(resident) || is.null(resident$backend) || is.null(resident$resident_key)) {
    return(invisible(FALSE))
  }

  backend <- tryCatch(.amatrix_get_backend(resident$backend), error = function(e) NULL)
  if (is.null(backend) || !is.function(backend$resident_has) || !is.function(backend$resident_drop)) {
    return(invisible(FALSE))
  }

  if (isTRUE(backend$resident_has(resident$resident_key))) {
    backend$resident_drop(resident$resident_key)
    return(invisible(TRUE))
  }

  invisible(FALSE)
}

.amatrix_product_plan_sparse_host_matrix <- function(lhs_bound, y, op, lhs_host = NULL) {
  if (!inherits(lhs_bound, "adgCMatrix")) {
    return(NULL)
  }

  rhs <- if (is.numeric(y) && is.null(dim(y))) matrix(y, ncol = 1L) else .amatrix_host_arg(y)
  rhs_mat <- if (is.matrix(rhs)) rhs else as.matrix(rhs)
  if (!is.double(rhs_mat)) {
    storage.mode(rhs_mat) <- "double"
  }

  if (is.null(lhs_host)) {
    lhs_host <- amatrix_materialize_host(lhs_bound)
  }
  value <- switch(
    op,
    matmul = lhs_host %*% rhs_mat,
    crossprod = Matrix::crossprod(lhs_host, rhs_mat),
    tcrossprod = Matrix::tcrossprod(lhs_host, rhs_mat)
  )
  as.matrix(value)
}

.amatrix_product_plan_matrix_result <- function(lhs_bound, y, op, backend_name, lhs_host = NULL) {
  if (is.null(backend_name) || identical(backend_name, "cpu") || !inherits(lhs_bound, "aMatrix")) {
    return(NULL)
  }

  backend <- tryCatch(.amatrix_get_backend(backend_name), error = function(e) NULL)
  if (is.null(backend)) {
    return(NULL)
  }
  if (!.amatrix_backend_supports_resident_op(backend, op, x = lhs_bound, y = y)) {
    return(.amatrix_product_plan_sparse_host_matrix(lhs_bound, y, op, lhs_host = lhs_host))
  }

  if (inherits(lhs_bound, "adgCMatrix") && is.function(backend$spmm_resident)) {
    lhs <- .amatrix_prepare_resident_arg(lhs_bound, backend_name, promote_amatrix = TRUE)
    if (!is.null(lhs) && isTRUE(lhs$sparse)) {
      rhs <- if (is.numeric(y) && is.null(dim(y))) matrix(y, ncol = 1L) else .amatrix_host_arg(y)
      rhs_mat <- if (is.matrix(rhs)) rhs else as.matrix(rhs)
      if (!is.double(rhs_mat)) {
        storage.mode(rhs_mat) <- "double"
      }
      trans_lhs <- identical(op, "crossprod")
      if (identical(op, "tcrossprod")) {
        rhs_mat <- t(rhs_mat)
      }

      value <- tryCatch(
        backend$spmm_resident(lhs$key, rhs_mat, trans_lhs = trans_lhs),
        error = function(e) NULL
      )
      .amatrix_cleanup_temp_resident(list(lhs), backend_name)
      if (!is.null(value)) {
        out <- if (is.matrix(value)) value else as.matrix(.amatrix_host_arg(value))
        if (!is.double(out)) {
          storage.mode(out) <- "double"
        }
        return(out)
      }
    } else {
      .amatrix_cleanup_temp_resident(list(lhs), backend_name)
    }
  }

  old_defer <- getOption("amatrix.defer_host", FALSE)
  options(amatrix.defer_host = FALSE)
  on.exit(options(amatrix.defer_host = old_defer), add = TRUE)

  resident <- switch(
    op,
    matmul = .amatrix_try_resident_matmul(lhs_bound, y, backend_name),
    crossprod = .amatrix_try_resident_crossprod(lhs_bound, y, backend_name),
    tcrossprod = .amatrix_try_resident_tcrossprod(lhs_bound, y, backend_name)
  )
  if (is.null(resident)) {
    return(NULL)
  }

  value <- resident$value
  if (is.null(value) && !isTRUE(resident$host_only) &&
      !is.null(resident$backend) && !is.null(resident$resident_key)) {
    backend <- tryCatch(.amatrix_get_backend(resident$backend), error = function(e) NULL)
    if (!is.null(backend) && is.function(backend$resident_materialize)) {
      value <- tryCatch(backend$resident_materialize(resident$resident_key), error = function(e) NULL)
    }
  }
  .amatrix_product_plan_drop_output(resident)

  if (is.null(value)) {
    return(NULL)
  }
  out <- if (is.matrix(value)) value else as.matrix(.amatrix_host_arg(value))
  if (!is.double(out)) {
    storage.mode(out) <- "double"
  }
  out
}

#' Compile a reusable matrix-product plan
#'
#' Prepares a fixed left operand for repeated products, choosing and binding a
#' resident accelerator backend when beneficial. The returned object is a
#' callable function, so repeated right-hand sides can be applied without
#' rethinking backend selection each time.
#'
#' @param x Fixed left operand.
#' @param op Product primitive: \code{"matmul"}, \code{"crossprod"}, or
#'   \code{"tcrossprod"}.
#' @param backend Backend name or \code{"auto"}.
#' @param precision Precision to use when wrapping base matrices.
#' @param policy Policy to use when wrapping base matrices.
#' @return A callable \code{am_product_plan}. Calling the plan with
#'   \code{materialize = "matrix"} returns a base matrix directly, which is
#'   useful for internal algorithms that immediately need host matrix data.
#' @export
amatrix_compile_product <- function(
  x,
  op = c("matmul", "crossprod", "tcrossprod"),
  backend = "auto",
  precision = amatrix_default_precision(),
  policy = amatrix_default_policy()
) {
  precision_missing <- missing(precision)
  op <- match.arg(op)
  precision <- .amatrix_resolve_backend_precision(
    backend,
    precision,
    precision_missing = precision_missing
  )

  lhs <- .amatrix_prepare_binary_operand(
    x,
    preferred_backend = backend,
    precision = precision,
    policy = policy
  )

  if (!inherits(lhs, "aMatrix")) {
    stop("x must be matrix-like", call. = FALSE)
  }

  live_backend_before <- .amatrix_live_resident_backend(lhs)
  backend_name <- if (!is.null(backend) && !identical(backend, "auto")) {
    backend
  } else {
    amatrix_resident_backend_for(lhs, op = op)
  }
  if (!is.null(backend_name) && !identical(backend_name, "cpu")) {
    .amatrix_check_backend_precision(backend_name, lhs@precision)
  }

  lhs_bound <- if (!is.null(backend_name) && !identical(backend_name, "cpu")) {
    amatrix_bind_resident(lhs, backend = backend_name, op = op)
  } else {
    lhs
  }
  lhs_host_cache <- if (inherits(lhs_bound, "adgCMatrix")) {
    amatrix_materialize_host(lhs_bound)
  } else {
    NULL
  }
  live_backend_after <- .amatrix_live_resident_backend(lhs_bound)
  owned_resident <- !is.null(backend_name) &&
    !identical(backend_name, "cpu") &&
    identical(live_backend_after, backend_name) &&
    !identical(live_backend_before, backend_name)

  plan_fun <- function(y, materialize = c("amatrix", "matrix")) {
    materialize <- match.arg(materialize)
    if (identical(materialize, "matrix")) {
      matrix_value <- .amatrix_product_plan_matrix_result(
        lhs_bound,
        y,
        op,
        backend_name,
        lhs_host = lhs_host_cache
      )
      if (!is.null(matrix_value)) {
        return(matrix_value)
      }
    }

    # Hot sparse-lhs plans can call the public wrappers directly. They already
    # know how to temp-upload and clean up a varying dense rhs while reusing the
    # bound sparse lhs resident key.
    value <- if (inherits(lhs_bound, "adgCMatrix")) {
      switch(
        op,
        matmul = matmul(lhs_bound, y),
        crossprod = am_crossprod(lhs_bound, y),
        tcrossprod = am_tcrossprod(lhs_bound, y)
      )
    } else {
      # Dense lhs with sparse rhs still needs op-aware lowering and backend
      # alignment, so keep using the higher-level preparer in that case.
      prep <- amatrix_prepare_operands(
        lhs_bound,
        y,
        op = op,
        backend = if (is.null(backend_name)) "auto" else backend_name,
        precision = lhs_bound@precision,
        policy = lhs_bound@policy
      )

      switch(
        op,
        matmul = matmul(prep$x, prep$y),
        crossprod = am_crossprod(prep$x, prep$y),
        tcrossprod = am_tcrossprod(prep$x, prep$y)
      )
    }

    if (identical(materialize, "matrix")) {
      return(as.matrix(value))
    }
    value
  }

  meta <- list(
    op = op,
    backend = backend_name,
    lhs_class = class(lhs_bound)[[1L]],
    precision = lhs_bound@precision,
    policy = lhs_bound@policy,
    resident = !is.null(.amatrix_live_resident_backend(lhs_bound)),
    owned_resident = owned_resident
  )

  attr(plan_fun, "amatrix_plan_meta") <- meta
  class(plan_fun) <- c("am_product_plan", "function")
  plan_fun
}

#' @export
print.am_product_plan <- function(x, ...) {
  meta <- .amatrix_product_plan_meta(x)
  backend_label <- if (is.null(meta$backend)) "cpu/host" else meta$backend
  resident_label <- if (isTRUE(meta$resident)) "resident" else "host"

  cat(sprintf(
    "am_product_plan [%s | lhs=%s | backend=%s | %s | precision=%s]\n",
    meta$op,
    meta$lhs_class,
    backend_label,
    resident_label,
    meta$precision
  ))
  invisible(x)
}

#' @export
predict.am_product_plan <- function(object, newdata, ...) {
  object(newdata)
}
