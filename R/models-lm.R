.amatrix_model_sparse_or_dense_arg <- function(x, template = NULL) {
  if (inherits(x, "adgCMatrix")) return(x)
  if (inherits(x, "dgCMatrix") || inherits(x, "sparseMatrix")) {
    preferred_backend <- if (is.null(template)) "cpu" else template@preferred_backend
    policy <- if (is.null(template)) amatrix_default_policy() else template@policy
    precision <- if (is.null(template)) amatrix_default_precision() else template@precision
    return(new_adgCMatrix(x, preferred_backend = preferred_backend,
                          policy = policy, precision = precision))
  }
  .amatrix_model_dense_arg(x, template)
}

.amatrix_model_dense_arg <- function(x, template = NULL) {
  if (inherits(x, "adgeMatrix")) {
    return(x)
  }

  if (inherits(x, "aMatrix")) {
    return(as_adgeMatrix(
      as.matrix(amatrix_materialize_host(x)),
      preferred_backend = x@preferred_backend,
      policy = x@policy,
      precision = x@precision
    ))
  }

  if (inherits(x, "Matrix") || is.matrix(x)) {
    preferred_backend <- if (is.null(template)) "cpu" else template@preferred_backend
    policy <- if (is.null(template)) amatrix_default_policy() else template@policy
    precision <- if (is.null(template)) amatrix_default_precision() else template@precision
    return(as_adgeMatrix(
      as.matrix(x),
      preferred_backend = preferred_backend,
      policy = policy,
      precision = precision
    ))
  }

  stop("x must be a dense matrix-like object")
}

.amatrix_model_response_arg <- function(y, template) {
  if (is.null(dim(y))) {
    y <- matrix(y, ncol = 1L)
  }

  .amatrix_model_dense_arg(y, template = template)
}

.amatrix_model_response_host <- function(y) {
  value <- .amatrix_host_arg(y)
  if (is.null(dim(value))) {
    return(matrix(value, ncol = 1L))
  }
  as.matrix(value)
}

.amatrix_response_layout <- function(Y) {
  dims <- dim(Y)

  if (is.null(dims)) {
    return(list(
      y_matrix = matrix(Y, ncol = 1L),
      observations = length(Y),
      response_dims = 1L,
      responses = 1L,
      is_array = FALSE
    ))
  }

  if (length(dims) == 2L) {
    return(list(
      y_matrix = Y,
      observations = dims[[1]],
      response_dims = dims[[2]],
      responses = dims[[2]],
      is_array = FALSE
    ))
  }

  observations <- dims[[1]]
  response_dims <- dims[-1]
  responses <- prod(response_dims)
  y_matrix <- matrix(Y, nrow = observations, ncol = responses)

  list(
    y_matrix = y_matrix,
    observations = observations,
    response_dims = response_dims,
    responses = responses,
    is_array = TRUE
  )
}

.amatrix_response_summary <- function(Y) {
  dims <- dim(Y)

  if (is.null(dims)) {
    return(list(observations = as.integer(length(Y)), responses = 1L))
  }

  if (length(dims) == 1L) {
    return(list(observations = as.integer(dims[[1]]), responses = 1L))
  }

  list(observations = as.integer(dims[[1]]), responses = as.integer(prod(dims[-1])))
}

.amatrix_restore_response_layout <- function(value, response_dims) {
  if (is.null(value) || is.null(response_dims)) {
    return(value)
  }

  value_matrix <- as.matrix(amatrix_materialize_host(value))
  array(value_matrix, dim = c(nrow(value_matrix), response_dims))
}

.amatrix_model_design_arg <- function(X, intercept = FALSE) {
  X_arg <- .amatrix_model_dense_arg(X)

  if (!isTRUE(intercept)) {
    return(X_arg)
  }

  intercept_col <- matrix(1, nrow = nrow(X_arg), ncol = 1L)
  X_host <- cbind(intercept_col, as.matrix(amatrix_materialize_host(X_arg)))
  adgeMatrix(
    X_host,
    preferred_backend = X_arg@preferred_backend,
    policy = X_arg@policy,
    precision = X_arg@precision
  )
}

.amatrix_dense_first_column <- function(X_arg) {
  stopifnot(inherits(X_arg, "adgeMatrix"))

  nr <- nrow(X_arg)
  if (nr < 1L) {
    return(numeric())
  }

  fenv <- X_arg@finalizer_env
  if (isTRUE(fenv$host_deferred)) {
    if (!is.null(fenv$host_x)) {
      return(as.matrix(fenv$host_x)[, 1L, drop = TRUE])
    }
  } else {
    return(X_arg@x[seq_len(nr)])
  }

  backend_name <- .amatrix_live_resident_backend(X_arg)
  if (!is.null(backend_name)) {
    selector <- matrix(0, nrow = ncol(X_arg), ncol = 1L)
    selector[1L, 1L] <- 1
    resident <- .amatrix_try_resident_matmul(X_arg, selector, backend_name)
    if (!is.null(resident)) {
      backend <- .amatrix_get_backend(backend_name)
      on.exit({
        if (!is.null(resident$resident_key) && isTRUE(backend$resident_has(resident$resident_key))) {
          backend$resident_drop(resident$resident_key)
        }
      }, add = TRUE)
      value <- resident$value
      if (is.null(value)) {
        value <- backend$resident_materialize(resident$resident_key)
      }
      return(drop(as.matrix(value)))
    }
  }

  as.matrix(amatrix_materialize_host(X_arg))[, 1L, drop = TRUE]
}

.amatrix_has_explicit_intercept_column <- function(X_arg) {
  stopifnot(inherits(X_arg, "adgeMatrix"))

  if (ncol(X_arg) < 1L) {
    return(FALSE)
  }

  first_col <- .amatrix_dense_first_column(X_arg)
  isTRUE(all(first_col == 1))
}

.amatrix_lm_cache_key <- function(X_arg, extra = NULL, object_id = NULL) {
  stopifnot(inherits(X_arg, "adgeMatrix"))
  paste(
    "lm",
    if (is.null(object_id)) X_arg@object_id else object_id,
    X_arg@preferred_backend,
    X_arg@policy,
    X_arg@precision,
    if (is.null(extra)) "" else extra,
    sep = "|"
  )
}

.amatrix_qr_cache_signature <- function(X_arg) {
  stopifnot(inherits(X_arg, "adgeMatrix"))

  backend <- X_arg@preferred_backend
  dims <- dim(X_arg)

  if (identical(backend, "mlx") && requireNamespace("amatrix.mlx", quietly = TRUE)) {
    sig_fun <- get("amatrix_mlx_qr_cache_signature", envir = asNamespace("amatrix.mlx"), inherits = FALSE)
    return(sig_fun(dims))
  }

  paste("qr", backend, X_arg@precision, sep = ":")
}

.amatrix_lm_cache_get <- function(cache_key) {
  .amatrix_cache_get(cache_key)
}

.amatrix_lm_cache_set <- function(cache_key, value) {
  .amatrix_cache_set(cache_key, value)
}

.amatrix_plan_prefers_cpu_solve <- function(lhs, rhs = NULL) {
  identical(amatrix_backend_plan(lhs, "solve", y = rhs)$chosen, "cpu")
}

.amatrix_host_solve_rewrap <- function(lhs, rhs = NULL, template = lhs) {
  lhs_host <- as.matrix(amatrix_materialize_host(lhs))
  if (is.null(rhs)) {
    return(.amatrix_rewrap_like(template, base::solve(lhs_host)))
  }

  rhs_host <- as.matrix(amatrix_materialize_host(rhs))
  .amatrix_rewrap_like(template, base::solve(lhs_host, rhs_host))
}

.amatrix_ridge_penalized_xtx <- function(XtX, lambda, penalize_intercept = FALSE, has_intercept = FALSE) {
  stopifnot(inherits(XtX, "adgeMatrix"))

  xtx_host <- as.matrix(amatrix_materialize_host(XtX))
  diag_len <- min(nrow(xtx_host), ncol(xtx_host))
  diag_index <- cbind(seq_len(diag_len), seq_len(diag_len))
  xtx_host[diag_index] <- xtx_host[diag_index] + as.double(lambda)

  if (isTRUE(has_intercept) && !isTRUE(penalize_intercept) && diag_len >= 1L) {
    xtx_host[1L, 1L] <- xtx_host[1L, 1L] - as.double(lambda)
  }

  .amatrix_rewrap_like(XtX, xtx_host)
}

.amatrix_lm_cache_value <- function(X_arg, cache = TRUE, need_xtx = TRUE, need_qr = FALSE, cache_key = NULL) {
  stopifnot(inherits(X_arg, "adgeMatrix"))

  cache_key <- if (is.null(cache_key)) .amatrix_lm_cache_key(X_arg) else cache_key
  cached <- if (isTRUE(cache)) .amatrix_lm_cache_get(cache_key) else NULL
  cache_reused <- FALSE
  if (is.null(cached)) {
    cached <- list()
  } else {
    cache_reused <- TRUE
  }

  if (isTRUE(need_xtx) && is.null(cached$xtx)) {
    cached$xtx <- am_crossprod(X_arg)
  }

  if (isTRUE(need_qr) && is.null(cached$qr)) {
    cached$qr <- am_qr(X_arg)
  }

  if (is.null(cached$rank)) {
    if (!is.null(cached$qr)) {
      cached$rank <- .amatrix_qr_rank(cached$qr)
    } else {
      qr_x <- base::qr(as.matrix(amatrix_materialize_host(X_arg)))
      cached$rank <- qr_x$rank
    }
  }

  if (isTRUE(cache)) {
    .amatrix_lm_cache_set(cache_key, cached)
  }

  c(cached, list(cache_key = cache_key, cache_reused = cache_reused))
}

.amatrix_model_outputs <- function(X_arg, Y_arg, coefficients, include_fitted = TRUE, include_residuals = TRUE) {
  fitted_values <- NULL
  if (isTRUE(include_fitted) || isTRUE(include_residuals)) {
    fitted_values <- matmul(X_arg, coefficients)
  }

  residuals <- NULL
  if (isTRUE(include_residuals)) {
    residuals <- ewise("-", Y_arg, fitted_values)
  }

  list(fitted_values = fitted_values, residuals = residuals)
}

.amatrix_qr_fit_metadata <- function(qr_fit = NULL) {
  if (is.null(qr_fit)) {
    return(list(
      qr_representation = NULL,
      qr_helper_path = NULL,
      qr_compact_factor_available = FALSE,
      qr_compact_factor_source = NULL
    ))
  }

  list(
    qr_representation = .amatrix_qr_kind(qr_fit),
    qr_helper_path = .amatrix_qr_helper_path(qr_fit),
    qr_compact_factor_available = isTRUE(.amatrix_qr_compact_available(qr_fit)),
    qr_compact_factor_source = qr_fit$state$factor_source
  )
}

.amatrix_profile_many_lm_qr <- function(X_arg, Y, cache = TRUE) {
  stopifnot(inherits(X_arg, "adgeMatrix"))

  Y_host <- .amatrix_model_response_host(Y)
  cache_key <- .amatrix_lm_cache_key(X_arg, extra = .amatrix_qr_cache_signature(X_arg))

  t0 <- proc.time()[["elapsed"]]
  cache_value <- .amatrix_lm_cache_value(X_arg, cache = cache, need_xtx = FALSE, need_qr = TRUE, cache_key = cache_key)
  t_cache <- proc.time()[["elapsed"]] - t0

  qr_fit <- cache_value$qr
  t0 <- proc.time()[["elapsed"]]
  coef_host <- .amatrix_qr_solve_value(qr_fit, b = Y_host)
  t_solve <- proc.time()[["elapsed"]] - t0

  t0 <- proc.time()[["elapsed"]]
  coefficients <- .amatrix_rewrap_like(X_arg, coef_host)
  fit <- .amatrix_make_many_lm_fit(
    list(
      coefficients = coefficients,
      fitted.values = NULL,
      residuals = NULL,
      xtx = NULL,
      xty = NULL,
      method = "qr",
      rank = cache_value$rank,
      cache = isTRUE(cache),
      cache_key = cache_value$cache_key,
      cache_reused = cache_value$cache_reused,
      qr_representation = .amatrix_qr_kind(qr_fit),
      qr_helper_path = .amatrix_qr_helper_path(qr_fit),
      qr_compact_factor_available = isTRUE(.amatrix_qr_compact_available(qr_fit)),
      qr_compact_factor_source = qr_fit$state$factor_source
    ),
    X_arg,
    Y,
    weights = NULL,
    call = NULL
  )
  t_wrap <- proc.time()[["elapsed"]] - t0

  list(
    fit = fit,
    timings = c(cache = t_cache, solve = t_solve, assemble = t_wrap),
    cache_key = cache_value$cache_key,
    cache_reused = cache_value$cache_reused,
    qr_representation = .amatrix_qr_kind(qr_fit),
    qr_helper_path = .amatrix_qr_helper_path(qr_fit),
    qr_compact_factor_source = qr_fit$state$factor_source
  )
}

.amatrix_many_lm_qr_hot <- function(X_arg, Y, cache = TRUE, include_fitted = FALSE, include_residuals = FALSE) {
  stopifnot(inherits(X_arg, "adgeMatrix"))

  Y_host <- .amatrix_model_response_host(Y)
  cache_key <- .amatrix_lm_cache_key(X_arg, extra = .amatrix_qr_cache_signature(X_arg))
  cache_value <- .amatrix_lm_cache_value(X_arg, cache = cache, need_xtx = FALSE, need_qr = TRUE, cache_key = cache_key)
  qr_fit <- cache_value$qr
  coefficients <- .amatrix_rewrap_like(X_arg, .amatrix_qr_solve_value(qr_fit, b = Y_host))
  fitted_values <- if (isTRUE(include_fitted) || isTRUE(include_residuals)) .amatrix_rewrap_like(X_arg, .amatrix_qr_fitted_value(qr_fit, Y_host)) else NULL
  residuals <- if (isTRUE(include_residuals)) .amatrix_rewrap_like(X_arg, .amatrix_qr_resid_value(qr_fit, Y_host)) else NULL

  .amatrix_make_many_lm_fit(
    list(
      coefficients = coefficients,
      fitted.values = fitted_values,
      residuals = residuals,
      xtx = NULL,
      xty = NULL,
      method = "qr",
      rank = cache_value$rank,
      cache = isTRUE(cache),
      cache_key = cache_value$cache_key,
      cache_reused = cache_value$cache_reused,
      qr_representation = .amatrix_qr_kind(qr_fit),
      qr_helper_path = .amatrix_qr_helper_path(qr_fit),
      qr_compact_factor_available = isTRUE(.amatrix_qr_compact_available(qr_fit)),
      qr_compact_factor_source = qr_fit$state$factor_source
    ),
    X_arg,
    Y,
    weights = NULL,
    call = match.call()
  )
}

.amatrix_lm_core <- function(
    X_arg,
    Y,
    cache = TRUE,
    include_fitted = TRUE,
    include_residuals = TRUE,
    method = c("normal", "qr")) {
  method <- match.arg(method)

  if (identical(method, "qr")) {
    Y_host <- .amatrix_model_response_host(Y)
    cache_value <- .amatrix_lm_cache_value(
      X_arg,
      cache = cache,
      need_xtx = FALSE,
      need_qr = TRUE,
      cache_key = .amatrix_lm_cache_key(X_arg, extra = .amatrix_qr_cache_signature(X_arg))
    )
    qr_fit <- cache_value$qr
    coefficients <- .amatrix_rewrap_like(X_arg, .amatrix_qr_solve_value(qr_fit, b = Y_host))
    fitted_values <- if (isTRUE(include_fitted) || isTRUE(include_residuals)) .amatrix_rewrap_like(X_arg, .amatrix_qr_fitted_value(qr_fit, Y_host)) else NULL
    residuals <- if (isTRUE(include_residuals)) .amatrix_rewrap_like(X_arg, .amatrix_qr_resid_value(qr_fit, Y_host)) else NULL
    qr_meta <- .amatrix_qr_fit_metadata(qr_fit)
    XtX <- NULL
    XtY <- NULL
  } else if (inherits(X_arg, "adgCMatrix")) {
    # Sparse "normal" path: bypass cache (which requires adgeMatrix)
    Y_arg <- .amatrix_model_response_arg(Y, template = X_arg)
    XtX <- am_crossprod(X_arg)
    XtY <- am_crossprod(X_arg, Y_arg)
    coefficients <- am_solve(XtX, XtY)
    outputs <- .amatrix_model_outputs(
      X_arg,
      Y_arg,
      coefficients,
      include_fitted = include_fitted,
      include_residuals = include_residuals
    )
    fitted_values <- outputs$fitted_values
    residuals <- outputs$residuals
    qr_meta <- .amatrix_qr_fit_metadata(NULL)
    rank_qr <- base::qr(as.matrix(amatrix_materialize_host(X_arg)))
    cache_value <- list(rank = rank_qr$rank, cache_key = NA_character_,
                        cache_reused = FALSE)
  } else {
    Y_arg <- .amatrix_model_response_arg(Y, template = X_arg)
    cache_value <- .amatrix_lm_cache_value(X_arg, cache = cache, need_xtx = TRUE, need_qr = FALSE)
    XtX <- cache_value$xtx
    XtY <- am_crossprod(X_arg, Y_arg)
    coefficients <- if (.amatrix_plan_prefers_cpu_solve(XtX, XtY)) {
      .amatrix_host_solve_rewrap(XtX, XtY, template = X_arg)
    } else {
      am_solve(XtX, XtY)
    }
    outputs <- .amatrix_model_outputs(
      X_arg,
      Y_arg,
      coefficients,
      include_fitted = include_fitted,
      include_residuals = include_residuals
    )
    fitted_values <- outputs$fitted_values
    residuals <- outputs$residuals
    qr_meta <- .amatrix_qr_fit_metadata(NULL)
  }

  list(
    coefficients = coefficients,
    fitted.values = fitted_values,
    residuals = residuals,
    xtx = XtX,
    xty = XtY,
    method = method,
    rank = cache_value$rank,
    cache = isTRUE(cache),
    cache_key = cache_value$cache_key,
    cache_reused = cache_value$cache_reused,
    qr_representation = qr_meta$qr_representation,
    qr_helper_path = qr_meta$qr_helper_path,
    qr_compact_factor_available = qr_meta$qr_compact_factor_available,
    qr_compact_factor_source = qr_meta$qr_compact_factor_source
  )
}

.amatrix_make_lm_fit <- function(core, X_arg, intercept = FALSE, call = NULL) {
  structure(
    c(
      core,
      list(
        df.residual = nrow(X_arg) - core$rank,
        intercept = isTRUE(intercept),
        precision = X_arg@precision,
        backend = X_arg@preferred_backend,
        call = call
      )
    ),
    class = "lm_fit"
  )
}

.amatrix_make_many_lm_fit <- function(core, X_arg, Y, weights = NULL, call = NULL) {
  y_info <- .amatrix_response_summary(Y)
  metrics <- .amatrix_many_lm_metrics(core$residuals, nrow(X_arg) - core$rank, weights = weights)

  structure(
    c(
      core,
      list(
        weights = if (is.null(weights)) NULL else as.double(weights),
        responses = y_info$responses,
        observations = y_info$observations,
        rss = metrics$rss,
        sigma2 = metrics$sigma2,
        df.residual = nrow(X_arg) - core$rank,
        precision = X_arg@precision,
        backend = X_arg@preferred_backend,
        call = call
      )
    ),
    class = "am_many_lm_fit"
  )
}

.amatrix_penalty_matrix <- function(X_arg, lambda, penalize_intercept = FALSE, has_intercept = FALSE) {
  stopifnot(inherits(X_arg, "adgeMatrix"))
  penalty <- diag(as.double(lambda), ncol(X_arg))
  if (isTRUE(has_intercept) && !isTRUE(penalize_intercept) && ncol(X_arg) >= 1L) {
    penalty[1, 1] <- 0
  }
  adgeMatrix(
    penalty,
    preferred_backend = X_arg@preferred_backend,
    policy = X_arg@policy,
    precision = X_arg@precision
  )
}

.amatrix_weights_signature <- function(weights) {
  paste(
    "w",
    length(weights),
    sprintf("%.8f", sum(weights)),
    sprintf("%.8f", sum(weights^2)),
    sep = ":"
  )
}

.amatrix_apply_row_weights <- function(X_arg, weights) {
  stopifnot(inherits(X_arg, "adgeMatrix"))
  sqrt_w <- sqrt(as.double(weights))

  # GPU-resident path: scale rows on GPU via broadcast_ewise, no host round-trip
  backend_name <- .amatrix_live_resident_backend(X_arg)
  if (!is.null(backend_name)) {
    backend <- .amatrix_get_backend(backend_name)
    if (.amatrix_backend_supports_resident_op(backend, "broadcast_ewise")) {
      lhs <- .amatrix_prepare_resident_arg(X_arg, backend_name)
      if (!is.null(lhs)) {
        out_key <- .amatrix_next_resident_key(backend_name)
        result <- tryCatch({
          val <- backend$broadcast_ewise_resident(lhs$key, sqrt_w, 1L, "*", out_key)
          .amatrix_cleanup_temp_resident(list(lhs), backend_name)
          val
        }, error = function(e) NULL)
        if (!is.null(result)) {
          value <- .amatrix_rewrap_like(X_arg, result)
          return(.amatrix_bind_resident(value, backend_name, out_key))
        }
        try(backend$resident_drop(out_key), silent = TRUE)
        .amatrix_cleanup_temp_resident(list(lhs), backend_name)
      }
    }
  }

  # CPU fallback
  x_host <- as.matrix(amatrix_materialize_host(X_arg))
  weighted <- x_host * sqrt_w
  adgeMatrix(
    weighted,
    preferred_backend = X_arg@preferred_backend,
    policy = X_arg@policy,
    precision = X_arg@precision
  )
}

.amatrix_validate_weights <- function(weights, n_obs) {
  if (is.null(weights)) {
    return(NULL)
  }

  if (!is.numeric(weights) || anyNA(weights) || any(weights < 0)) {
    stop("weights must be a numeric vector of non-missing non-negative values")
  }

  if (length(weights) != n_obs) {
    stop("weights must have length equal to nrow(X)")
  }

  if (sum(weights) <= 0) {
    stop("weights must have positive total weight")
  }

  as.double(weights)
}

.amatrix_centered_design <- function(X_arg, center = TRUE, weights = NULL) {
  stopifnot(inherits(X_arg, "adgeMatrix"))

  x_host <- as.matrix(amatrix_materialize_host(X_arg))

  if (!isTRUE(center)) {
    return(X_arg)
  }

  centers <- if (is.null(weights)) {
    colMeans(x_host)
  } else {
    colSums(x_host * as.double(weights)) / sum(weights)
  }

  centered <- sweep(x_host, 2L, centers, FUN = "-")
  adgeMatrix(
    centered,
    preferred_backend = X_arg@preferred_backend,
    policy = X_arg@policy,
    precision = X_arg@precision
  )
}

.amatrix_validate_block_size <- function(block_size, n_cols) {
  if (is.null(block_size)) {
    return(NULL)
  }

  if (!is.numeric(block_size) || length(block_size) != 1L || is.na(block_size)) {
    stop("block_size must be NULL or a single positive integer")
  }

  block_size <- as.integer(block_size)
  if (block_size <= 0L) {
    stop("block_size must be NULL or a single positive integer")
  }

  min(block_size, as.integer(n_cols))
}

.amatrix_dense_like <- function(x_host, template) {
  adgeMatrix(
    x_host,
    preferred_backend = template@preferred_backend,
    policy = template@policy,
    precision = template@precision
  )
}

.amatrix_covariance_blockwise <- function(X_centered, block_size, weights = NULL) {
  stopifnot(inherits(X_centered, "adgeMatrix"))

  x_host <- as.matrix(amatrix_materialize_host(X_centered))
  n_cols <- ncol(x_host)
  block_size <- .amatrix_validate_block_size(block_size, n_cols)
  if (is.null(block_size) || block_size >= n_cols) {
    return(NULL)
  }

  if (!is.null(weights)) {
    x_host <- x_host * sqrt(as.double(weights))
  }

  cov_host <- matrix(0, nrow = n_cols, ncol = n_cols)
  starts <- seq.int(1L, n_cols, by = block_size)

  for (i_start in starts) {
    i_end <- min(i_start + block_size - 1L, n_cols)
    i_idx <- i_start:i_end
    Xi <- .amatrix_dense_like(x_host[, i_idx, drop = FALSE], X_centered)

    for (j_start in starts[starts >= i_start]) {
      j_end <- min(j_start + block_size - 1L, n_cols)
      j_idx <- j_start:j_end
      Xj <- .amatrix_dense_like(x_host[, j_idx, drop = FALSE], X_centered)
      block <- as.matrix(amatrix_materialize_host(am_crossprod(Xi, Xj)))

      cov_host[i_idx, j_idx] <- block
      if (j_start != i_start) {
        cov_host[j_idx, i_idx] <- t(block)
      }
    }
  }

  .amatrix_dense_like(cov_host, X_centered)
}

covariance <- function(X, center = TRUE, sample = TRUE, weights = NULL, block_size = NULL) {
  # Sparse path: avoid densifying large sparse matrices.
  # Formula: cov(X) = (X^TX - n * mu * mu^T) / (n-1)
  # X^TX via Matrix::crossprod (sparse BLAS, O(NNZ*p)); mu via Matrix::colMeans (O(NNZ)).
  if (is.null(weights) && is.null(block_size) &&
      (inherits(X, "sparseMatrix") || inherits(X, "adgCMatrix"))) {
    n <- nrow(X)
    denom <- if (isTRUE(sample)) n - 1L else n
    if (denom <= 0L) stop("effective denominator must be positive")
    XtX <- as.matrix(Matrix::crossprod(X))   # p x p dense result
    if (isTRUE(center)) {
      mu  <- Matrix::colMeans(X)             # length-p vector
      XtX <- XtX - n * tcrossprod(mu)
    }
    return(as_adgeMatrix(XtX / denom))
  }

  X_arg <- .amatrix_model_dense_arg(X)
  weights <- .amatrix_validate_weights(weights, nrow(X_arg))

  denom <- if (is.null(weights)) {
    if (isTRUE(sample)) nrow(X_arg) - 1 else nrow(X_arg)
  } else {
    if (isTRUE(sample)) sum(weights) - 1 else sum(weights)
  }

  if (denom <= 0) {
    stop("effective denominator must be positive")
  }

  # Fused GPU path: center + XtX + scale in one MLX lazy graph.
  # Only available when weights and block_size are NULL (unweighted, non-blockwise).
  if (is.null(weights) && is.null(block_size)) {
    bk_name <- X_arg@preferred_backend
    bk <- tryCatch(.amatrix_get_backend(bk_name), error = function(e) NULL)
    if (!is.null(bk) && is.function(bk[["covariance"]]) && isTRUE(bk$supports("covariance", X_arg))) {
      x_mat <- as.matrix(amatrix_materialize_host(X_arg))
      cov_mat <- bk$covariance(x_mat, center = center, denom = denom)
      return(adgeMatrix(
        cov_mat,
        preferred_backend = X_arg@preferred_backend,
        policy = X_arg@policy,
        precision = X_arg@precision
      ))
    }
  }

  # Fallback: CPU centering + GPU XtX + CPU scaling.
  X_centered <- .amatrix_centered_design(X_arg, center = center, weights = weights)
  gram_arg <- .amatrix_covariance_blockwise(X_centered, block_size = block_size, weights = weights)
  if (is.null(gram_arg)) {
    gram_arg <- if (is.null(weights)) {
      am_crossprod(X_centered)
    } else {
      am_crossprod(.amatrix_apply_row_weights(X_centered, weights))
    }
  }

  cov_host <- as.matrix(amatrix_materialize_host(gram_arg)) / denom
  adgeMatrix(
    cov_host,
    preferred_backend = X_arg@preferred_backend,
    policy = X_arg@policy,
    precision = X_arg@precision
  )
}

correlation <- function(X, center = TRUE, weights = NULL, block_size = NULL) {
  cov_arg <- covariance(X, center = center, sample = TRUE, weights = weights, block_size = block_size)
  cov_host <- as.matrix(amatrix_materialize_host(cov_arg))
  sds <- sqrt(diag(cov_host))
  scale_mat <- outer(sds, sds)
  cor_host <- cov_host / scale_mat
  cor_host[!is.finite(cor_host)] <- NA_real_
  diag(cor_host) <- 1

  adgeMatrix(
    cor_host,
    preferred_backend = cov_arg@preferred_backend,
    policy = cov_arg@policy,
    precision = cov_arg@precision
  )
}

lm_fit <- function(
  X,
  Y,
  intercept = FALSE,
  include_fitted = TRUE,
  include_residuals = TRUE,
  cache = TRUE,
  method = c("normal", "qr")
) {
  method <- match.arg(method)

  # For "normal" method, preserve sparsity if X is sparse (no intercept —

  # intercept column forces dense anyway).
  if (identical(method, "normal") && !isTRUE(intercept)) {
    X_arg <- .amatrix_model_sparse_or_dense_arg(X)
  } else {
    X_arg <- .amatrix_model_design_arg(X, intercept = intercept)
  }

  core <- .amatrix_lm_core(
    X_arg,
    Y,
    cache = cache,
    include_fitted = include_fitted,
    include_residuals = include_residuals,
    method = method
  )

  .amatrix_make_lm_fit(core, X_arg, intercept = intercept, call = match.call())
}

ridge_fit <- function(
  X,
  Y,
  lambda,
  intercept = FALSE,
  penalize_intercept = FALSE,
  include_fitted = TRUE,
  include_residuals = TRUE,
  cache = TRUE
) {
  if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda) || lambda < 0) {
    stop("lambda must be a single non-negative numeric value")
  }

  X_arg <- .amatrix_model_design_arg(X, intercept = intercept)
  Y_arg <- .amatrix_model_response_arg(Y, template = X_arg)

  cache_value <- .amatrix_lm_cache_value(X_arg, cache = cache)
  XtX <- cache_value$xtx
  XtY <- am_crossprod(X_arg, Y_arg)
  has_intercept_col <- isTRUE(intercept) || .amatrix_has_explicit_intercept_column(X_arg)
  if (.amatrix_plan_prefers_cpu_solve(XtX, XtY)) {
    penalized_xtx <- .amatrix_ridge_penalized_xtx(
      XtX,
      lambda,
      penalize_intercept = penalize_intercept,
      has_intercept = has_intercept_col
    )
    coefficients <- .amatrix_host_solve_rewrap(penalized_xtx, XtY, template = X_arg)
  } else {
    penalty <- .amatrix_penalty_matrix(
      X_arg,
      lambda,
      penalize_intercept = penalize_intercept,
      has_intercept = has_intercept_col
    )
    penalized_xtx <- ewise("+", XtX, penalty)
    coefficients <- am_solve(penalized_xtx, XtY)
  }
  outputs <- .amatrix_model_outputs(
    X_arg,
    Y_arg,
    coefficients,
    include_fitted = include_fitted,
    include_residuals = include_residuals
  )

  rank <- cache_value$rank

  structure(
    list(
      coefficients = coefficients,
      fitted.values = outputs$fitted_values,
      residuals = outputs$residuals,
      xtx = XtX,
      penalized_xtx = penalized_xtx,
      xty = XtY,
      lambda = as.double(lambda),
      rank = rank,
      df.residual = nrow(X_arg) - rank,
      intercept = isTRUE(intercept),
      penalize_intercept = isTRUE(penalize_intercept),
      cache = isTRUE(cache),
      cache_key = cache_value$cache_key,
      cache_reused = cache_value$cache_reused,
      precision = X_arg@precision,
      backend = X_arg@preferred_backend,
      call = match.call()
    ),
    class = "ridge_fit"
  )
}

wls_fit <- function(
  X,
  Y,
  weights,
  intercept = FALSE,
  include_fitted = TRUE,
  include_residuals = TRUE,
  cache = TRUE,
  method = c("normal", "qr")
) {
  X_arg <- .amatrix_model_design_arg(X, intercept = intercept)
  Y_arg <- .amatrix_model_response_arg(Y, template = X_arg)
  method <- match.arg(method)
  weights <- .amatrix_validate_weights(weights, nrow(X_arg))

  cache_extra <- .amatrix_weights_signature(weights)
  if (identical(method, "qr")) {
    weighted_X <- .amatrix_apply_row_weights(X_arg, weights)
    cache_extra <- paste(cache_extra, .amatrix_qr_cache_signature(weighted_X), sep = "|")
    cache_key <- .amatrix_lm_cache_key(
      weighted_X,
      extra = cache_extra,
      object_id = X_arg@object_id
    )
    weighted_Y <- .amatrix_model_response_host(Y_arg) * sqrt(as.double(weights))
    cache_value <- .amatrix_lm_cache_value(weighted_X, cache = cache, need_xtx = FALSE, need_qr = TRUE, cache_key = cache_key)
    qr_fit <- cache_value$qr
    coefficients <- .amatrix_rewrap_like(weighted_X, .amatrix_qr_solve_value(qr_fit, b = weighted_Y))
    qr_meta <- .amatrix_qr_fit_metadata(qr_fit)
  } else {
    cache_key <- .amatrix_lm_cache_key(
      X_arg,
      extra = cache_extra,
      object_id = X_arg@object_id
    )
    xtx <- crossprod_weighted(X_arg, weights)
    xty <- xty_weighted(X_arg, weights, Y_arg)
    coefficients <- if (.amatrix_plan_prefers_cpu_solve(xtx, xty)) {
      .amatrix_host_solve_rewrap(xtx, xty, template = X_arg)
    } else {
      am_solve(xtx, xty)
    }
    cache_value <- list(rank = ncol(X_arg), cache_key = cache_key,
                        cache_reused = FALSE)
    qr_meta <- .amatrix_qr_fit_metadata(NULL)
  }

  outputs <- .amatrix_model_outputs(
    X_arg,
    Y_arg,
    coefficients,
    include_fitted = include_fitted,
    include_residuals = include_residuals
  )

  rank <- cache_value$rank

  structure(
    list(
      coefficients = coefficients,
      fitted.values = outputs$fitted_values,
      residuals = outputs$residuals,
      weights = as.double(weights),
      method = method,
      rank = rank,
      df.residual = nrow(X_arg) - rank,
      intercept = isTRUE(intercept),
      cache = isTRUE(cache),
      cache_key = cache_value$cache_key,
      cache_reused = cache_value$cache_reused,
      qr_representation = qr_meta$qr_representation,
      qr_helper_path = qr_meta$qr_helper_path,
      qr_compact_factor_available = qr_meta$qr_compact_factor_available,
      qr_compact_factor_source = qr_meta$qr_compact_factor_source,
      precision = X_arg@precision,
      backend = X_arg@preferred_backend,
      call = match.call()
    ),
    class = "wls_fit"
  )
}

.amatrix_many_lm_metrics <- function(residuals_obj, df_residual, weights = NULL) {
  if (is.null(residuals_obj)) {
    return(list(rss = NULL, sigma2 = NULL))
  }

  residuals_mat <- as.matrix(amatrix_materialize_host(residuals_obj))
  if (is.null(weights)) {
    rss <- colSums(residuals_mat^2)
  } else {
    rss <- colSums(residuals_mat^2 * as.double(weights))
  }
  sigma2 <- rss / df_residual

  list(rss = rss, sigma2 = sigma2)
}

many_lm <- function(
  X,
  Y,
  weights = NULL,
  intercept = FALSE,
  include_fitted = FALSE,
  include_residuals = FALSE,
  cache = TRUE,
  method = c("normal", "qr")
) {
  if (is.null(weights)) {
    X_arg <- .amatrix_model_design_arg(X, intercept = intercept)
    if (identical(match.arg(method), "qr")) {
      return(.amatrix_many_lm_qr_hot(
        X_arg,
        Y,
        cache = cache,
        include_fitted = include_fitted,
        include_residuals = include_residuals
      ))
    }
    core <- .amatrix_lm_core(
      X_arg,
      Y,
      cache = cache,
      include_fitted = include_fitted,
      include_residuals = include_residuals,
      method = method
    )
    return(.amatrix_make_many_lm_fit(core, X_arg, Y, weights = NULL, call = match.call()))
  } else {
    fit <- wls_fit(
      X,
      Y,
      weights = weights,
      intercept = intercept,
      include_fitted = include_fitted,
      include_residuals = include_residuals,
      cache = cache,
      method = method
    )
    y_info <- .amatrix_response_summary(Y)
    metrics <- .amatrix_many_lm_metrics(fit$residuals, fit$df.residual, weights = weights)

    return(structure(
      list(
        coefficients = fit$coefficients,
        fitted.values = fit$fitted.values,
        residuals = fit$residuals,
        weights = as.double(weights),
        responses = y_info$responses,
        observations = y_info$observations,
        rss = metrics$rss,
        sigma2 = metrics$sigma2,
        method = fit$method,
        rank = fit$rank,
        df.residual = fit$df.residual,
        cache = fit$cache,
        cache_reused = fit$cache_reused,
        cache_key = fit$cache_key,
        qr_representation = fit$qr_representation,
        qr_helper_path = fit$qr_helper_path,
        qr_compact_factor_available = fit$qr_compact_factor_available,
        qr_compact_factor_source = fit$qr_compact_factor_source,
        precision = fit$precision,
        backend = fit$backend,
        call = match.call()
      ),
      class = "am_many_lm_fit"
    ))
  }
}

array_lm <- function(
  X,
  Y,
  weights = NULL,
  intercept = FALSE,
  include_fitted = FALSE,
  include_residuals = FALSE,
  cache = TRUE,
  method = c("normal", "qr"),
  restore_array = TRUE
) {
  layout <- .amatrix_response_layout(Y)
  fit <- many_lm(
    X,
    layout$y_matrix,
    weights = weights,
    intercept = intercept,
    include_fitted = include_fitted,
    include_residuals = include_residuals,
    cache = cache,
    method = method
  )

  structure(
    list(
      fit = fit,
      coefficients = fit$coefficients,
      fitted.values = if (isTRUE(restore_array)) .amatrix_restore_response_layout(fit$fitted.values, layout$response_dims) else fit$fitted.values,
      residuals = if (isTRUE(restore_array)) .amatrix_restore_response_layout(fit$residuals, layout$response_dims) else fit$residuals,
      responses = fit$responses,
      observations = fit$observations,
      response_dims = layout$response_dims,
      restore_array = isTRUE(restore_array),
      weights = fit$weights,
      rss = if (!is.null(fit$rss) && length(layout$response_dims) > 1L) array(fit$rss, dim = layout$response_dims) else fit$rss,
      sigma2 = if (!is.null(fit$sigma2) && length(layout$response_dims) > 1L) array(fit$sigma2, dim = layout$response_dims) else fit$sigma2,
      method = fit$method,
      rank = fit$rank,
      df.residual = fit$df.residual,
      cache = fit$cache,
      cache_reused = fit$cache_reused,
      cache_key = fit$cache_key,
      qr_representation = fit$qr_representation,
      qr_helper_path = fit$qr_helper_path,
      qr_compact_factor_available = fit$qr_compact_factor_available,
      qr_compact_factor_source = fit$qr_compact_factor_source,
      precision = fit$precision,
      backend = fit$backend,
      call = match.call()
    ),
    class = "am_array_lm_fit"
  )
}

coef.lm_fit <- function(object, ...) {
  object$coefficients
}

fitted.lm_fit <- function(object, ...) {
  object$fitted.values
}

residuals.lm_fit <- function(object, ...) {
  object$residuals
}

coef.am_many_lm_fit <- function(object, ...) {
  object$coefficients
}

fitted.am_many_lm_fit <- function(object, ...) {
  object$fitted.values
}

residuals.am_many_lm_fit <- function(object, ...) {
  object$residuals
}

coef.wls_fit <- function(object, ...) {
  object$coefficients
}

fitted.wls_fit <- function(object, ...) {
  object$fitted.values
}

residuals.wls_fit <- function(object, ...) {
  object$residuals
}

coef.am_array_lm_fit <- function(object, ...) {
  object$coefficients
}

fitted.am_array_lm_fit <- function(object, ...) {
  object$fitted.values
}

residuals.am_array_lm_fit <- function(object, ...) {
  object$residuals
}

print.lm_fit <- function(x, ...) {
  cat(sprintf(
    "am_lm_fit [backend=%s|precision=%s|method=%s|rank=%d|df.residual=%d|cache=%s",
    x$backend,
    x$precision,
    x$method,
    x$rank,
    x$df.residual,
    if (isTRUE(x$cache_reused)) "hit" else "miss"
  ))
  if (identical(x$method, "qr")) {
    cat(sprintf("|qr=%s|helper=%s|compact=%s", x$qr_representation, .amatrix_qr_or(x$qr_helper_path, "none"), .amatrix_qr_or(x$qr_compact_factor_source, "none")))
  }
  cat("]\n")
  print(coef(x), ...)
  invisible(x)
}

print.am_many_lm_fit <- function(x, ...) {
  cat(sprintf(
    "am_many_lm_fit [backend=%s|precision=%s|method=%s|responses=%d|rank=%d|df.residual=%d|cache=%s",
    x$backend,
    x$precision,
    x$method,
    x$responses,
    x$rank,
    x$df.residual,
    if (isTRUE(x$cache_reused)) "hit" else "miss"
  ))
  if (identical(x$method, "qr")) {
    cat(sprintf("|qr=%s|helper=%s|compact=%s", x$qr_representation, .amatrix_qr_or(x$qr_helper_path, "none"), .amatrix_qr_or(x$qr_compact_factor_source, "none")))
  }
  cat("]\n")
  print(coef(x), ...)
  invisible(x)
}

print.wls_fit <- function(x, ...) {
  cat(sprintf(
    "am_wls_fit [backend=%s|precision=%s|method=%s|rank=%d|df.residual=%d|cache=%s",
    x$backend,
    x$precision,
    x$method,
    x$rank,
    x$df.residual,
    if (isTRUE(x$cache_reused)) "hit" else "miss"
  ))
  if (identical(x$method, "qr")) {
    cat(sprintf("|qr=%s|helper=%s|compact=%s", x$qr_representation, .amatrix_qr_or(x$qr_helper_path, "none"), .amatrix_qr_or(x$qr_compact_factor_source, "none")))
  }
  cat("]\n")
  print(coef(x), ...)
  invisible(x)
}

print.am_array_lm_fit <- function(x, ...) {
  cat(sprintf(
    "am_array_lm_fit [backend=%s|precision=%s|method=%s|responses=%d|response_dims=%s|rank=%d|df.residual=%d|cache=%s",
    x$backend,
    x$precision,
    x$method,
    x$responses,
    paste(x$response_dims, collapse = "x"),
    x$rank,
    x$df.residual,
    if (isTRUE(x$cache_reused)) "hit" else "miss"
  ))
  if (identical(x$method, "qr")) {
    cat(sprintf("|qr=%s|helper=%s|compact=%s", x$qr_representation, .amatrix_qr_or(x$qr_helper_path, "none"), .amatrix_qr_or(x$qr_compact_factor_source, "none")))
  }
  cat("]\n")
  print(coef(x), ...)
  invisible(x)
}

coef.ridge_fit <- function(object, ...) {
  object$coefficients
}

fitted.ridge_fit <- function(object, ...) {
  object$fitted.values
}

residuals.ridge_fit <- function(object, ...) {
  object$residuals
}

print.ridge_fit <- function(x, ...) {
  cat(sprintf(
    "am_ridge_fit [backend=%s|precision=%s|lambda=%g|rank=%d|df.residual=%d|cache=%s]\n",
    x$backend,
    x$precision,
    x$lambda,
    x$rank,
    x$df.residual,
    if (isTRUE(x$cache_reused)) "hit" else "miss"
  ))
  print(coef(x), ...)
  invisible(x)
}
