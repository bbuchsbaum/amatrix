setClass(
  "amSVD",
  slots = list(
    u = "matrix",
    d = "numeric",
    v = "matrix",
    k = "integer",
    source_id = "character",
    precision = "character",
    backend = "character",
    ut_am = "ANY",  # adgeMatrix of t(u) for GPU matmul routing (matmul_min_dim=128); NULL on CPU
    v_am = "ANY"    # adgeMatrix of v for GPU matmul routing; NULL on CPU
  ),
  validity = function(object) {
    if (!identical(ncol(object@u), length(object@d))) {
      return("ncol(u) must equal length(d)")
    }
    if (!identical(ncol(object@v), length(object@d))) {
      return("ncol(v) must equal length(d)")
    }
    if (!identical(length(object@d), as.integer(object@k))) {
      return("length(d) must equal k")
    }
    TRUE
  }
)

setMethod("show", "amSVD", function(object) {
  cat(sprintf(
    "amSVD [%dx%d -> rank %d | %s | source: %s]\n  d[1:min(3,k)]: %s\n",
    nrow(object@u), nrow(object@v), object@k, object@precision, object@source_id,
    paste(round(object@d[seq_len(min(3L, object@k))], 4), collapse = ", ")
  ))
  invisible(object)
})

.amatrix_svd_factor_methods <- c("auto", "exact", "rsvd")

.amatrix_svd_factor_default_oversamples <- function(k) {
  as.integer(max(5L, min(10L, as.integer(k))))
}

.amatrix_svd_factor_rsvd_min_dim <- function() {
  as.integer(getOption("amatrix.svd_factor.rsvd_min_dim", 400L))
}

.amatrix_svd_factor_rsvd_max_rank_ratio <- function() {
  as.numeric(getOption("amatrix.svd_factor.rsvd_max_rank_ratio", 0.2))
}

.amatrix_backend_supports_capability <- function(name, capability, precision = NULL) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    return(FALSE)
  }

  if (!(name %in% amatrix_backend_names())) {
    return(FALSE)
  }

  backend <- tryCatch(.amatrix_get_backend(name), error = function(e) NULL)
  if (is.null(backend) || !isTRUE(backend$available())) {
    return(FALSE)
  }

  if (!is.null(precision) && !(precision %in% backend$precision_modes())) {
    return(FALSE)
  }

  capability %in% backend$capabilities()
}

.amatrix_svd_factor_rsvd_backend <- function(X) {
  if (!inherits(X, "adgeMatrix") || !identical(X@precision, "fast")) {
    return(NULL)
  }

  requested <- unique(c(X@preferred_backend, X@policy, amatrix_default_policy()))
  requested <- requested[nzchar(requested) & !(requested %in% c("auto", "cpu"))]

  for (backend_name in requested) {
    if (.amatrix_backend_supports_capability(backend_name, "rsvd", precision = X@precision)) {
      return(backend_name)
    }
  }

  # When the requested GPU backend lacks a truncated-SVD kernel, prefer MLX on
  # Apple Silicon if it is available rather than falling straight back to CPU.
  if (length(requested) > 0L && .amatrix_backend_supports_capability("mlx", "rsvd", precision = X@precision)) {
    return("mlx")
  }

  NULL
}

.amatrix_svd_factor_should_use_rsvd <- function(X, k, backend_name) {
  if (is.null(backend_name)) {
    return(FALSE)
  }

  dims <- dim(X)
  if (is.null(dims) || length(dims) != 2L) {
    return(FALSE)
  }

  min_dim <- min(dims)
  if (min_dim < .amatrix_svd_factor_rsvd_min_dim() || k >= min_dim) {
    return(FALSE)
  }

  (k / min_dim) <= .amatrix_svd_factor_rsvd_max_rank_ratio()
}

.amatrix_svd_factor_plan <- function(X, k, method, n_oversamples, n_iter) {
  method <- match.arg(method, .amatrix_svd_factor_methods)
  rsvd_backend <- .amatrix_svd_factor_rsvd_backend(X)

  if (identical(method, "auto")) {
    method <- if (.amatrix_svd_factor_should_use_rsvd(X, k, rsvd_backend)) "rsvd" else "exact"
  }

  factor_backend <- X@preferred_backend
  if (identical(method, "rsvd") && !is.null(rsvd_backend)) {
    factor_backend <- rsvd_backend
  }

  list(
    method = method,
    factor_backend = factor_backend,
    rsvd_backend = if (identical(method, "rsvd")) rsvd_backend else NULL,
    n_oversamples = as.integer(n_oversamples),
    n_iter = as.integer(n_iter)
  )
}

.amatrix_svd_cache_key <- function(X, k, plan) {
  paste(
    "svd",
    X@object_id,
    as.integer(k),
    plan$method,
    plan$factor_backend,
    if (identical(plan$method, "rsvd")) plan$n_oversamples else "exact",
    if (identical(plan$method, "rsvd")) plan$n_iter else "exact",
    sep = ":"
  )
}

am_svd_factor <- function(X,
                          k = min(dim(X)),
                          method = c("auto", "exact", "rsvd"),
                          n_oversamples = .amatrix_svd_factor_default_oversamples(k),
                          n_iter = 2L) {
  if (!inherits(X, "aMatrix")) {
    stop("X must be an aMatrix", call. = FALSE)
  }
  k <- as.integer(k)
  if (is.na(k) || k < 1L) {
    stop("k must be a positive integer", call. = FALSE)
  }
  max_k <- as.integer(min(dim(X)))
  if (k > max_k) {
    stop(sprintf("k (%d) cannot exceed min(dim(X)) (%d)", k, max_k), call. = FALSE)
  }

  n_oversamples <- as.integer(n_oversamples)
  n_iter <- as.integer(n_iter)
  if (is.na(n_oversamples) || n_oversamples < 0L) {
    stop("n_oversamples must be a non-negative integer", call. = FALSE)
  }
  if (is.na(n_iter) || n_iter < 0L) {
    stop("n_iter must be a non-negative integer", call. = FALSE)
  }

  plan <- .amatrix_svd_factor_plan(
    X = X,
    k = k,
    method = method,
    n_oversamples = n_oversamples,
    n_iter = n_iter
  )

  cache_key <- .amatrix_svd_cache_key(X, k, plan)
  cached <- get0(cache_key, envir = .amatrix_state$model_cache, inherits = FALSE)
  if (!is.null(cached)) {
    return(cached)
  }

  svd_result <- if (identical(plan$method, "rsvd")) {
    work_x <- X
    target_backend <- plan$rsvd_backend

    if (!inherits(work_x, "adgeMatrix")) {
      work_x <- adgeMatrix(
        amatrix_materialize_host(X),
        preferred_backend = if (is.null(target_backend)) X@preferred_backend else target_backend,
        policy = X@policy,
        precision = X@precision
      )
    } else if (!is.null(target_backend) && !identical(work_x@preferred_backend, target_backend)) {
      work_x <- adgeMatrix(
        amatrix_materialize_host(X),
        preferred_backend = target_backend,
        policy = X@policy,
        precision = X@precision
      )
    }

    am_rsvd(
      work_x,
      k = k,
      n_oversamples = plan$n_oversamples,
      n_iter = plan$n_iter
    )
  } else {
    am_svd(X, nu = k, nv = k)
  }

  d <- as.numeric(svd_result$d)
  if (length(d) > k) {
    d <- d[seq_len(k)]
  }

  u <- as.matrix(svd_result$u)
  v <- as.matrix(svd_result$v)

  backend <- plan$factor_backend
  use_gpu <- nzchar(backend) && backend != "cpu" && X@precision == "fast"
  # Pre-transpose u so am_matmul(ut_am, Y) routes through the matmul path
  # (min_dim threshold 128) instead of crossprod (threshold 2048).
  ut_am <- if (use_gpu) adgeMatrix(t(u), preferred_backend = backend, precision = X@precision) else NULL
  v_am  <- if (use_gpu) adgeMatrix(v,    preferred_backend = backend, precision = X@precision) else NULL

  factor <- methods::new(
    "amSVD",
    u = u,
    d = d,
    v = v,
    k = k,
    source_id = X@object_id,
    precision = X@precision,
    backend = backend,
    ut_am = ut_am,
    v_am = v_am
  )

  assign(cache_key, factor, envir = .amatrix_state$model_cache)
  factor
}

am_svd_project <- function(factor, Y) {
  if (!methods::is(factor, "amSVD")) {
    stop("factor must be an amSVD object", call. = FALSE)
  }
  Y_mat <- if (is.matrix(Y)) Y else matrix(as.numeric(Y), nrow = NROW(Y))
  ut_am <- factor@ut_am
  if (!is.null(ut_am)) {
    # GPU path: ut_am = t(u) as adgeMatrix; am_matmul routes through the
    # matmul path (min_dim threshold 128, vs 2048 for crossprod). t(u) is
    # stored once at factor creation time — no per-call transpose.
    if (nrow(Y_mat) != ncol(ut_am)) {
      stop(sprintf(
        "Y has %d rows but factor@u has %d rows",
        nrow(Y_mat), ncol(ut_am)
      ), call. = FALSE)
    }
    as.matrix(am_matmul(ut_am, Y_mat))
  } else {
    u <- factor@u
    if (nrow(Y_mat) != nrow(u)) {
      stop(sprintf(
        "Y has %d rows but factor@u has %d rows",
        nrow(Y_mat), nrow(u)
      ), call. = FALSE)
    }
    base::crossprod(u, Y_mat)
  }
}

am_svd_reconstruct <- function(factor, Z) {
  if (!methods::is(factor, "amSVD")) {
    stop("factor must be an amSVD object", call. = FALSE)
  }
  d <- factor@d
  Z_mat <- if (is.matrix(Z)) Z else matrix(as.numeric(Z), nrow = NROW(Z))
  if (nrow(Z_mat) != factor@k) {
    stop(sprintf(
      "Z has %d rows but factor@k is %d",
      nrow(Z_mat), factor@k
    ), call. = FALSE)
  }
  v_am <- factor@v_am
  if (!is.null(v_am)) {
    # GPU path: V %*% (Z / d) — scale Z columns by 1/d then matmul through backend
    as.matrix(am_matmul(v_am, Z_mat / d))
  } else {
    factor@v %*% (Z_mat / d)
  }
}

am_pca_coef <- function(factor, Y) {
  am_svd_reconstruct(factor, am_svd_project(factor, Y))
}
