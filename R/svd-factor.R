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

.amatrix_svd_cache_key <- function(X, k) {
  paste0("svd:", X@object_id, ":", as.integer(k))
}

am_svd_factor <- function(X, k = min(dim(X))) {
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

  cache_key <- .amatrix_svd_cache_key(X, k)
  cached <- get0(cache_key, envir = .amatrix_state$model_cache, inherits = FALSE)
  if (!is.null(cached)) {
    return(cached)
  }

  svd_result <- am_svd(X, nu = k, nv = k)

  d <- as.numeric(svd_result$d)
  if (length(d) > k) {
    d <- d[seq_len(k)]
  }

  u <- as.matrix(svd_result$u)
  v <- as.matrix(svd_result$v)

  backend <- X@preferred_backend
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
