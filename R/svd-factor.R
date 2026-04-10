setClass(
  "amSVD",
  slots = list(
    u = "matrix",
    d = "numeric",
    v = "matrix",
    k = "integer",
    method = "character",
    engine = "character",
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
    "amSVD [%dx%d -> rank %d | %s | %s/%s@%s | source: %s]\n  d[1:min(3,k)]: %s\n",
    nrow(object@u), nrow(object@v), object@k, object@precision,
    object@method, object@engine, object@backend, object@source_id,
    paste(round(object@d[seq_len(min(3L, object@k))], 4), collapse = ", ")
  ))
  invisible(object)
})

.amatrix_svd_factor_methods <- c("auto", "exact", "rsvd", "subspace")

.amatrix_svd_factor_default_oversamples <- function(k) {
  as.integer(max(5L, min(10L, as.integer(k))))
}

.amatrix_svd_factor_rsvd_min_dim <- function() {
  as.integer(getOption("amatrix.svd_factor.rsvd_min_dim", 400L))
}

.amatrix_svd_factor_rsvd_max_rank_ratio <- function() {
  as.numeric(getOption("amatrix.svd_factor.rsvd_max_rank_ratio", 0.2))
}

.amatrix_svd_factor_subspace_min_dim <- function() {
  as.integer(getOption("amatrix.svd_factor.subspace_min_dim", 400L))
}

.amatrix_svd_factor_subspace_min_rank_ratio <- function() {
  as.numeric(getOption("amatrix.svd_factor.subspace_min_rank_ratio", 0.1))
}

.amatrix_svd_factor_subspace_max_rank_ratio <- function() {
  as.numeric(getOption("amatrix.svd_factor.subspace_max_rank_ratio", 0.65))
}

.amatrix_backend_supports_capability <- function(name, capability, precision = NULL) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
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

.amatrix_svd_factor_subspace_backend <- function(X) {
  if (!inherits(X, "aMatrix")) {
    return(NULL)
  }

  requested <- unique(c(X@preferred_backend, X@policy, amatrix_default_policy()))
  requested <- requested[nzchar(requested) & !(requested %in% c("auto", "cpu"))]

  if (inherits(X, "adgeMatrix")) {
    for (backend_name in requested) {
      if (.amatrix_backend_supports_capability(backend_name, "rsvd", precision = X@precision)) {
        return(backend_name)
      }
    }
  }

  candidates <- .amatrix_resident_backend_candidates(X, op = "matmul")
  for (backend_name in candidates) {
    backend <- tryCatch(.amatrix_get_backend(backend_name), error = function(e) NULL)
    if (is.null(backend) ||
        !isTRUE(backend$available()) ||
        !.amatrix_backend_residency_capable(backend) ||
        !(X@precision %in% unique(backend$precision_modes()))) {
      next
    }

    if (.amatrix_backend_supports_resident_op(backend, "matmul", x = X) &&
        .amatrix_backend_supports_resident_op(backend, "crossprod", x = X)) {
      return(backend_name)
    }
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

.amatrix_svd_factor_should_use_subspace <- function(X, k, backend_name) {
  if (is.null(backend_name) || !inherits(X, c("adgeMatrix", "adgCMatrix"))) {
    return(FALSE)
  }

  if (!identical(X@precision, "fast")) {
    return(FALSE)
  }

  dims <- dim(X)
  if (is.null(dims) || length(dims) != 2L) {
    return(FALSE)
  }

  min_dim <- min(dims)
  if (min_dim < .amatrix_svd_factor_subspace_min_dim() || k >= min_dim) {
    return(FALSE)
  }

  rank_ratio <- k / min_dim
  rank_ratio >= .amatrix_svd_factor_subspace_min_rank_ratio() &&
    rank_ratio <= .amatrix_svd_factor_subspace_max_rank_ratio()
}

.amatrix_svd_factor_prepare_work_x <- function(X, target_backend = NULL) {
  work_x <- X

  if (inherits(work_x, "adgCMatrix")) {
    if (!is.null(target_backend) && !identical(work_x@preferred_backend, target_backend)) {
      work_x <- as_adgCMatrix(
        amatrix_materialize_host(X),
        preferred_backend = target_backend,
        policy = X@policy,
        precision = X@precision
      )
    }
  } else if (!inherits(work_x, "adgeMatrix")) {
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

  work_x
}

.amatrix_subspace_compile_operator <- function(work_x, op, target_backend = NULL) {
  backend_name <- if (is.null(target_backend) || !nzchar(target_backend)) "auto" else target_backend
  tryCatch(
    amatrix_compile_product(
      work_x,
      op = op,
      backend = backend_name,
      precision = work_x@precision,
      policy = work_x@policy
    ),
    error = function(e) NULL
  )
}

.amatrix_subspace_rademacher <- function(n, r) {
  draws <- ifelse(stats::runif(n * r) < 0.5, -1.0, 1.0)
  matrix(draws, nrow = n, ncol = r)
}

.amatrix_subspace_qr <- function(x) {
  qr_obj <- qr(x)
  list(
    q = as.matrix(qr.Q(qr_obj, complete = FALSE)),
    r = as.matrix(qr.R(qr_obj, complete = FALSE))
  )
}

.amatrix_subspace_trim_q <- function(q, r, k, eps_rank) {
  diag_r <- if (length(r) > 0L) abs(base::diag(as.matrix(r))) else numeric(0L)
  keep <- which(diag_r > eps_rank)

  if (length(keep) == 0L) {
    keep <- seq_len(min(k, ncol(q)))
  }

  list(
    q = q[, keep, drop = FALSE],
    diag_history = if (length(diag_r) > 0L) diag_r[keep] else numeric(0L)
  )
}

.amatrix_subspace_core_solver <- function(k_eff, rank_discovered, min_dim) {
  if (rank_discovered < k_eff && rank_discovered <= 64L) {
    return("svd_core")
  }

  if (k_eff >= 128L || (k_eff / min_dim) >= 0.12) {
    return("gram")
  }

  "qr"
}

.amatrix_subspace_svd <- function(X,
                                  k,
                                  n_oversamples = .amatrix_svd_factor_default_oversamples(k),
                                  n_iter = 2L,
                                  target_backend = NULL,
                                  eps_rank = 1e-8) {
  work_x <- .amatrix_svd_factor_prepare_work_x(X, target_backend = target_backend)
  left_plan <- .amatrix_subspace_compile_operator(work_x, op = "matmul", target_backend = target_backend)
  right_plan <- .amatrix_subspace_compile_operator(work_x, op = "crossprod", target_backend = target_backend)
  on.exit(.amatrix_release_product_plan(left_plan), add = TRUE)
  on.exit(.amatrix_release_product_plan(right_plan), add = TRUE)

  left_apply <- function(rhs) {
    if (inherits(left_plan, "am_product_plan")) {
      return(as.matrix(left_plan(rhs)))
    }
    as.matrix(work_x %*% rhs)
  }

  right_apply <- function(rhs) {
    if (inherits(right_plan, "am_product_plan")) {
      return(as.matrix(right_plan(rhs)))
    }
    as.matrix(crossprod(work_x, rhs))
  }

  dims <- dim(work_x)
  n <- dims[[2L]]
  min_dim <- min(dims)

  target_rank <- min(as.integer(k + n_oversamples), min_dim)

  if (inherits(work_x, "adgeMatrix") && !is.null(target_backend)) {
    backend <- tryCatch(.amatrix_get_backend(target_backend), error = function(e) NULL)
    if (!is.null(backend) && isTRUE(backend$available()) && is.function(backend$rsvd)) {
      fast_res <- backend$rsvd(
        work_x,
        k = as.integer(k),
        n_oversamples = as.integer(n_oversamples),
        n_iter = as.integer(n_iter)
      )
      return(list(
        u = fast_res$u,
        d = as.numeric(fast_res$d),
        v = fast_res$v,
        rank_discovered = length(fast_res$d),
        core_solver = "backend_rsvd",
        diag_history = as.numeric(fast_res$d)
      ))
    }
  }

  Omega <- .amatrix_subspace_rademacher(n = n, r = target_rank)

  phase <- .amatrix_subspace_qr(left_apply(Omega))
  q_basis <- phase$q
  r_basis <- phase$r

  if (n_iter > 0L) {
    for (iter in seq_len(n_iter)) {
      right_phase <- .amatrix_subspace_qr(right_apply(q_basis))
      phase <- .amatrix_subspace_qr(left_apply(right_phase$q))
      q_basis <- phase$q
      r_basis <- phase$r
    }
  }

  trimmed <- .amatrix_subspace_trim_q(
    q = q_basis,
    r = r_basis,
    k = k,
    eps_rank = eps_rank
  )
  q_basis <- trimmed$q
  rank_discovered <- ncol(q_basis)
  if (rank_discovered < 1L) {
    stop("subspace SVD did not discover a usable range space", call. = FALSE)
  }

  t_mat <- right_apply(q_basis)
  k_eff <- min(as.integer(k), rank_discovered)
  core_solver <- .amatrix_subspace_core_solver(
    k_eff = k_eff,
    rank_discovered = rank_discovered,
    min_dim = min_dim
  )

  if (identical(core_solver, "gram")) {
    gram <- base::crossprod(t_mat)
    eig <- eigen((gram + base::t(gram)) / 2, symmetric = TRUE)
    ord <- order(eig$values, decreasing = TRUE)
    lam <- pmax(eig$values[ord], 0)
    u_core <- eig$vectors[, ord, drop = FALSE][, seq_len(k_eff), drop = FALSE]
    d <- sqrt(lam[seq_len(k_eff)])
    d[d < 1e-12] <- 1e-12
    u <- q_basis %*% u_core
    v <- t_mat %*% u_core
    v <- sweep(v, 2L, d, "/", check.margin = FALSE)
  } else if (identical(core_solver, "svd_core")) {
    b_core <- t(t_mat)
    sv <- base::svd(b_core, nu = k_eff, nv = k_eff)
    u <- q_basis %*% sv$u[, seq_len(k_eff), drop = FALSE]
    d <- sv$d[seq_len(k_eff)]
    v <- sv$v[, seq_len(k_eff), drop = FALSE]
  } else {
    qr_t <- qr(t_mat)
    rank_t <- qr_t$rank
    if (rank_t < 1L) {
      stop("projected core is numerically rank-deficient", call. = FALSE)
    }

    pivot <- qr_t$pivot
    if (is.null(pivot)) {
      pivot <- seq_len(ncol(t_mat))
    }

    q_b <- qr.Q(qr_t, complete = FALSE)[, seq_len(rank_t), drop = FALSE]
    r_full <- qr.R(qr_t, complete = FALSE)
    r_reduced <- r_full[seq_len(rank_t), , drop = FALSE]
    tiny <- t(r_reduced)
    sv <- base::svd(tiny, nu = k_eff, nv = k_eff)

    u <- q_basis[, pivot, drop = FALSE] %*% sv$u[, seq_len(k_eff), drop = FALSE]
    d <- sv$d[seq_len(k_eff)]
    v <- q_b %*% sv$v[, seq_len(k_eff), drop = FALSE]
  }

  list(
    u = u,
    d = as.numeric(d),
    v = v,
    rank_discovered = rank_discovered,
    core_solver = core_solver,
    diag_history = trimmed$diag_history
  )
}

.amatrix_svd_factor_plan <- function(X, k, method, n_oversamples, n_iter) {
  method <- match.arg(method, .amatrix_svd_factor_methods)
  rsvd_backend <- .amatrix_svd_factor_rsvd_backend(X)
  subspace_backend <- .amatrix_svd_factor_subspace_backend(X)

  if (identical(method, "auto")) {
    if (.amatrix_svd_factor_should_use_rsvd(X, k, rsvd_backend)) {
      method <- "rsvd"
    } else if (.amatrix_svd_factor_should_use_subspace(X, k, subspace_backend)) {
      method <- "subspace"
    } else {
      method <- "exact"
    }
  }

  factor_backend <- X@preferred_backend
  if (identical(method, "rsvd") && !is.null(rsvd_backend)) {
    factor_backend <- rsvd_backend
  } else if (identical(method, "subspace") && !is.null(subspace_backend)) {
    factor_backend <- subspace_backend
  }

  list(
    method = method,
    factor_backend = factor_backend,
    rsvd_backend = if (identical(method, "rsvd")) rsvd_backend else NULL,
    subspace_backend = if (identical(method, "subspace")) subspace_backend else NULL,
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
    if (identical(plan$method, "exact")) "exact" else plan$n_oversamples,
    if (identical(plan$method, "exact")) "exact" else plan$n_iter,
    sep = ":"
  )
}

.amatrix_svd_factor_engine <- function(plan, svd_result) {
  if (identical(plan$method, "exact")) {
    return("exact_svd")
  }

  if (identical(plan$method, "rsvd")) {
    if (!is.null(plan$rsvd_backend)) {
      return("backend_rsvd")
    }
    if (requireNamespace("irlba", quietly = TRUE)) {
      return("irlba_svdr")
    }
    return("base_svd")
  }

  engine <- svd_result$core_solver
  if (is.null(engine) || !nzchar(as.character(engine[[1L]]))) {
    return("subspace_fallback")
  }

  as.character(engine[[1L]])
}

.amatrix_svd_factor_backend <- function(X, plan, svd_result) {
  if (identical(plan$method, "exact")) {
    return(.amatrix_backend_for(X, "svd")$name)
  }

  if (identical(plan$method, "rsvd")) {
    if (!is.null(plan$rsvd_backend)) {
      return(plan$rsvd_backend)
    }
    return("cpu")
  }

  if (!is.null(plan$subspace_backend) && nzchar(plan$subspace_backend)) {
    return(plan$subspace_backend)
  }

  if (!is.null(svd_result$core_solver) && identical(as.character(svd_result$core_solver[[1L]]), "backend_rsvd")) {
    return("cpu")
  }

  "cpu"
}

svd_factor <- function(X,
                          k = min(dim(X)),
                          method = c("auto", "exact", "rsvd", "subspace"),
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
  cached <- .amatrix_cache_get(cache_key)
  if (!is.null(cached)) {
    return(cached)
  }

  svd_result <- if (identical(plan$method, "rsvd")) {
    work_x <- .amatrix_svd_factor_prepare_work_x(
      X,
      target_backend = plan$rsvd_backend
    )
    rsvd(
      work_x,
      k = k,
      n_oversamples = plan$n_oversamples,
      n_iter = plan$n_iter
    )
  } else if (identical(plan$method, "subspace")) {
    .amatrix_subspace_svd(
      X,
      k = k,
      n_oversamples = plan$n_oversamples,
      n_iter = plan$n_iter,
      target_backend = plan$subspace_backend
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

  engine <- .amatrix_svd_factor_engine(plan, svd_result)
  backend <- .amatrix_svd_factor_backend(X, plan, svd_result)
  use_gpu <- nzchar(backend) && backend != "cpu" && X@precision == "fast"
  # Pre-transpose u so matmul(ut_am, Y) routes through the matmul path
  # (min_dim threshold 128) instead of crossprod (threshold 2048).
  ut_am <- if (use_gpu) adgeMatrix(t(u), preferred_backend = backend, precision = X@precision) else NULL
  v_am  <- if (use_gpu) adgeMatrix(v,    preferred_backend = backend, precision = X@precision) else NULL

  factor <- methods::new(
    "amSVD",
    u = u,
    d = d,
    v = v,
    k = k,
    method = plan$method,
    engine = engine,
    source_id = X@object_id,
    precision = X@precision,
    backend = backend,
    ut_am = ut_am,
    v_am = v_am
  )

  .amatrix_cache_set(cache_key, factor)
  factor
}

svd_project <- function(factor, Y) {
  if (!methods::is(factor, "amSVD")) {
    stop("factor must be an amSVD object", call. = FALSE)
  }
  Y_mat <- if (is.matrix(Y)) Y else matrix(as.numeric(Y), nrow = NROW(Y))
  ut_am <- factor@ut_am
  if (!is.null(ut_am)) {
    # GPU path: ut_am = t(u) as adgeMatrix; matmul routes through the
    # matmul path (min_dim threshold 128, vs 2048 for crossprod). t(u) is
    # stored once at factor creation time — no per-call transpose.
    if (nrow(Y_mat) != ncol(ut_am)) {
      stop(sprintf(
        "Y has %d rows but factor@u has %d rows",
        nrow(Y_mat), ncol(ut_am)
      ), call. = FALSE)
    }
    as.matrix(matmul(ut_am, Y_mat))
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

svd_reconstruct <- function(factor, Z) {
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
    as.matrix(matmul(v_am, Z_mat / d))
  } else {
    factor@v %*% (Z_mat / d)
  }
}

pca_coef <- function(factor, Y) {
  svd_reconstruct(factor, svd_project(factor, Y))
}
