amatrix_mlx_capabilities <- function() {
  c("matmul", "crossprod", "tcrossprod", "ewise", "broadcast_ewise", "argmax", "scatter_mean", "segment_sum", "segment_mean", "addmm",
    "rowSums", "colSums",
    "qr", "svd", "rsvd", "chol", "chol_gpu", "batched_trsm", "eigen", "covariance")
}

amatrix_mlx_features <- function() {
  c("dense_f32", "resident_dense", "unified_memory", "custom_ops",
    "qr", "rsvd", "chol_gpu", "batched_trsm", "eigen_sym", "sparse_spmm")
}

amatrix_mlx_precision_modes <- function() {
  "fast"
}

amatrix_mlx_native_available <- function() {
  .Call("amatrix_mlx_native_available_bridge")
}

# Activate the Metal GPU probe for the current session.  Safe to call from
# Rscript -e / interactive / testthat contexts.  Do NOT call from the body
# of a plain `Rscript file.R` script — Metal device init crashes in that
# launch mode (upstream MLX bug: https://github.com/ml-explore/mlx/issues/2691).
amatrix_mlx_enable_gpu_probe <- function() {
  Sys.setenv(AMATRIX_MLX_PROBE_GPU = "1")
  invisible(amatrix_mlx_native_available())
}

amatrix_mlx_is_available <- function() {
  isTRUE(getOption("amatrix.mlx.available", FALSE)) || isTRUE(amatrix_mlx_native_available())
}

amatrix_mlx_bridge_info <- function() {
  info <- .Call("amatrix_mlx_bridge_info_bridge")
  info$available <- amatrix_mlx_is_available()
  info$capabilities <- amatrix_mlx_capabilities()
  info
}

amatrix_mlx_matmul <- function(x, y) {
  x_mat <- as.matrix(x)
  y_mat <- as.matrix(y)

  if (!is.double(x_mat)) {
    storage.mode(x_mat) <- "double"
  }

  if (!is.double(y_mat)) {
    storage.mode(y_mat) <- "double"
  }

  .Call("amatrix_mlx_matmul_bridge", x_mat, y_mat)
}

amatrix_mlx_crossprod <- function(x, y = NULL) {
  x_mat <- as.matrix(x)
  y_mat <- if (is.null(y)) NULL else as.matrix(y)

  if (!is.double(x_mat)) {
    storage.mode(x_mat) <- "double"
  }

  if (!is.null(y_mat) && !is.double(y_mat)) {
    storage.mode(y_mat) <- "double"
  }

  .Call("amatrix_mlx_crossprod_bridge", x_mat, y_mat)
}

amatrix_mlx_tcrossprod <- function(x, y = NULL) {
  x_mat <- as.matrix(x)
  y_mat <- if (is.null(y)) NULL else as.matrix(y)

  if (!is.double(x_mat)) {
    storage.mode(x_mat) <- "double"
  }

  if (!is.null(y_mat) && !is.double(y_mat)) {
    storage.mode(y_mat) <- "double"
  }

  .Call("amatrix_mlx_tcrossprod_bridge", x_mat, y_mat)
}

amatrix_mlx_spmm <- function(x_sp, y, trans_lhs = FALSE) {
  # x_sp: dgCMatrix (materialized from adgCMatrix via amatrix_materialize_host)
  # y:    dense host matrix
  # trans_lhs=TRUE: compute t(x_sp) %*% y; FALSE: x_sp %*% y
  y_mat <- if (is.matrix(y)) y else as.matrix(y)
  if (!is.double(y_mat)) storage.mode(y_mat) <- "double"
  .Call("amatrix_mlx_spmm_bridge",
        as.double(x_sp@x), as.integer(x_sp@p), as.integer(x_sp@i),
        as.integer(x_sp@Dim), y_mat, as.logical(trans_lhs),
        PACKAGE = "amatrix.mlx")
}

amatrix_mlx_ewise <- function(lhs, rhs = NULL, op) {
  lhs_mat <- as.matrix(lhs)
  rhs_mat <- if (is.null(rhs)) NULL else as.matrix(rhs)

  if (!is.double(lhs_mat)) {
    storage.mode(lhs_mat) <- "double"
  }

  if (!is.null(rhs_mat) && !is.double(rhs_mat)) {
    storage.mode(rhs_mat) <- "double"
  }

  .Call("amatrix_mlx_ewise_bridge", lhs_mat, rhs_mat, op)
}

amatrix_mlx_axis_sums <- function(x, axis) {
  x_mat <- as.matrix(x)
  if (!is.double(x_mat)) {
    storage.mode(x_mat) <- "double"
  }
  .Call("amatrix_mlx_sum_axis_bridge", x_mat, as.integer(axis))
}

.amatrix_mlx_qr_rank_from_r <- function(r_mat) {
  stopifnot(is.matrix(r_mat))
  diag_len <- min(dim(r_mat))
  if (diag_len == 0L) {
    return(0L)
  }
  diag_abs <- abs(diag(r_mat[seq_len(diag_len), seq_len(diag_len), drop = FALSE]))
  scale <- if (length(diag_abs) == 0L) 0 else max(diag_abs)
  tol <- max(dim(r_mat)) * .Machine$double.eps * scale
  as.integer(sum(diag_abs > tol))
}

.amatrix_mlx_qr_explicit_raw <- function(x_mat) {
  q_key <- amatrix:::.amatrix_next_resident_key("mlx")
  qr_raw <- .Call("amatrix_mlx_qr_bridge", x_mat, q_key)
  qr_raw$representation <- "explicit_qr"
  qr_raw$rank <- .amatrix_mlx_qr_rank_from_r(qr_raw$r)
  qr_raw$backend_ops <- "mlx"
  qr_raw
}

.amatrix_mlx_host_qr_factor_builder <- function(x_mat) {
  force(x_mat)
  function() base::qr(x_mat)
}

.amatrix_mlx_qr_explicit <- function(x_mat, include_factor = FALSE, lazy_factor = FALSE) {
  qr_raw <- .amatrix_mlx_qr_explicit_raw(x_mat)
  c(
    qr_raw,
    list(
      factor = if (isTRUE(include_factor)) base::qr(x_mat) else NULL,
      factor_builder = if (!isTRUE(include_factor) && isTRUE(lazy_factor)) .amatrix_mlx_host_qr_factor_builder(x_mat) else NULL,
      factor_source = if (isTRUE(include_factor) || isTRUE(lazy_factor)) "host_compact" else "reconstructable"
    )
  )
}

.amatrix_mlx_qr_compact_method <- function() {
  mode <- getOption("amatrix.mlx.qr_compact_method", "auto")
  match.arg(mode, c("auto", "bridge", "tsqr"))
}

.amatrix_mlx_qr_block_rows <- function(n, p) {
  override <- getOption("amatrix.mlx.qr_tsqr_block_rows", NULL)
  if (!is.null(override)) {
    return(max(1L, min(as.integer(override), as.integer(n))))
  }
  as.integer(max(2L * p, ceiling(n / 4L)))
}

.amatrix_mlx_qr_use_tsqr <- function(x_mat) {
  mode <- .amatrix_mlx_qr_compact_method()
  if (identical(mode, "bridge")) {
    return(FALSE)
  }
  if (identical(mode, "tsqr")) {
    return(TRUE)
  }
  n <- nrow(x_mat)
  p <- ncol(x_mat)
  isTRUE(n >= (4L * p) && n > p)
}

amatrix_mlx_qr_cache_signature <- function(x) {
  dims <- if (length(x) == 2L && is.numeric(x) && is.null(dim(x))) {
    as.integer(x)
  } else {
    as.integer(dim(x))
  }

  if (length(dims) != 2L || anyNA(dims)) {
    stop("x must supply two dimensions", call. = FALSE)
  }

  helper_mode <- .amatrix_mlx_qr_helper_mode()
  if (identical(helper_mode, "native")) {
    return("mlx:native")
  }

  compact_method <- .amatrix_mlx_qr_compact_method()
  n <- dims[[1]]
  p <- dims[[2]]
  use_tsqr <- identical(compact_method, "tsqr") || (identical(compact_method, "auto") && isTRUE(n >= (4L * p) && n > p))
  if (!use_tsqr) {
    return("mlx:compact:bridge")
  }

  sprintf("mlx:compact:tsqr:%d", .amatrix_mlx_qr_block_rows(n, p))
}

.amatrix_mlx_tsqr_block_factor <- function(x_mat, block_rows = NULL) {
  stopifnot(is.matrix(x_mat), is.double(x_mat))
  n <- nrow(x_mat)
  p <- ncol(x_mat)
  if (is.null(block_rows)) {
    block_rows <- .amatrix_mlx_qr_block_rows(n, p)
  }
  nblocks <- ceiling(n / block_rows)
  block_q_keys <- vapply(seq_len(nblocks), function(...) amatrix:::.amatrix_next_resident_key("mlx"), character(1))
  top_q_key <- amatrix:::.amatrix_next_resident_key("mlx")
  top_r_key <- amatrix:::.amatrix_next_resident_key("mlx")
  r_stack_key <- amatrix:::.amatrix_next_resident_key("mlx")
  payload <- .Call(
    "amatrix_mlx_tsqr_build_bridge",
    x_mat,
    as.integer(block_rows),
    block_q_keys,
    top_q_key,
    top_r_key,
    r_stack_key
  )
  factor_cache <- new.env(parent = emptyenv())
  factor_cache$block_factors <- vector("list", length(payload$block_rows))
  factor_cache$top_factor <- NULL
  factor_cache$r_stack <- NULL
  factor_cache$r <- NULL
  structure(
    list(
      kind = "tsqr_blocked",
      source_dim = c(n, p),
      block_rows = as.integer(payload$block_rows),
      cache_env = factor_cache,
      block_q_keys = as.character(payload$block_q_keys),
      r_stack_key = as.character(payload$r_stack_key),
      top_q_key = as.character(payload$top_q_key),
      top_r_key = as.character(payload$top_r_key),
      rank = as.integer(payload$rank),
      thin = TRUE,
      pivot = NULL,
      pivoted = FALSE
    ),
    class = "amatrix_mlx_tsqr_factor"
  )
}

.amatrix_mlx_tsqr_block_slices <- function(factor) {
  widths <- as.integer(factor$block_rows)
  ends <- cumsum(widths)
  starts <- ends - widths + 1L
  Map(seq.int, starts, ends)
}

.amatrix_mlx_tsqr_head_slices <- function(factor) {
  p <- factor$source_dim[[2]]
  count <- length(factor$block_rows)
  starts <- ((seq_len(count) - 1L) * p) + 1L
  ends <- starts + p - 1L
  Map(seq.int, starts, ends)
}

.amatrix_mlx_tsqr_get_block_factor <- function(factor, idx) {
  cached <- factor$cache_env$block_factors[[idx]]
  if (!is.null(cached)) {
    return(cached)
  }

  q_mat <- amatrix_mlx_resident_materialize(factor$block_q_keys[[idx]])
  r_rows <- .amatrix_mlx_tsqr_head_slices(factor)[[idx]]
  fac <- base::qr(q_mat %*% .amatrix_mlx_tsqr_r_stack(factor)[r_rows, , drop = FALSE])
  factor$cache_env$block_factors[[idx]] <- fac
  fac
}

.amatrix_mlx_tsqr_r_stack <- function(factor) {
  cached <- factor$cache_env$r_stack
  if (!is.null(cached)) {
    return(cached)
  }

  r_stack <- amatrix_mlx_resident_materialize(factor$r_stack_key)
  factor$cache_env$r_stack <- r_stack
  r_stack
}

.amatrix_mlx_tsqr_r <- function(factor) {
  cached <- factor$cache_env$r
  if (!is.null(cached)) {
    return(cached)
  }

  r <- amatrix_mlx_resident_materialize(factor$top_r_key)
  factor$cache_env$r <- r
  r
}

.amatrix_mlx_tsqr_get_top_factor <- function(factor) {
  cached <- factor$cache_env$top_factor
  if (!is.null(cached)) {
    return(cached)
  }

  q_top <- amatrix_mlx_resident_materialize(factor$top_q_key)
  fac <- base::qr(q_top %*% .amatrix_mlx_tsqr_r(factor))
  factor$cache_env$top_factor <- fac
  fac
}

.amatrix_mlx_tsqr_local_head_qty <- function(factor, y_mat) {
  row_slices <- .amatrix_mlx_tsqr_block_slices(factor)
  p <- factor$source_dim[[2]]
  Map(
    function(rows, q_key, idx) {
      y_block <- y_mat[rows, , drop = FALSE]
      if (!is.null(q_key) && nzchar(q_key)) {
        return(amatrix_mlx_qr_qty_key(q_key, y_block))
      }
      base::qr.qty(.amatrix_mlx_tsqr_get_block_factor(factor, idx), y_block)[seq_len(p), , drop = FALSE]
    },
    row_slices,
    as.list(factor$block_q_keys),
    seq_along(row_slices)
  )
}

.amatrix_mlx_tsqr_top_qty <- function(factor, head_stack) {
  base::qr.qty(.amatrix_mlx_tsqr_get_top_factor(factor), head_stack)
}

amatrix_mlx_tsqr_qty <- function(factor, y) {
  y_mat <- as.matrix(y)
  row_slices <- .amatrix_mlx_tsqr_block_slices(factor)
  head_slices <- .amatrix_mlx_tsqr_head_slices(factor)
  local_qty <- Map(
    function(rows, idx) {
      base::qr.qty(.amatrix_mlx_tsqr_get_block_factor(factor, idx), y_mat[rows, , drop = FALSE])
    },
    row_slices,
    seq_along(row_slices)
  )
  head_stack <- do.call(
    rbind,
    lapply(local_qty, function(qi) qi[seq_len(factor$source_dim[[2]]), , drop = FALSE])
  )
  top_qty <- .amatrix_mlx_tsqr_top_qty(factor, head_stack)
  for (idx in seq_along(local_qty)) {
    local_qty[[idx]][seq_len(factor$source_dim[[2]]), ] <- top_qty[head_slices[[idx]], , drop = FALSE]
  }
  do.call(rbind, local_qty)
}

amatrix_mlx_tsqr_qy <- function(factor, y) {
  y_mat <- as.matrix(y)
  row_slices <- .amatrix_mlx_tsqr_block_slices(factor)
  head_slices <- .amatrix_mlx_tsqr_head_slices(factor)
  local_in <- lapply(row_slices, function(rows) y_mat[rows, , drop = FALSE])
  head_stack <- do.call(
    rbind,
    lapply(local_in, function(qi) qi[seq_len(factor$source_dim[[2]]), , drop = FALSE])
  )
  top_qy <- base::qr.qy(.amatrix_mlx_tsqr_get_top_factor(factor), head_stack)
  out <- vector("list", length(local_in))
  for (idx in seq_along(local_in)) {
    local_in[[idx]][seq_len(factor$source_dim[[2]]), ] <- top_qy[head_slices[[idx]], , drop = FALSE]
    out[[idx]] <- base::qr.qy(.amatrix_mlx_tsqr_get_block_factor(factor, idx), local_in[[idx]])
  }
  do.call(rbind, out)
}

amatrix_mlx_tsqr_coef <- function(factor, y) {
  y_mat <- as.matrix(y)
  p <- factor$source_dim[[2]]
  rank <- factor$rank
  top_q_key <- factor$top_q_key
  top_r_key <- factor$top_r_key
  if (!is.null(top_q_key) && nzchar(top_q_key) && !is.null(top_r_key) && nzchar(top_r_key) && identical(rank, p)) {
    return(.Call(
      "amatrix_mlx_tsqr_coef_resident_bridge",
      as.character(factor$block_q_keys),
      as.integer(factor$block_rows),
      as.character(top_q_key),
      as.character(top_r_key),
      y_mat
    ))
  }
  head_stack <- do.call(rbind, .amatrix_mlx_tsqr_local_head_qty(factor, y_mat))
  top_qty <- .amatrix_mlx_tsqr_top_qty(factor, head_stack)
  coef <- matrix(NA_real_, nrow = p, ncol = ncol(y_mat))
  if (rank > 0L) {
    r_top <- .amatrix_mlx_tsqr_r(factor)[seq_len(rank), seq_len(rank), drop = FALSE]
    qty_top <- top_qty[seq_len(rank), , drop = FALSE]
    coef[seq_len(rank), ] <- amatrix_mlx_solve_triangular(r_top, qty_top, upper = TRUE)
  }
  coef
}

amatrix_mlx_tsqr_solve <- function(factor, b = NULL, tol = 1e-07) {
  source_dim <- factor$source_dim
  p <- source_dim[[2]]
  rank <- as.integer(factor$rank)

  if (is.null(b)) {
    if (source_dim[[1]] != source_dim[[2]]) {
      stop("only square matrices can be inverted", call. = FALSE)
    }
    b <- diag(p)
  } else {
    b <- as.matrix(b)
  }

  if (rank < p) {
    stop("singular matrix 'a' in solve", call. = FALSE)
  }

  amatrix_mlx_tsqr_coef(factor, b)
}

amatrix_mlx_tsqr_fitted <- function(factor, y, k = NULL) {
  y_mat <- as.matrix(y)
  qty <- amatrix_mlx_tsqr_qty(factor, y_mat)
  if (is.null(k)) {
    k <- factor$rank
  }
  if (k <= 0L) {
    return(matrix(0, nrow = factor$source_dim[[1]], ncol = ncol(y_mat)))
  }
  qty_trunc <- qty
  if (k < nrow(qty_trunc)) {
    qty_trunc[(k + 1L):nrow(qty_trunc), ] <- 0
  }
  amatrix_mlx_tsqr_qy(factor, qty_trunc)
}

amatrix_mlx_tsqr_resid <- function(factor, y) {
  y_mat <- as.matrix(y)
  y_mat - amatrix_mlx_tsqr_fitted(factor, y_mat)
}

amatrix_mlx_tsqr_q <- function(factor, complete = FALSE) {
  n <- factor$source_dim[[1]]
  p <- factor$source_dim[[2]]
  seed <- if (isTRUE(complete)) {
    diag(n)
  } else {
    rbind(diag(p), matrix(0, nrow = n - p, ncol = p))
  }
  amatrix_mlx_tsqr_qy(factor, seed)
}

amatrix_mlx_qr <- function(x) {
  x_mat <- as.matrix(x)
  if (!is.double(x_mat)) {
    storage.mode(x_mat) <- "double"
  }
  helper_mode <- .amatrix_mlx_qr_helper_mode()
  if (identical(helper_mode, "compact") && .amatrix_mlx_qr_use_tsqr(x_mat)) {
    factor <- .amatrix_mlx_tsqr_block_factor(x_mat)
    return(list(
      representation = "mlx_compact_qr",
      r_key = factor$top_r_key,
      rank = factor$rank,
      factor = factor,
      factor_source = "tsqr_blocked",
      backend_ops = "mlx"
    ))
  }

  qr_raw <- .amatrix_mlx_qr_explicit(
    x_mat,
    include_factor = FALSE,
    lazy_factor = identical(helper_mode, "compact")
  )
  qr_raw$representation <- if (identical(helper_mode, "compact")) "mlx_compact_qr" else "explicit_qr"
  qr_raw
}

.amatrix_mlx_qr_helper_mode <- function() {
  mode <- getOption("amatrix.mlx.qr_helper_mode", "native")
  match.arg(mode, c("native", "compact"))
}

amatrix_mlx_solve_triangular <- function(a, b, upper = TRUE) {
  a_mat <- as.matrix(a)
  b_mat <- as.matrix(b)

  if (!is.double(a_mat)) {
    storage.mode(a_mat) <- "double"
  }
  if (!is.double(b_mat)) {
    storage.mode(b_mat) <- "double"
  }

  .Call("amatrix_mlx_solve_triangular_bridge", a_mat, b_mat, as.logical(upper))
}

amatrix_mlx_qr_qty_key <- function(q_key, y) {
  y_mat <- as.matrix(y)
  if (!is.double(y_mat)) {
    storage.mode(y_mat) <- "double"
  }
  .Call("amatrix_mlx_qr_qty_key_bridge", as.character(q_key), y_mat)
}

amatrix_mlx_qr_qy_key <- function(q_key, y) {
  y_mat <- as.matrix(y)
  if (!is.double(y_mat)) {
    storage.mode(y_mat) <- "double"
  }
  .Call("amatrix_mlx_qr_qy_key_bridge", as.character(q_key), y_mat)
}

amatrix_mlx_qr_coef_key <- function(q_key, r, y) {
  r_mat <- as.matrix(r)
  y_mat <- as.matrix(y)

  if (!is.double(r_mat)) {
    storage.mode(r_mat) <- "double"
  }
  if (!is.double(y_mat)) {
    storage.mode(y_mat) <- "double"
  }

  .Call("amatrix_mlx_qr_coef_key_bridge", as.character(q_key), r_mat, y_mat)
}

amatrix_mlx_tsqr_coef_key <- function(q_keys, block_rows, top_q_key, r, y) {
  r_mat <- as.matrix(r)
  y_mat <- as.matrix(y)

  if (!is.double(r_mat)) {
    storage.mode(r_mat) <- "double"
  }
  if (!is.double(y_mat)) {
    storage.mode(y_mat) <- "double"
  }

  .Call(
    "amatrix_mlx_tsqr_coef_key_bridge",
    as.character(q_keys),
    as.integer(block_rows),
    as.character(top_q_key),
    r_mat,
    y_mat
  )
}

amatrix_mlx_qr_qty <- function(q, y) {
  if (is.character(q) && length(q) == 1L) {
    return(amatrix_mlx_qr_qty_key(q, y))
  }
  amatrix_mlx_crossprod(q, y)
}

amatrix_mlx_qr_qy <- function(q, y) {
  if (is.character(q) && length(q) == 1L) {
    return(amatrix_mlx_qr_qy_key(q, y))
  }
  amatrix_mlx_matmul(q, y)
}

amatrix_mlx_qr_coef <- function(q, r, y, rank = NULL) {
  r_mat <- as.matrix(r)
  y_mat <- as.matrix(y)
  p <- ncol(r_mat)
  rank <- if (is.null(rank)) p else as.integer(rank)
  if (is.character(q) && length(q) == 1L && identical(rank, p)) {
    return(amatrix_mlx_qr_coef_key(q, r_mat, y_mat))
  }
  q_mat <- as.matrix(q)
  qty <- amatrix_mlx_qr_qty(q_mat, y_mat)
  coef <- matrix(NA_real_, nrow = p, ncol = ncol(y_mat))

  if (rank > 0L) {
    r_top <- r_mat[seq_len(rank), seq_len(rank), drop = FALSE]
    qty_top <- qty[seq_len(rank), , drop = FALSE]
    coef[seq_len(rank), ] <- amatrix_mlx_solve_triangular(r_top, qty_top, upper = TRUE)
  }

  coef
}

amatrix_mlx_qr_fitted <- function(q, y, rank = NULL) {
  q_mat <- if (is.character(q) && length(q) == 1L) NULL else as.matrix(q)
  y_mat <- as.matrix(y)
  if (is.null(q_mat) && is.null(rank)) {
    stop("rank must be supplied when q is an MLX resident key")
  }
  p <- if (is.null(q_mat)) ncol(y_mat) else ncol(q_mat)
  rank <- if (is.null(rank)) p else as.integer(rank)
  qty <- amatrix_mlx_qr_qty(if (is.null(q_mat)) q else q_mat, y)

  if (rank <= 0L) {
    return(matrix(0, nrow = if (is.null(q_mat)) nrow(y_mat) else nrow(q_mat), ncol = ncol(y_mat)))
  }

  if (is.character(q) && length(q) == 1L && is.null(q_mat) && identical(rank, p)) {
    return(amatrix_mlx_qr_qy(q, qty))
  }

  amatrix_mlx_qr_qy(q_mat[, seq_len(rank), drop = FALSE], qty[seq_len(rank), , drop = FALSE])
}

amatrix_mlx_qr_resid <- function(q, y, rank = NULL) {
  y_mat <- as.matrix(y)
  y_mat - amatrix_mlx_qr_fitted(q, y_mat, rank = rank)
}

amatrix_mlx_resident_has <- function(key) {
  .Call("amatrix_mlx_resident_has_bridge", as.character(key))
}

amatrix_mlx_resident_store <- function(key, x) {
  x_mat <- as.matrix(x)
  if (!is.double(x_mat)) {
    storage.mode(x_mat) <- "double"
  }
  invisible(.Call("amatrix_mlx_resident_store_bridge", as.character(key), x_mat))
}

amatrix_mlx_resident_drop <- function(key) {
  invisible(.Call("amatrix_mlx_resident_drop_bridge", as.character(key)))
}

amatrix_mlx_resident_materialize <- function(key) {
  .Call("amatrix_mlx_resident_materialize_bridge", as.character(key))
}

amatrix_mlx_transpose_resident <- function(x_key, out_key) {
  invisible(.Call("amatrix_mlx_transpose_resident_bridge", as.character(x_key), as.character(out_key)))
}

amatrix_mlx_matmul_resident <- function(x_key, y_key, out_key) {
  .Call("amatrix_mlx_matmul_resident_bridge", as.character(x_key), as.character(y_key), as.character(out_key))
}

amatrix_mlx_matmul_resident_host <- function(x_key, y) {
  y_mat <- as.matrix(y)
  if (!is.double(y_mat)) {
    storage.mode(y_mat) <- "double"
  }
  .Call("amatrix_mlx_matmul_resident_host_bridge", as.character(x_key), y_mat)
}

amatrix_mlx_crossprod_resident <- function(x_key, y_key = NULL, out_key) {
  rhs_key <- if (is.null(y_key)) NULL else as.character(y_key)
  .Call("amatrix_mlx_crossprod_resident_bridge", as.character(x_key), rhs_key, as.character(out_key))
}

amatrix_mlx_tcrossprod_resident <- function(x_key, y_key = NULL, out_key) {
  rhs_key <- if (is.null(y_key)) NULL else as.character(y_key)
  .Call("amatrix_mlx_tcrossprod_resident_bridge", as.character(x_key), rhs_key, as.character(out_key))
}

amatrix_mlx_qr_Q_resident <- function(x_key, q_out_key) {
  invisible(.Call("amatrix_mlx_qr_Q_resident_bridge",
                  as.character(x_key), as.character(q_out_key)))
}

# ts_svd: QR on CPU + two GPU matmuls for tall-skinny matrices.
# For A [m×n] with m >> n, QR→SVD(R) reduces the problem to an n×n SVD on CPU:
#   A = Q R  (CPU QR, Q [m×n])
#   SVD(R)   (CPU, R is n×n)
#   U = Q U_R  (GPU matmul, result [m×k])
# This avoids sending a large [m×n] matrix through mlx_linalg_svd and returns
# the back-transform matmul to the GPU where it scales well.
amatrix_mlx_ts_svd <- function(x, nu, nv) {
  mat <- as.matrix(x)
  if (!is.double(mat)) storage.mode(mat) <- "double"
  m <- nrow(mat); n <- ncol(mat); k <- min(m, n)
  # Step 1: thin QR on CPU
  Q <- qr.Q(qr(mat))          # m×k orthonormal
  storage.mode(Q) <- "double"
  # Step 2: B = Q^T A  [k×n] on GPU
  B <- .Call("amatrix_mlx_crossprod_bridge", Q, mat, PACKAGE = "amatrix.mlx")
  # Step 3: exact SVD of small B on CPU
  nu_eff <- min(as.integer(nu), k)
  nv_eff <- min(as.integer(nv), k)
  sv_B <- base::svd(B, nu = nu_eff, nv = nv_eff)
  # Step 4: U = Q * U_B  [m×nu_eff] on GPU
  U_B <- sv_B$u
  storage.mode(U_B) <- "double"
  U <- .Call("amatrix_mlx_matmul_bridge", Q, U_B, PACKAGE = "amatrix.mlx")
  list(u = U, d = sv_B$d, v = sv_B$v)
}

amatrix_mlx_svd <- function(x, nu, nv) {
  x_mat <- as.matrix(x)
  if (!is.double(x_mat)) storage.mode(x_mat) <- "double"
  m <- nrow(x_mat); n <- ncol(x_mat)
  # Route tall-skinny matrices (aspect ratio > 4, narrow dim ≤ 512) to ts_svd.
  # For m >> n, QR→SVD(R) reduces the SVD from O(mn²) to O(n³) on CPU plus
  # two GPU matmuls — significantly cheaper than full mlx_linalg_svd.
  ts_max_n <- getOption("amatrix.mlx.ts_svd_max_n", 512L)
  if (m > 4L * n && n <= ts_max_n) {
    return(amatrix_mlx_ts_svd(x_mat, nu = nu, nv = nv))
  }
  .Call("amatrix_mlx_svd_bridge", x_mat, as.integer(nu), as.integer(nv))
}

amatrix_mlx_rowSums_resident <- function(x_key, na.rm = FALSE, dims = 1L) {
  x_host <- amatrix_mlx_resident_materialize(x_key)
  base::rowSums(x_host, na.rm = na.rm, dims = dims)
}

amatrix_mlx_colSums_resident <- function(x_key, na.rm = FALSE, dims = 1L) {
  x_host <- amatrix_mlx_resident_materialize(x_key)
  base::colSums(x_host, na.rm = na.rm, dims = dims)
}

amatrix_mlx_solve_resident <- function(a_key, b_key = NULL, out_key) {
  a_host <- amatrix_mlx_resident_materialize(a_key)
  result <- if (is.null(b_key)) {
    # Matrix inverse (small p×p): keep on CPU — result is also small.
    base::solve(a_host)
  } else {
    b_host <- amatrix_mlx_resident_materialize(b_key)
    # A x = B solve: use GPU Cholesky when mlx-c is available (vectorised over
    # all columns of B — critical for large q in am_lm_fit / am_ridge_fit).
    tryCatch(
      .Call("amatrix_mlx_chol_solve_bridge", a_host, b_host,
            PACKAGE = "amatrix.mlx"),
      error = function(e) base::solve(a_host, b_host)
    )
  }
  amatrix_mlx_resident_store(out_key, result)
  result
}

amatrix_mlx_chol_resident <- function(x_key, out_key) {
  x_host <- amatrix_mlx_resident_materialize(x_key)
  result <- tryCatch(
    {
      storage.mode(x_host) <- "double"
      R <- .Call("amatrix_mlx_chol_factor_bridge", x_host, PACKAGE = "amatrix.mlx")
      # MLX does not zero the non-triangular half; enforce upper triangular
      # convention to match base::chol and satisfy the amChol validator.
      R[lower.tri(R)] <- 0
      R
    },
    error = function(e) base::chol(x_host)
  )
  amatrix_mlx_resident_store(out_key, result)
  result
}

# GPU-accelerated triangular solve for am_chol_solve: two trsm calls on the
# cached upper-triangular R factor.  Both calls stay on the MLX GPU stream.
amatrix_mlx_chol_solve_factor <- function(R, B) {
  R_mat <- as.matrix(R)
  B_mat <- as.matrix(B)
  storage.mode(R_mat) <- "double"
  storage.mode(B_mat) <- "double"
  # forward substitution: L z = B  where L = R^T
  Rt <- t(R_mat)
  z <- .Call("amatrix_mlx_solve_triangular_bridge", Rt, B_mat, FALSE,
             PACKAGE = "amatrix.mlx")
  # backward substitution: R x = z
  .Call("amatrix_mlx_solve_triangular_bridge", R_mat, z, TRUE,
        PACKAGE = "amatrix.mlx")
}

amatrix_mlx_scatter_mean <- function(x_key, labels, K) {
  .Call("amatrix_mlx_scatter_mean_bridge",
        as.character(x_key), as.integer(labels), as.integer(K),
        PACKAGE = "amatrix.mlx")
}

amatrix_mlx_segment_sum <- function(x_key, labels, K, out_key) {
  .Call("amatrix_mlx_segment_sum_bridge",
        as.character(x_key), as.integer(labels), as.integer(K),
        as.character(out_key), PACKAGE = "amatrix.mlx")
}

amatrix_mlx_segment_mean <- function(x_key, labels, K, out_key) {
  .Call("amatrix_mlx_segment_mean_bridge",
        as.character(x_key), as.integer(labels), as.integer(K),
        as.character(out_key), PACKAGE = "amatrix.mlx")
}

amatrix_mlx_addmm <- function(a_key, b_r, c_r, alpha, beta, out_key) {
  .Call("amatrix_mlx_addmm_bridge",
        as.character(a_key), b_r, c_r,
        as.double(alpha)[1L], as.double(beta)[1L],
        as.character(out_key), PACKAGE = "amatrix.mlx")
}

amatrix_mlx_argreduce <- function(x_key, axis, is_max) {
  .Call("amatrix_mlx_argreduce_bridge",
        as.character(x_key), as.integer(axis), as.logical(is_max),
        PACKAGE = "amatrix.mlx")
}

amatrix_mlx_broadcast_ewise_resident <- function(lhs_key, v, margin, op, out_key) {
  .Call("amatrix_mlx_broadcast_ewise_resident_bridge",
        as.character(lhs_key), as.double(v), as.integer(margin),
        as.character(op), as.character(out_key),
        PACKAGE = "amatrix.mlx")
}

amatrix_mlx_ewise_resident <- function(lhs_key, rhs, op, out_key) {
  rhs_arg <- rhs
  if (is.character(rhs_arg)) {
    rhs_arg <- as.character(rhs_arg)
  } else if (is.numeric(rhs_arg) && length(rhs_arg) == 1L) {
    rhs_arg <- as.double(rhs_arg)
  } else if (!is.null(rhs_arg)) {
    stop("rhs must be NULL, a resident key, or a numeric scalar")
  }

  .Call("amatrix_mlx_ewise_resident_bridge", as.character(lhs_key), rhs_arg, as.character(op), as.character(out_key))
}

.amatrix_mlx_product_thresholds <- function() {
  list(
    matmul_min_dim = getOption("amatrix.mlx.matmul_min_dim", 128L),
    crossprod_min_dim = getOption("amatrix.mlx.crossprod_min_dim", 2048L),
    tcrossprod_min_dim = getOption("amatrix.mlx.tcrossprod_min_dim", 2048L),
    qr_min_dim = getOption("amatrix.mlx.qr_min_dim", 256L)
  )
}

.amatrix_mlx_forced_available <- function() {
  isTRUE(getOption("amatrix.mlx.available", FALSE))
}

.amatrix_mlx_meets_threshold <- function(x, threshold) {
  dims <- dim(x)
  !is.null(dims) && length(dims) == 2L && max(dims) >= threshold
}

amatrix_mlx_rsvd <- function(x, k, n_oversamples = 10L, n_iter = 2L) {
  mat <- if (is(x, "adgeMatrix")) as.matrix(x) else x
  .Call("amatrix_mlx_rsvd_bridge",
        mat,
        as.integer(k),
        as.integer(n_oversamples),
        as.integer(n_iter))
}

amatrix_mlx_eigh <- function(x) {
  mat <- if (is(x, "adgeMatrix") || is(x, "dgeMatrix")) as.matrix(x) else x
  storage.mode(mat) <- "double"
  .Call("amatrix_mlx_eigh_bridge", mat, PACKAGE = "amatrix.mlx")
}

amatrix_mlx_covariance <- function(x, center = TRUE, denom) {
  mat <- as.matrix(x)
  if (!is.double(mat)) storage.mode(mat) <- "double"
  .Call("amatrix_mlx_covariance_bridge",
        mat,
        as.logical(center),
        as.double(denom),
        PACKAGE = "amatrix.mlx")
}

amatrix_mlx_backend <- function() {
  cpu <- amatrix:::.amatrix_cpu_backend()
  capabilities <- amatrix_mlx_capabilities()
  features <- amatrix_mlx_features()
  precision_modes <- amatrix_mlx_precision_modes()
  thresholds <- .amatrix_mlx_product_thresholds()
  sparse_cache <- new.env(parent = emptyenv())

  list(
    capabilities = function() {
      capabilities
    },
    features = function() {
      features
    },
    precision_modes = function() {
      precision_modes
    },
    available = function() {
      amatrix_mlx_is_available()
    },
    supports = function(op, x, y = NULL) {
      # ── Sparse SpMM path ─────────────────────────────────────────────────
      if (is(x, "adgCMatrix")) {
        if (!(op %in% c("matmul", "crossprod", "tcrossprod"))) return(FALSE)
        if (op %in% c("crossprod", "tcrossprod") && is.null(y)) return(FALSE)
        nnz <- length(x@x)
        return(nnz >= getOption("amatrix.mlx.spmm_min_nnz", 10000L))
      }

      if (!is(x, "adgeMatrix") || !(op %in% capabilities)) {
        return(FALSE)
      }

      if (!(x@precision %in% precision_modes)) {
        return(FALSE)
      }

      if (.amatrix_mlx_forced_available()) {
        return(TRUE)
      }

      if (identical(op, "matmul")) {
        return(.amatrix_mlx_meets_threshold(x, thresholds$matmul_min_dim))
      }

      if (identical(op, "crossprod")) {
        return(.amatrix_mlx_meets_threshold(x, thresholds$crossprod_min_dim))
      }

      if (identical(op, "tcrossprod")) {
        return(.amatrix_mlx_meets_threshold(x, thresholds$tcrossprod_min_dim))
      }

      if (identical(op, "qr")) {
        return(.amatrix_mlx_meets_threshold(x, thresholds$qr_min_dim))
      }

      if (identical(op, "eigen")) {
        dims <- dim(x)
        return(!is.null(dims) && length(dims) == 2L && dims[1L] == dims[2L] &&
               dims[1L] >= getOption("amatrix.mlx.eigen_min_dim", 200L))
      }

      if (identical(op, "broadcast_ewise")) {
        return(.amatrix_mlx_meets_threshold(x, thresholds$matmul_min_dim))
      }

      if (identical(op, "argmax")) {
        return(TRUE)
      }

      if (identical(op, "scatter_mean")) {
        return(TRUE)
      }

      if (identical(op, "segment_sum") || identical(op, "segment_mean")) {
        return(TRUE)
      }

      TRUE
    },
    matmul = function(x, y) {
      if (inherits(x, "dgCMatrix"))
        return(amatrix_mlx_spmm(x, y, trans_lhs = FALSE))
      amatrix_mlx_matmul(x, y)
    },
    crossprod = function(x, y = NULL, ...) {
      if (inherits(x, "dgCMatrix")) {
        y_mat <- if (is.matrix(y)) y else as.matrix(y)
        if (!is.double(y_mat)) storage.mode(y_mat) <- "double"
        return(amatrix_mlx_spmm(x, y_mat, trans_lhs = TRUE))
      }
      amatrix_mlx_crossprod(x, y = y)
    },
    tcrossprod = function(x, y = NULL, ...) {
      if (inherits(x, "dgCMatrix")) {
        y_mat <- if (is.matrix(y)) y else as.matrix(y)
        if (!is.double(y_mat)) storage.mode(y_mat) <- "double"
        return(amatrix_mlx_spmm(x, t(y_mat), trans_lhs = FALSE))
      }
      amatrix_mlx_tcrossprod(x, y = y)
    },
    ewise = function(x, lhs, rhs = NULL, op, ...) {
      amatrix_mlx_ewise(lhs = lhs, rhs = rhs, op = op)
    },
    broadcast_ewise = function(x, lhs, v, margin, op, ...) {
      base::sweep(as.matrix(lhs), MARGIN = margin, STATS = v, FUN = op)
    },
    broadcast_ewise_resident = function(lhs_key, v, margin, op, out_key) {
      amatrix_mlx_broadcast_ewise_resident(lhs_key, v, margin, op, out_key)
    },
    scatter_mean_resident = function(x_key, labels, K) {
      amatrix_mlx_scatter_mean(x_key, labels, K)
    },
    segment_sum_resident = function(x_key, labels, K, out_key) {
      amatrix_mlx_segment_sum(x_key, labels, K, out_key)
    },
    segment_mean_resident = function(x_key, labels, K, out_key) {
      amatrix_mlx_segment_mean(x_key, labels, K, out_key)
    },
    addmm_resident = function(a_key, b_r, c_r = NULL, alpha = 1.0, beta = 1.0, out_key) {
      amatrix_mlx_addmm(a_key, b_r, c_r, alpha, beta, out_key)
    },
    rowargmax_resident = function(x_key) amatrix_mlx_argreduce(x_key, 1L, TRUE),
    rowargmin_resident = function(x_key) amatrix_mlx_argreduce(x_key, 1L, FALSE),
    colargmax_resident = function(x_key) amatrix_mlx_argreduce(x_key, 0L, TRUE),
    colargmin_resident = function(x_key) amatrix_mlx_argreduce(x_key, 0L, FALSE),
    rowSums = function(x, na.rm = FALSE, dims = 1L) {
      if (isTRUE(na.rm) || !identical(dims, 1L)) {
        return(cpu$rowSums(x, na.rm = na.rm, dims = dims))
      }
      amatrix_mlx_axis_sums(x, axis = 1L)
    },
    colSums = function(x, na.rm = FALSE, dims = 1L) {
      if (isTRUE(na.rm) || !identical(dims, 1L)) {
        return(cpu$colSums(x, na.rm = na.rm, dims = dims))
      }
      amatrix_mlx_axis_sums(x, axis = 0L)
    },
    qr = function(x, ...) {
      amatrix_mlx_qr(x)
    },
    resident_has = function(key) {
      isTRUE(amatrix_mlx_resident_has(key))
    },
    resident_store = function(key, x) {
      amatrix_mlx_resident_store(key, x)
    },
    resident_drop = function(key) {
      amatrix_mlx_resident_drop(key)
    },
    resident_materialize = function(key) {
      amatrix_mlx_resident_materialize(key)
    },
    transpose_resident = function(x_key, out_key) {
      amatrix_mlx_transpose_resident(x_key, out_key)
    },
    matmul_resident = function(x_key, y_key, out_key) {
      amatrix_mlx_matmul_resident(x_key, y_key, out_key)
    },
    matmul_resident_host = function(x_key, y) {
      amatrix_mlx_matmul_resident_host(x_key, y)
    },
    crossprod_resident = function(x_key, y_key = NULL, out_key) {
      amatrix_mlx_crossprod_resident(x_key, y_key = y_key, out_key = out_key)
    },
    tcrossprod_resident = function(x_key, y_key = NULL, out_key) {
      amatrix_mlx_tcrossprod_resident(x_key, y_key = y_key, out_key = out_key)
    },
    ewise_resident = function(lhs_key, rhs, op, out_key) {
      amatrix_mlx_ewise_resident(lhs_key, rhs, op, out_key)
    },
    rowSums_resident = function(x_key, na.rm = FALSE, dims = 1L) {
      amatrix_mlx_rowSums_resident(x_key, na.rm = na.rm, dims = dims)
    },
    colSums_resident = function(x_key, na.rm = FALSE, dims = 1L) {
      amatrix_mlx_colSums_resident(x_key, na.rm = na.rm, dims = dims)
    },
    solve_resident = function(a_key, b_key = NULL, out_key) {
      amatrix_mlx_solve_resident(a_key, b_key = b_key, out_key = out_key)
    },
    chol_resident = function(x_key, out_key) {
      amatrix_mlx_chol_resident(x_key, out_key = out_key)
    },
    qr_Q_resident = function(x_key, q_out_key) {
      amatrix_mlx_qr_Q_resident(x_key, q_out_key)
    },
    chol_solve_factor = function(R, B) {
      amatrix_mlx_chol_solve_factor(R, B)
    },
    svd = function(x, nu, nv, ...) {
      amatrix_mlx_svd(x, nu = nu, nv = nv)
    },
    rsvd = function(x, k, n_oversamples = 10L, n_iter = 2L, ...) {
      amatrix_mlx_rsvd(x, k = k, n_oversamples = n_oversamples, n_iter = n_iter)
    },
    eigen = function(x, symmetric, only.values = FALSE, EISPACK = FALSE) {
      # x is already a host matrix (materialized by amatrix_dispatch_op)
      if (!isTRUE(symmetric)) {
        return(base::eigen(x, symmetric = FALSE, only.values = only.values))
      }
      res <- amatrix_mlx_eigh(x)
      if (isTRUE(only.values)) {
        list(values = res$values, vectors = NULL)
      } else {
        res
      }
    },
    covariance = function(x, center = TRUE, denom) {
      amatrix_mlx_covariance(x, center = center, denom = denom)
    },
    sparse_resident_store = function(key, x_sp) {
      assign(key, x_sp, envir = sparse_cache)
      invisible(TRUE)
    },
    sparse_resident_has = function(key) {
      exists(key, envir = sparse_cache, inherits = FALSE)
    },
    sparse_resident_drop = function(key) {
      if (exists(key, envir = sparse_cache, inherits = FALSE))
        rm(list = key, envir = sparse_cache)
      invisible(TRUE)
    },
    spmm_resident = function(sp_key, B, trans_lhs = FALSE) {
      x_sp <- get0(sp_key, envir = sparse_cache, inherits = FALSE)
      if (is.null(x_sp)) stop("mlx sparse resident key not found: ", sp_key)
      amatrix_mlx_spmm(x_sp, B, trans_lhs = trans_lhs)
    }
  )
}

amatrix_mlx_register <- function(overwrite = TRUE) {
  amatrix_register_backend("mlx", amatrix_mlx_backend(), overwrite = overwrite)
  invisible("mlx")
}
