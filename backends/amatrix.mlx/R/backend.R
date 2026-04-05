amatrix_mlx_capabilities <- function() {
  c("matmul", "crossprod", "tcrossprod", "ewise", "rowSums", "colSums", "qr", "rsvd")
}

amatrix_mlx_features <- function() {
  c("dense_f32", "resident_dense", "unified_memory", "custom_ops", "qr", "rsvd")
}

amatrix_mlx_precision_modes <- function() {
  "fast"
}

amatrix_mlx_native_available <- function() {
  .Call("amatrix_mlx_native_available_bridge")
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

.amatrix_mlx_qr_explicit <- function(x_mat, include_factor = TRUE) {
  qr_raw <- .amatrix_mlx_qr_explicit_raw(x_mat)
  qr_raw$factor <- if (isTRUE(include_factor)) base::qr(x_mat) else NULL
  qr_raw$factor_source <- if (isTRUE(include_factor)) "bridge_compact" else "bridge_raw"
  qr_raw
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

  qr_raw <- .amatrix_mlx_qr_explicit(x_mat)
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

amatrix_mlx_matmul_resident <- function(x_key, y_key, out_key) {
  .Call("amatrix_mlx_matmul_resident_bridge", as.character(x_key), as.character(y_key), as.character(out_key))
}

amatrix_mlx_crossprod_resident <- function(x_key, y_key = NULL, out_key) {
  rhs_key <- if (is.null(y_key)) NULL else as.character(y_key)
  .Call("amatrix_mlx_crossprod_resident_bridge", as.character(x_key), rhs_key, as.character(out_key))
}

amatrix_mlx_tcrossprod_resident <- function(x_key, y_key = NULL, out_key) {
  rhs_key <- if (is.null(y_key)) NULL else as.character(y_key)
  .Call("amatrix_mlx_tcrossprod_resident_bridge", as.character(x_key), rhs_key, as.character(out_key))
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
    base::solve(a_host)
  } else {
    b_host <- amatrix_mlx_resident_materialize(b_key)
    base::solve(a_host, b_host)
  }
  amatrix_mlx_resident_store(out_key, result)
  result
}

amatrix_mlx_chol_resident <- function(x_key, out_key) {
  x_host <- amatrix_mlx_resident_materialize(x_key)
  result <- base::chol(x_host)
  amatrix_mlx_resident_store(out_key, result)
  result
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

amatrix_mlx_backend <- function() {
  cpu <- amatrix:::.amatrix_cpu_backend()
  capabilities <- amatrix_mlx_capabilities()
  features <- amatrix_mlx_features()
  precision_modes <- amatrix_mlx_precision_modes()
  thresholds <- .amatrix_mlx_product_thresholds()

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

      TRUE
    },
    matmul = function(x, y) {
      amatrix_mlx_matmul(x, y)
    },
    crossprod = function(x, y = NULL, ...) {
      amatrix_mlx_crossprod(x, y = y)
    },
    tcrossprod = function(x, y = NULL, ...) {
      amatrix_mlx_tcrossprod(x, y = y)
    },
    ewise = function(x, lhs, rhs = NULL, op, ...) {
      amatrix_mlx_ewise(lhs = lhs, rhs = rhs, op = op)
    },
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
    matmul_resident = function(x_key, y_key, out_key) {
      amatrix_mlx_matmul_resident(x_key, y_key, out_key)
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
    rsvd = function(x, k, n_oversamples = 10L, n_iter = 2L, ...) {
      amatrix_mlx_rsvd(x, k = k, n_oversamples = n_oversamples, n_iter = n_iter)
    }
  )
}

amatrix_mlx_register <- function(overwrite = TRUE) {
  amatrix_register_backend("mlx", amatrix_mlx_backend(), overwrite = overwrite)
  invisible("mlx")
}
