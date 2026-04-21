 .amatrix_arrayfire_probe_var <- "AMATRIX_ARRAYFIRE_PROBE_GPU"
 .amatrix_arrayfire_state <- new.env(parent = emptyenv())

 .amatrix_arrayfire_probe_enabled <- function() {
   identical(Sys.getenv(.amatrix_arrayfire_probe_var, unset = ""), "1")
 }

 .amatrix_arrayfire_probe_cache_get <- function() {
   get0("native_available", envir = .amatrix_arrayfire_state, inherits = FALSE)
 }

 .amatrix_arrayfire_probe_cache_set <- function(value) {
   assign("native_available", isTRUE(value), envir = .amatrix_arrayfire_state)
   invisible(isTRUE(value))
 }

.amatrix_arrayfire_probe_cache_clear <- function() {
  if (exists("native_available", envir = .amatrix_arrayfire_state, inherits = FALSE)) {
    rm("native_available", envir = .amatrix_arrayfire_state)
  }
  invisible(NULL)
}

.amatrix_arrayfire_running_on_apple_silicon <- function() {
  identical(Sys.info()[["sysname"]], "Darwin") &&
    grepl("arm64|aarch64", R.version$arch, ignore.case = TRUE)
}

.amatrix_arrayfire_configured_runtime_backend <- function() {
  env_backend <- trimws(Sys.getenv("AMATRIX_ARRAYFIRE_BACKEND", unset = ""))
  opt_backend <- trimws(as.character(getOption("amatrix.arrayfire.backend", "")))
  backend <- if (nzchar(env_backend)) {
    env_backend
  } else if (nzchar(opt_backend)) {
    opt_backend
  } else if (.amatrix_arrayfire_running_on_apple_silicon()) {
    # ArrayFire defaults to OpenCL on Apple, which aborts inside clGetDeviceIDs
    # on this host. Pin the runtime to ArrayFire's CPU backend unless the user
    # explicitly requests a different runtime.
    "cpu"
  } else {
    ""
  }

  if (!nzchar(backend)) {
    return(NULL)
  }

  backend <- tolower(backend)
  allowed <- c("cpu", "opencl", "cuda", "oneapi")
  if (!(backend %in% allowed)) {
    warning(
      sprintf(
        "Ignoring unsupported AMATRIX_ARRAYFIRE_BACKEND/amatrix.arrayfire.backend value '%s'",
        backend
      ),
      call. = FALSE,
      immediate. = TRUE
    )
    return(NULL)
  }

  backend
}

.amatrix_arrayfire_configure_runtime_backend <- function(quiet = FALSE) {
  backend <- .amatrix_arrayfire_configured_runtime_backend()
  if (is.null(backend)) {
    return(invisible(NA_character_))
  }

  ok <- tryCatch({
    amatrix_arrayfire_set_backend(backend)
    TRUE
  }, error = function(e) {
    if (!isTRUE(quiet)) {
      warning(
        sprintf("Failed to configure ArrayFire runtime backend '%s': %s", backend, conditionMessage(e)),
        call. = FALSE,
        immediate. = TRUE
      )
    }
    FALSE
  })

  if (isTRUE(ok)) backend else FALSE
}

amatrix_arrayfire_capabilities <- function() {
  c("matmul", "crossprod", "tcrossprod", "ewise", "broadcast_ewise", "argmax", "scatter_mean", "segment_sum", "segment_mean",
    "rowSums", "colSums",
    "qr", "rsvd", "chol", "solve", "covariance", "svd", "kernel_resident")
}

amatrix_arrayfire_features <- function() {
  c("dense_f32", "qr", "op_resident", "rsvd", "chol_gpu", "solve_gpu", "sparse_spmm")
}

amatrix_arrayfire_lapack_available <- function() {
  diag <- amatrix_arrayfire_diagnostics()
  isTRUE(diag$lapack_available)
}

amatrix_arrayfire_precision_modes <- function() {
  "fast"
}

amatrix_arrayfire_native_available <- function(force = FALSE) {
  if (!force) {
    cached <- .amatrix_arrayfire_probe_cache_get()
    if (!is.null(cached)) {
      return(cached)
    }
  }

  if (!.amatrix_arrayfire_probe_enabled()) {
    return(FALSE)
  }

  available <- isTRUE(.Call("amatrix_arrayfire_native_available_bridge"))
  .amatrix_arrayfire_probe_cache_set(available)
}

amatrix_arrayfire_enable_probe <- function(register = TRUE) {
  Sys.setenv(AMATRIX_ARRAYFIRE_PROBE_GPU = "1")
  options(amatrix.enable_arrayfire = TRUE)
  .amatrix_arrayfire_probe_cache_clear()

  available <- amatrix_arrayfire_native_available(force = TRUE)
  if (isTRUE(register)) {
    try(amatrix_arrayfire_register(overwrite = TRUE), silent = TRUE)
  }

  invisible(available)
}

amatrix_arrayfire_is_available <- function() {
  isTRUE(getOption("amatrix.arrayfire.available", FALSE)) || isTRUE(amatrix_arrayfire_native_available())
}

amatrix_arrayfire_bridge_info <- function() {
  info <- .Call("amatrix_arrayfire_bridge_info_bridge")
  info$available <- amatrix_arrayfire_is_available()
  info$capabilities <- amatrix_arrayfire_capabilities()
  info
}

amatrix_arrayfire_diagnostics <- function() {
  .Call("amatrix_arrayfire_diagnostics_bridge")
}

amatrix_arrayfire_active_backend <- function() {
  amatrix_arrayfire_diagnostics()$active_backend
}

amatrix_arrayfire_set_backend <- function(backend = c("cpu", "opencl", "cuda", "oneapi")) {
  backend <- match.arg(backend)
  backend_id <- switch(
    backend,
    cpu = 1L,
    cuda = 2L,
    opencl = 4L,
    oneapi = 8L
  )
  invisible(.Call("amatrix_arrayfire_set_backend_bridge", as.integer(backend_id)))
}

.amatrix_arrayfire_qr_experimental <- function() {
  isTRUE(getOption("amatrix.arrayfire.experimental_qr", FALSE))
}

.amatrix_arrayfire_qr_safe <- function() {
  identical(amatrix_arrayfire_active_backend(), 1L) || .amatrix_arrayfire_qr_experimental()
}

.amatrix_arrayfire_solve_safe <- function() {
  identical(amatrix_arrayfire_active_backend(), 1L) ||
    isTRUE(getOption("amatrix.arrayfire.experimental_solve", FALSE))
}

.amatrix_arrayfire_sparse_safe <- function() {
  identical(amatrix_arrayfire_active_backend(), 1L) ||
    isTRUE(getOption("amatrix.arrayfire.experimental_sparse", FALSE))
}

amatrix_arrayfire_matmul <- function(x, y) {
  x_mat <- as.matrix(x)
  y_mat <- as.matrix(y)

  if (!is.double(x_mat)) {
    storage.mode(x_mat) <- "double"
  }

  if (!is.double(y_mat)) {
    storage.mode(y_mat) <- "double"
  }

  .Call("amatrix_arrayfire_matmul_bridge", x_mat, y_mat)
}

amatrix_arrayfire_crossprod <- function(x, y = NULL) {
  x_mat <- as.matrix(x)
  y_mat <- if (is.null(y)) x_mat else as.matrix(y)

  if (!is.double(x_mat)) storage.mode(x_mat) <- "double"
  if (!is.double(y_mat)) storage.mode(y_mat) <- "double"

  # Use the correct column-major bridge (amatrix_af_from_r + AF_MAT_TRANS, AF_MAT_NONE).
  # The old crossprod_bridge used a row-major staging trick that only works for square
  # matrices; for non-square (n≠p), AF_MAT_TRANS on the y vector causes a dimension
  # mismatch.  The "correct" bridge handles all shapes correctly.
  .Call("amatrix_arrayfire_crossprod_correct_bridge", x_mat, y_mat,
        PACKAGE = "amatrix.arrayfire")
}

amatrix_arrayfire_tcrossprod <- function(x, y = NULL) {
  x_mat <- as.matrix(x)
  y_mat <- if (is.null(y)) x_mat else as.matrix(y)

  if (!is.double(x_mat)) storage.mode(x_mat) <- "double"
  if (!is.double(y_mat)) storage.mode(y_mat) <- "double"

  .Call("amatrix_arrayfire_tcrossprod_correct_bridge", x_mat, y_mat,
        PACKAGE = "amatrix.arrayfire")
}

amatrix_arrayfire_spmm <- function(x_sp, y, trans_lhs = FALSE) {
  # x_sp: dgCMatrix (materialized from adgCMatrix via amatrix_materialize_host)
  # y:    dense host matrix
  # trans_lhs=TRUE: compute t(x_sp) %*% y; FALSE: x_sp %*% y
  y_mat <- if (is.matrix(y)) y else as.matrix(y)
  if (!is.double(y_mat)) storage.mode(y_mat) <- "double"
  .Call("amatrix_arrayfire_spmm_bridge",
        as.double(x_sp@x), as.integer(x_sp@p), as.integer(x_sp@i),
        as.integer(x_sp@Dim), y_mat, as.logical(trans_lhs),
        PACKAGE = "amatrix.arrayfire")
}

amatrix_arrayfire_sparse_store <- function(key, x_sp) {
  # x_sp: dgCMatrix (materialized from adgCMatrix via amatrix_materialize_host)
  invisible(.Call("amatrix_arrayfire_sparse_store_bridge",
                  as.character(key),
                  as.double(x_sp@x), as.integer(x_sp@p), as.integer(x_sp@i),
                  as.integer(x_sp@Dim), NULL,
                  PACKAGE = "amatrix.arrayfire"))
}

amatrix_arrayfire_sparse_has <- function(key) {
  isTRUE(.Call("amatrix_arrayfire_sparse_has_bridge", as.character(key),
               PACKAGE = "amatrix.arrayfire"))
}

amatrix_arrayfire_sparse_drop <- function(key) {
  invisible(.Call("amatrix_arrayfire_sparse_drop_bridge", as.character(key),
                  PACKAGE = "amatrix.arrayfire"))
}

amatrix_arrayfire_spmm_resident <- function(sp_key, B, trans_lhs = FALSE) {
  B_mat <- if (is.matrix(B)) B else as.matrix(B)
  if (!is.double(B_mat)) storage.mode(B_mat) <- "double"
  .Call("amatrix_arrayfire_spmm_resident_bridge",
        as.character(sp_key), B_mat, as.logical(trans_lhs),
        PACKAGE = "amatrix.arrayfire")
}

amatrix_arrayfire_ewise <- function(lhs, rhs = NULL, op) {
  lhs_mat <- as.matrix(lhs)
  rhs_arg <- rhs

  if (!is.double(lhs_mat)) {
    storage.mode(lhs_mat) <- "double"
  }

  if (is.matrix(rhs_arg)) {
    rhs_arg <- as.matrix(rhs_arg)
    if (!is.double(rhs_arg)) {
      storage.mode(rhs_arg) <- "double"
    }
  } else if (is.numeric(rhs_arg) && length(rhs_arg) == 1L) {
    rhs_arg <- as.double(rhs_arg)
  } else if (!is.null(rhs_arg)) {
    stop("rhs must be NULL, a scalar, or a matrix")
  }

  .Call("amatrix_arrayfire_ewise_bridge", lhs_mat, rhs_arg, as.character(op))
}

amatrix_arrayfire_axis_sums <- function(x, axis) {
  x_mat <- as.matrix(x)
  if (!is.double(x_mat)) {
    storage.mode(x_mat) <- "double"
  }
  .Call("amatrix_arrayfire_sum_axis_bridge", x_mat, as.integer(axis))
}

amatrix_arrayfire_qr <- function(x) {
  x_mat <- as.matrix(x)
  if (!is.double(x_mat)) {
    storage.mode(x_mat) <- "double"
  }
  .Call("amatrix_arrayfire_qr_bridge", x_mat)
}

amatrix_arrayfire_resident_store <- function(key, x) {
  x_mat <- as.matrix(x)
  if (!is.double(x_mat)) storage.mode(x_mat) <- "double"
  invisible(.Call("amatrix_arrayfire_resident_store_bridge", as.character(key), x_mat))
}

amatrix_arrayfire_resident_has <- function(key) {
  isTRUE(.Call("amatrix_arrayfire_resident_has_bridge", as.character(key)))
}

amatrix_arrayfire_resident_drop <- function(key) {
  invisible(.Call("amatrix_arrayfire_resident_drop_bridge", as.character(key)))
}

amatrix_arrayfire_resident_materialize <- function(key) {
  .Call("amatrix_arrayfire_resident_materialize_bridge", as.character(key))
}

amatrix_arrayfire_matmul_resident <- function(x_key, y_key, out_key, defer = FALSE) {
  .Call("amatrix_arrayfire_matmul_resident_bridge",
        as.character(x_key), as.character(y_key), as.character(out_key))
  if (defer) NULL else amatrix_arrayfire_resident_materialize(out_key)
}

amatrix_arrayfire_crossprod_resident <- function(x_key, y_key = NULL, out_key, defer = FALSE) {
  rhs_key <- if (is.null(y_key)) NULL else as.character(y_key)
  .Call("amatrix_arrayfire_crossprod_resident_bridge",
        as.character(x_key), rhs_key, as.character(out_key))
  if (defer) NULL else amatrix_arrayfire_resident_materialize(out_key)
}

amatrix_arrayfire_tcrossprod_resident <- function(x_key, y_key = NULL, out_key, defer = FALSE) {
  rhs_key <- if (is.null(y_key)) NULL else as.character(y_key)
  .Call("amatrix_arrayfire_tcrossprod_resident_bridge",
        as.character(x_key), rhs_key, as.character(out_key))
  if (defer) NULL else amatrix_arrayfire_resident_materialize(out_key)
}

amatrix_arrayfire_kernel_resident <- function(out_key, X_mat, Y_mat = NULL,
                                               kernel, sigma, degree, coef,
                                               zero_diag = FALSE) {
  .Call("am_af_kernel_resident_bridge",
        as.character(out_key), X_mat, Y_mat, as.character(kernel),
        as.double(sigma), as.integer(degree), as.double(coef),
        as.logical(zero_diag),
        PACKAGE = "amatrix.arrayfire")
  amatrix_arrayfire_resident_materialize(out_key)
}

amatrix_arrayfire_scatter_mean <- function(x_key, labels, K) {
  .Call("amatrix_arrayfire_scatter_mean_bridge",
        as.character(x_key), as.integer(labels), as.integer(K),
        PACKAGE = "amatrix.arrayfire")
}

amatrix_arrayfire_segment_sum <- function(x_key, labels, K, out_key) {
  .Call("amatrix_arrayfire_segment_sum_bridge",
        as.character(x_key), as.integer(labels), as.integer(K),
        as.character(out_key), PACKAGE = "amatrix.arrayfire")
}

amatrix_arrayfire_segment_mean <- function(x_key, labels, K, out_key) {
  .Call("amatrix_arrayfire_segment_mean_bridge",
        as.character(x_key), as.integer(labels), as.integer(K),
        as.character(out_key), PACKAGE = "amatrix.arrayfire")
}

amatrix_arrayfire_argreduce <- function(x_key, axis, is_max) {
  .Call("amatrix_arrayfire_argreduce_bridge",
        as.character(x_key), as.integer(axis), as.logical(is_max),
        PACKAGE = "amatrix.arrayfire")
}

amatrix_arrayfire_broadcast_ewise_resident <- function(lhs_key, v, margin, op, out_key,
                                                        defer = FALSE) {
  .Call("amatrix_arrayfire_broadcast_ewise_resident_bridge",
        as.character(lhs_key), as.double(v), as.integer(margin),
        as.character(op), as.character(out_key),
        PACKAGE = "amatrix.arrayfire")
  if (defer) NULL else amatrix_arrayfire_resident_materialize(out_key)
}

amatrix_arrayfire_ewise_resident <- function(lhs_key, rhs, op, out_key, defer = FALSE) {
  rhs_arg <- if (is.character(rhs)) as.character(rhs)
             else if (is.numeric(rhs) && length(rhs) == 1L) as.double(rhs)
             else stop("rhs must be a resident key or numeric scalar")
  .Call("amatrix_arrayfire_ewise_resident_bridge",
        as.character(lhs_key), rhs_arg, as.character(op), as.character(out_key))
  if (defer) NULL else amatrix_arrayfire_resident_materialize(out_key)
}

amatrix_arrayfire_rowSums_resident <- function(x_key, na.rm = FALSE, dims = 1L) {
  if (isTRUE(na.rm) || !identical(dims, 1L)) {
    return(base::rowSums(amatrix_arrayfire_resident_materialize(x_key), na.rm = na.rm, dims = dims))
  }
  # Resident arrays use column-major (amatrix_af_from_r): axis=1 sums along columns → rowSums
  .Call("amatrix_arrayfire_sum_axis_resident_bridge", as.character(x_key), 1L)
}

amatrix_arrayfire_colSums_resident <- function(x_key, na.rm = FALSE, dims = 1L) {
  if (isTRUE(na.rm) || !identical(dims, 1L)) {
    return(base::colSums(amatrix_arrayfire_resident_materialize(x_key), na.rm = na.rm, dims = dims))
  }
  # Resident arrays use column-major (amatrix_af_from_r): axis=0 sums along rows → colSums
  .Call("amatrix_arrayfire_sum_axis_resident_bridge", as.character(x_key), 0L)
}

amatrix_arrayfire_chol <- function(x) {
  x_mat <- as.matrix(x)
  if (!is.double(x_mat)) storage.mode(x_mat) <- "double"
  if (amatrix_arrayfire_lapack_available()) {
    .Call("amatrix_arrayfire_chol_bridge", x_mat)
  } else {
    base::chol(x_mat)
  }
}

amatrix_arrayfire_solve <- function(a, b = NULL) {
  a_mat <- as.matrix(a)
  if (!is.double(a_mat)) storage.mode(a_mat) <- "double"
  b_arg <- if (!is.null(b)) { b_mat <- as.matrix(b); if (!is.double(b_mat)) storage.mode(b_mat) <- "double"; b_mat } else NULL
  if (.amatrix_arrayfire_solve_safe() && amatrix_arrayfire_lapack_available()) {
    .Call("amatrix_arrayfire_solve_bridge", a_mat, b_arg)
  } else {
    if (is.null(b_arg)) base::solve(a_mat) else base::solve(a_mat, b_arg)
  }
}

amatrix_arrayfire_covariance <- function(x, center = TRUE, denom = NULL) {
  mat <- as.matrix(x)
  if (!is.double(mat)) storage.mode(mat) <- "double"
  n <- nrow(mat)
  if (isTRUE(center)) {
    col_means <- colMeans(mat)
    mat <- mat - matrix(col_means, nrow = n, ncol = ncol(mat), byrow = TRUE)
  }
  d <- if (is.null(denom)) n - 1L else as.integer(denom)
  amatrix_arrayfire_crossprod(mat) / d
}

amatrix_arrayfire_solve_resident <- function(a_key, b_key = NULL, out_key) {
  if (amatrix_arrayfire_lapack_available()) {
    .Call("amatrix_arrayfire_solve_resident_bridge",
          as.character(a_key),
          if (is.null(b_key)) NULL else as.character(b_key),
          as.character(out_key))
    amatrix_arrayfire_resident_materialize(out_key)
  } else {
    a_host <- amatrix_arrayfire_resident_materialize(a_key)
    result <- if (is.null(b_key)) base::solve(a_host)
              else base::solve(a_host, amatrix_arrayfire_resident_materialize(b_key))
    amatrix_arrayfire_resident_store(out_key, result)
    result
  }
}

amatrix_arrayfire_chol_resident <- function(x_key, out_key) {
  if (amatrix_arrayfire_lapack_available()) {
    .Call("amatrix_arrayfire_chol_resident_bridge",
          as.character(x_key), as.character(out_key))
    amatrix_arrayfire_resident_materialize(out_key)
  } else {
    result <- base::chol(amatrix_arrayfire_resident_materialize(x_key))
    amatrix_arrayfire_resident_store(out_key, result)
    result
  }
}

amatrix_arrayfire_qr_Q_resident <- function(x_key, q_out_key) {
  invisible(.Call("amatrix_arrayfire_qr_Q_resident_bridge",
                  as.character(x_key), as.character(q_out_key)))
}

.amatrix_arrayfire_forced_available <- function() {
  isTRUE(getOption("amatrix.arrayfire.available", FALSE))
}

.amatrix_arrayfire_product_thresholds <- function() {
  list(
    matmul_min_dim = getOption("amatrix.arrayfire.matmul_min_dim", 512L),
    gemv_min_dim = getOption("amatrix.arrayfire.gemv_min_dim", 2048L),
    crossprod_min_dim = getOption("amatrix.arrayfire.crossprod_min_dim", 2048L),
    tcrossprod_min_dim = getOption("amatrix.arrayfire.tcrossprod_min_dim", 2048L),
    ewise_min_dim = getOption("amatrix.arrayfire.ewise_min_dim", 4096L),
    sum_min_dim = getOption("amatrix.arrayfire.sum_min_dim", 4096L),
    qr_min_dim = getOption("amatrix.arrayfire.qr_min_dim", 512L)
  )
}

.amatrix_arrayfire_meets_threshold <- function(x, threshold) {
  dims <- dim(x)
  !is.null(dims) && length(dims) == 2L && max(dims) >= threshold
}

.amatrix_arrayfire_rhs_width <- function(y) {
  if (is.null(y)) {
    return(NA_integer_)
  }

  dims <- dim(y)
  if (is.null(dims)) {
    return(1L)
  }

  if (length(dims) == 1L) {
    return(1L)
  }

  as.integer(dims[[2L]])
}

# ── WY-blocked Householder application ──────────────────────────────────────
#
# Apply Q = H_1 H_2 ... H_nb to C[start_row:(start_row+m_sub-1), ] and return
# the updated C.  Uses the WY representation:  Q = I + Y T Y^T
#
# Arguments:
#   Y_sub    m_sub × nb Householder vectors; col p has zeros in rows 1..p-1
#            (lower-triangular structure from Golub-Kahan storage).
#   tau_sub  length-nb tau values corresponding to H_1 ... H_nb.
#   C        m_full × p matrix — rows start_row..(start_row+m_sub-1) updated.
#   start_row  1-based row offset of Y_sub within C.
#   bdc_min  minimum dimension for GPU GEMM dispatch.
#
# The big GEMM  Y_sub %*% V2  routes to GPU when both m_sub and ncol(C) ≥ bdc_min.
.bdc_wy_apply_left <- function(Y_sub, tau_sub, C, start_row, bdc_min) {
  nb_blk <- length(tau_sub)
  m_sub  <- nrow(Y_sub)
  rows   <- start_row:(start_row + m_sub - 1L)
  if (!any(tau_sub != 0)) return(C)

  # Build T (dlarft column-by-column): H = I + Y T Y^T
  T_mat <- matrix(0, nb_blk, nb_blk)
  for (p in seq_len(nb_blk)) {
    if (tau_sub[p] == 0) next
    T_mat[p, p] <- -tau_sub[p]
    if (p > 1L) {
      prev <- seq_len(p - 1L)
      # Only rows p..m_sub are nonzero in col p (lower-triangular structure)
      T_mat[prev, p] <- -tau_sub[p] *
        T_mat[prev, prev, drop = FALSE] %*%
        crossprod(Y_sub[p:m_sub, prev, drop = FALSE],
                  Y_sub[p:m_sub, p,    drop = FALSE])
    }
  }

  C_sub  <- C[rows, , drop = FALSE]
  V1     <- crossprod(Y_sub, C_sub)   # nb × p  — small
  V2     <- T_mat %*% V1              # nb × p  — small
  # Y_sub %*% V2 : m_sub × p — GPU when large
  update <- if (m_sub >= bdc_min && ncol(C_sub) >= bdc_min)
    .Call("amatrix_arrayfire_matmul_correct_bridge", Y_sub, V2,
          PACKAGE = "amatrix.arrayfire")
  else
    Y_sub %*% V2
  C[rows, ] <- C_sub + update
  C
}

# ── BDC SVD: Golub-Kahan bidiagonalization + CPU BDC + WY back-transform ────
#
# For large square / near-square matrices (aspect ratio < 4) where ts_svd
# (m >> n) gives no advantage.  GPU is used only in Phase 3 (2 GEMMs total):
#   Phase 1: LAPACK dgebrd (blocked, entirely on CPU) — replaces k unblocked
#            GPU round-trips with a single optimised LAPACK call.
#   Phase 2: base::svd on the small k×k bidiagonal B (CPU, O(k²)).
#   Phase 3: dorgbr forms Q (m×k) and P^T (k×n) on CPU, then
#            U = Q %*% sv_B$u and V = t(P^T) %*% sv_B$v via 2 GPU GEMMs.
#
# Returns list(u, d, v) matching base::svd convention.
amatrix_arrayfire_bdc_svd <- function(x, nu, nv) {
  mat <- if (is.matrix(x)) x else as.matrix(x)
  storage.mode(mat) <- "double"
  m <- nrow(mat); n <- ncol(mat); k <- min(m, n)
  nu_eff <- min(as.integer(nu), k)
  nv_eff <- min(as.integer(nv), k)

  # ── Phase 1: LAPACK dgebrd — blocked bidiagonalization, 0 GPU calls ────
  # Returns list(a, d, e, tauq, taup) where a holds packed reflectors.
  # Q^T * mat * P = B  (upper bidiagonal if m >= n, lower if m < n)
  brd <- .Call("amatrix_arrayfire_bdc_bidiag_bridge", mat,
               PACKAGE = "amatrix.arrayfire")

  # ── Phase 2: bidiagonal D&C SVD via dbdsdc ─────────────────────────────
  # dbdsdc operates directly on (d, e) vectors — no k×k dense B matrix needed.
  # This is 5-10× faster than base::svd() which calls dgesdd and treats the
  # bidiagonal matrix as a general dense matrix.
  #
  # uplo convention: dgebrd produces upper bidiagonal (m >= n) or lower (m < n).
  # dbdsdc returns singular values in decreasing order (matching base::svd).
  # sv_bd$vt is V^T (k×k): row j = j-th right singular vector.
  uplo <- if (m >= n) "U" else "L"
  sv_bd <- .Call("amatrix_arrayfire_bdc_dbdsdc_bridge", brd$d, brd$e, uplo,
                 PACKAGE = "amatrix.arrayfire")
  # Alias into the name convention used in Phase 3 (sv_B$u, sv_B$v)
  u_B <- if (nu_eff > 0L) sv_bd$u[, seq_len(nu_eff), drop = FALSE] else matrix(0, k, 0L)
  v_B <- if (nv_eff > 0L) t(sv_bd$vt)[, seq_len(nv_eff), drop = FALSE] else matrix(0, k, 0L)
  d_B <- sv_bd$d

  # ── Phase 3: Form Q and P^T via dorgbr, back-transform ──────────────────
  #
  # LAPACK dorgbr K parameter semantics (critical for wide matrices):
  #   vect="Q": K = number of COLUMNS of the original matrix passed to dgebrd = n
  #   vect="P": K = number of ROWS    of the original matrix passed to dgebrd = m
  #
  # For square/tall (m >= n), k = n so K_Q = n = k (coincides).
  # For wide (m < n), k = m so K_Q = n ≠ k — using k here takes the wrong LAPACK
  # branch (reads one extra unset reflector), producing garbage.
  #
  # dorgbr("Q", M=m, N=k, K=n) → m×k  Q
  # dorgbr("P", M=k, N=n, K=m) → k×n  P^T
  #
  # GPU GEMM threshold: for small matrices (m*nu_eff or n*nv_eff < gemm_min²)
  # fall back to CPU `%*%` to preserve float64 precision; large matrices use
  # GPU (float32, ~1e-7) which is well within the 1e-4 cross-backend tolerance.
  gemm_min <- as.integer(getOption("amatrix.arrayfire.bdc_gemm_min", 256L))

  U_out <- matrix(0, m, 0L)
  if (nu_eff > 0L) {
    Q_brd <- .Call("amatrix_arrayfire_bdc_orgbr_bridge",
                   "Q", brd$a, brd$tauq, m, k, n,   # K = n (cols of original)
                   PACKAGE = "amatrix.arrayfire")
    storage.mode(Q_brd) <- "double"
    if (m * nu_eff >= gemm_min * gemm_min) {
      U_out <- .Call("amatrix_arrayfire_matmul_correct_bridge",
                     Q_brd, u_B, PACKAGE = "amatrix.arrayfire")
    } else {
      U_out <- Q_brd %*% u_B
    }
  }

  V_out <- matrix(0, n, 0L)
  if (nv_eff > 0L) {
    Pt_brd <- .Call("amatrix_arrayfire_bdc_orgbr_bridge",
                    "P", brd$a, brd$taup, k, n, m,   # K = m (rows of original)
                    PACKAGE = "amatrix.arrayfire")
    storage.mode(Pt_brd) <- "double"
    # crossprod(A, B) = t(A) %*% B, so t(P^T) %*% v_B = P %*% v_B
    if (n * nv_eff >= gemm_min * gemm_min) {
      V_out <- .Call("amatrix_arrayfire_crossprod_correct_bridge",
                     Pt_brd, v_B, PACKAGE = "amatrix.arrayfire")
    } else {
      V_out <- crossprod(Pt_brd, v_B)
    }
  }

  list(u = U_out, d = d_B, v = V_out)
}

# Layer 1: R-QR + GPU matmul SVD — always safe.
# af_qr crashes on Metal/OpenCL for matrices >= ~90 rows, so QR is done on CPU.
# Only the two large matmuls (projection B = Q^T A, back-transform U = Q * U_R)
# are offloaded to GPU via the proven-stable correct bridges.
# Pattern mirrors amatrix_arrayfire_rsvd (which already uses CPU QR + GPU matmul).
amatrix_arrayfire_ts_svd <- function(x, nu, nv) {
  mat <- if (is.matrix(x)) x else as.matrix(x)
  storage.mode(mat) <- "double"
  m <- nrow(mat); n <- ncol(mat); k <- min(m, n)
  # Step 1: thin QR on CPU — safe for any size
  Q <- qr.Q(qr(mat))          # m×k orthonormal basis
  storage.mode(Q) <- "double"
  # Step 2: B = Q^T A  [k×n] on GPU (crossprod = t(Q) %*% mat)
  B <- .Call("amatrix_arrayfire_crossprod_correct_bridge", Q, mat,
             PACKAGE = "amatrix.arrayfire")
  # Step 3: exact SVD of B [k×n] on CPU
  nu_eff <- min(as.integer(nu), k)
  nv_eff <- min(as.integer(nv), k)
  sv_B <- base::svd(B, nu = nu_eff, nv = nv_eff)
  # Step 4: U = Q * U_B  [m×nu_eff] on GPU
  U_B <- sv_B$u
  storage.mode(U_B) <- "double"
  U <- .Call("amatrix_arrayfire_matmul_correct_bridge", Q, U_B,
             PACKAGE = "amatrix.arrayfire")
  list(u = U, d = sv_B$d, v = sv_B$v)
}

# ── GESVDA-style SVD (normal equations path) ─────────────────────────────────
#
# For very tall-skinny matrices (m >> n), the GESVDA approach from the NVIDIA
# GTC talk (Chien & Rodriguez Bernabeu) is faster than QR-based ts_svd:
#
#   Phase 1:  B = A^T A  [n×n] via GPU crossprod — pure BLAS-3, fully parallel
#   Phase 2:  (S², V) = eigen(B, symmetric=TRUE)  [n×n] on CPU
#   Phase 3:  U = A V S^{-1}  [m×nu_eff] via GPU matmul
#
# Rationale: CPU QR on a tall matrix (O(mn²) flops, BLAS-1 dominated) becomes
# the bottleneck when m >> n; replacing it with a GPU GEMM (BLAS-3) eliminates
# this. The paper shows 5-17× speedup for m/n ≥ 10.
#
# Accuracy trade-off: forming A^T A squares the condition number. Singular
# values are accurate to ~1e-6 (vs ~1e-15 for QR). Acceptable for PCA /
# data analytics; not suitable for ill-conditioned systems.
#
# Only invoked when:
#   aspect >= amatrix.arrayfire.gesvda_min_aspect  (default 8)
#   AND m  >= amatrix.arrayfire.gesvda_min_rows    (default 512)
amatrix_arrayfire_gesvda_svd <- function(x, nu, nv) {
  mat <- if (is.matrix(x)) x else as.matrix(x)
  storage.mode(mat) <- "double"
  m <- nrow(mat); n <- ncol(mat); k <- min(m, n)
  nu_eff <- min(as.integer(nu), k)
  nv_eff <- min(as.integer(nv), k)

  # Phase 1: B = A^T A [n×n] on GPU — BLAS-3, replaces CPU QR
  B <- .Call("amatrix_arrayfire_crossprod_correct_bridge", mat, mat,
             PACKAGE = "amatrix.arrayfire")
  storage.mode(B) <- "double"

  # Phase 2: eigendecomposition of n×n symmetric B on CPU
  # eigen() returns eigenvalues in *decreasing* order for symmetric matrices,
  # which matches base::svd's convention (singular values decreasing).
  ev <- eigen(B, symmetric = TRUE)
  # sigma = sqrt(lambda); clamp negatives (floating-point noise near zero)
  d_all <- sqrt(pmax(ev$values, 0))   # length n, decreasing
  V_all <- ev$vectors                 # n×n, columns = right singular vectors

  # Phase 3: U = A V S^{-1} [m×nu_eff] on GPU
  # Tolerance: singular values below eps_rel * sigma_1 are numerically zero
  eps_rel <- .Machine$double.eps * max(m, n)
  sigma_thresh <- d_all[1L] * eps_rel

  U_out <- matrix(0, m, 0L)
  if (nu_eff > 0L) {
    V_eff <- V_all[, seq_len(nu_eff), drop = FALSE]
    storage.mode(V_eff) <- "double"
    AV <- .Call("amatrix_arrayfire_matmul_correct_bridge", mat, V_eff,
                PACKAGE = "amatrix.arrayfire")
    # Scale columns: u_j = (A v_j) / sigma_j
    sigma_eff <- d_all[seq_len(nu_eff)]
    safe_inv  <- ifelse(sigma_eff > sigma_thresh, 1 / sigma_eff, 0)
    U_out <- AV %*% diag(safe_inv, nu_eff, nu_eff)
  }

  V_out <- if (nv_eff > 0L) V_all[, seq_len(nv_eff), drop = FALSE] else matrix(0, n, 0L)

  list(u = U_out, d = d_all, v = V_out)
}

# ── Persistent probe/quarantine ───────────────────────────────────────────────
#
# Probe result is keyed by:
#   package version + R arch + OS type + AF lib mtime (catches upgrades)
# Stored in tools::R_user_dir("amatrix","cache")/af-probe.json.
# Falls back to session-only caching if the cache dir is unwritable.

.amatrix_arrayfire_probe_cache_key <- function() {
  ver   <- tryCatch(as.character(utils::packageVersion("amatrix.arrayfire")),
                    error = function(e) "unknown")
  arch  <- R.version$arch
  os    <- .Platform$OS.type
  # Include the AF library mtime so a library upgrade invalidates the cache.
  af_mtime <- tryCatch({
    info <- .Call("amatrix_arrayfire_bridge_info_bridge",
                  PACKAGE = "amatrix.arrayfire")
    lib  <- info$lib_path
    if (!is.null(lib) && nzchar(lib) && file.exists(lib))
      format(file.mtime(lib), "%Y%m%d%H%M%S")
    else
      "nomtime"
  }, error = function(e) "nomtime")
  paste(ver, arch, os, af_mtime, sep = ":")
}

.amatrix_arrayfire_probe_cache_file <- function() {
  d <- tryCatch(
    tools::R_user_dir("amatrix", "cache"),
    error = function(e) NULL
  )
  if (is.null(d)) return(NULL)
  file.path(d, "af-probe.json")
}

.amatrix_arrayfire_probe_read_cache <- function(key) {
  path <- .amatrix_arrayfire_probe_cache_file()
  if (is.null(path) || !file.exists(path)) return(NULL)
  tryCatch({
    entries <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    entries[[key]]   # NULL if key absent
  }, error = function(e) NULL)
}

.amatrix_arrayfire_probe_write_cache <- function(key, value) {
  path <- .amatrix_arrayfire_probe_cache_file()
  if (is.null(path)) return(invisible(NULL))
  tryCatch({
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    entries <- if (file.exists(path)) {
      tryCatch(jsonlite::fromJSON(path, simplifyVector = FALSE),
               error = function(e) list())
    } else {
      list()
    }
    entries[[key]] <- value
    writeLines(jsonlite::toJSON(entries, auto_unbox = TRUE), path)
  }, error = function(e) invisible(NULL))
}

# Layer 2: subprocess probe — tests af_svd(100×100) in a fresh Rscript child.
# Non-zero exit (SIGSEGV → 139, abort → non-zero) → quarantined.
# Result persisted to ~/.cache/amatrix/af-probe.json and also session-cached.
.amatrix_arrayfire_probe_native_svd <- function() {
  # Fast path: session cache
  cached <- getOption("amatrix.arrayfire.native_svd_available")
  if (!is.null(cached)) return(isTRUE(cached))

  # Persistent cache (keyed by env signature)
  key <- .amatrix_arrayfire_probe_cache_key()
  persisted <- .amatrix_arrayfire_probe_read_cache(key)
  if (!is.null(persisted)) {
    options(amatrix.arrayfire.native_svd_available = isTRUE(persisted))
    return(isTRUE(persisted))
  }

  # Run probe subprocess
  runtime_backend <- .amatrix_arrayfire_configured_runtime_backend()
  env_setup <- c(
    'Sys.setenv(AMATRIX_ARRAYFIRE_PROBE_GPU="1")'
  )
  if (!is.null(runtime_backend)) {
    env_setup <- c(
      env_setup,
      sprintf('Sys.setenv(AMATRIX_ARRAYFIRE_BACKEND="%s")', runtime_backend)
    )
  }

  script <- paste0(
    paste(env_setup, collapse = ";"), ";",
    'suppressPackageStartupMessages(library(amatrix.arrayfire));',
    'x <- matrix(rnorm(100*100), 100L, 100L);',
    'storage.mode(x) <- "double";',
    '.Call("amatrix_arrayfire_svd_bridge", x, 5L, 5L,',
    ' PACKAGE="amatrix.arrayfire");',
    'cat("OK\\n")'
  )
  ok <- tryCatch({
    ret <- system2(
      file.path(R.home("bin"), "Rscript"),
      c("--vanilla", "-e", shQuote(script)),
      stdout = FALSE, stderr = FALSE
    )
    identical(ret, 0L)
  }, error = function(e) FALSE)

  .amatrix_arrayfire_probe_write_cache(key, ok)
  options(amatrix.arrayfire.native_svd_available = ok)
  ok
}

# Layer 3 helper: C bridge checks af_get_active_backend(); TRUE for CUDA/oneAPI.
.amatrix_arrayfire_svd_safe_backend <- function() {
  tryCatch(
    isTRUE(.Call("amatrix_arrayfire_svd_safe_bridge",
                 PACKAGE = "amatrix.arrayfire")),
    error = function(e) FALSE
  )
}

# Top-level SVD dispatcher.
#
# Priority:
#   1. Native af_svd when known safe (Layer 3: CUDA/oneAPI, or Layer 2: probe).
#   2. BDC path for large square-ish matrices (aspect < 4, min(m,n) >= bdc_min_n).
#   3. GESVDA for very tall-skinny matrices (aspect >= gesvda_min_aspect AND m >= gesvda_min_rows).
#      Replaces CPU QR with GPU A^T A GEMM + CPU eigen — 5-17× faster when m >> n.
#      Accuracy: ~1e-6 (condition squaring); suitable for PCA / data analytics.
#   4. ts_svd fallback (QR→SVD(R)) for moderate aspect or small matrices.
#
# Shape routing rationale:
#   BDC handles square/near-square where QR offers no dimension reduction.
#   GESVDA beats ts_svd when QR dominates (m/n ≥ 8, paper shows 5-17× speedup).
#   ts_svd remains the safe default for moderate aspect or accuracy-sensitive work.
amatrix_arrayfire_svd <- function(x, nu, nv) {
  use_native <- .amatrix_arrayfire_svd_safe_backend() ||
                .amatrix_arrayfire_probe_native_svd()
  if (use_native) {
    mat <- if (is.matrix(x)) x else as.matrix(x)
    storage.mode(mat) <- "double"
    .Call("amatrix_arrayfire_svd_bridge", mat,
          as.integer(nu), as.integer(nv),
          PACKAGE = "amatrix.arrayfire")
  } else {
    mat    <- if (is.matrix(x)) x else as.matrix(x)
    m      <- nrow(mat); n <- ncol(mat)
    aspect <- max(m, n) / max(min(m, n), 1L)
    bdc_n  <- as.integer(getOption("amatrix.arrayfire.bdc_min_n", 512L))
    gesvda_aspect <- as.numeric(getOption("amatrix.arrayfire.gesvda_min_aspect", 8))
    gesvda_rows   <- as.integer(getOption("amatrix.arrayfire.gesvda_min_rows",  512L))
    if (min(m, n) >= bdc_n && aspect < 4) {
      amatrix_arrayfire_bdc_svd(mat, nu = nu, nv = nv)
    } else if (m >= gesvda_rows && m >= n && aspect >= gesvda_aspect) {
      amatrix_arrayfire_gesvda_svd(mat, nu = nu, nv = nv)
    } else {
      amatrix_arrayfire_ts_svd(mat, nu = nu, nv = nv)
    }
  }
}

amatrix_arrayfire_rsvd <- function(x, k, n_oversamples = 10L, n_iter = 4L) {
  mat <- as.matrix(x)
  storage.mode(mat) <- "double"
  m <- nrow(mat)
  n <- ncol(mat)
  p <- min(k + n_oversamples, m, n)
  k <- min(k, p)
  # Sketch: Y = A * Omega (m x p) on GPU
  Omega <- matrix(stats::rnorm(n * p), n, p)
  storage.mode(Omega) <- "double"
  Y <- .Call("amatrix_arrayfire_matmul_correct_bridge", mat, Omega, PACKAGE = "amatrix.arrayfire")
  # Thin QR on CPU: p is small (<=60), so O(m*p^2) cost is negligible
  Q <- qr.Q(qr(Y))  # m x p orthonormal
  storage.mode(Q) <- "double"
  # Power iteration: big matmuls on GPU, thin QR on CPU
  for (i in seq_len(n_iter)) {
    Z <- .Call("amatrix_arrayfire_crossprod_correct_bridge", mat, Q, PACKAGE = "amatrix.arrayfire")  # n x p
    Q <- qr.Q(qr(Z))  # n x p orthonormal, CPU
    storage.mode(Q) <- "double"
    Z <- .Call("amatrix_arrayfire_matmul_correct_bridge", mat, Q, PACKAGE = "amatrix.arrayfire")  # m x p
    Q <- qr.Q(qr(Z))  # m x p orthonormal, CPU
    storage.mode(Q) <- "double"
  }
  # Project: B = Q^T * A  (p x n) on GPU, small SVD on CPU
  B <- .Call("amatrix_arrayfire_crossprod_correct_bridge", Q, mat, PACKAGE = "amatrix.arrayfire")
  svd_B <- base::svd(B, nu = k, nv = k)
  U_B <- svd_B$u
  storage.mode(U_B) <- "double"
  # Lift U back: U = Q * U_B (m x k) on GPU
  U <- .Call("amatrix_arrayfire_matmul_correct_bridge", Q, U_B, PACKAGE = "amatrix.arrayfire")
  list(u = U, d = svd_B$d[seq_len(k)], v = svd_B$v)
}

amatrix_arrayfire_backend <- function() {
  capabilities <- amatrix_arrayfire_capabilities()
  features <- amatrix_arrayfire_features()
  precision_modes <- amatrix_arrayfire_precision_modes()
  thresholds <- .amatrix_arrayfire_product_thresholds()

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
      amatrix_arrayfire_is_available()
    },
    supports = function(op, x, y = NULL) {
      # ── Sparse SpMM path ─────────────────────────────────────────────────
      if (is(x, "adgCMatrix")) {
        if (!(op %in% c("matmul", "crossprod", "tcrossprod"))) return(FALSE)
        if (!.amatrix_arrayfire_sparse_safe()) return(FALSE)
        # crossprod/tcrossprod without y: result is square dense — CPU sparse
        # BLAS (Matrix pkg) handles these well; no benefit to GPU round-trip.
        if (op %in% c("crossprod", "tcrossprod") && is.null(y)) return(FALSE)
        # Sparse routing stays disabled by default until it is calibrated with
        # shape-aware thresholds.
        nnz <- length(x@x)
        rhs_width <- .amatrix_arrayfire_rhs_width(y)
        min_nnz <- if (!is.na(rhs_width) && rhs_width <= 1L) {
          getOption("amatrix.arrayfire.spmv_min_nnz", Inf)
        } else {
          getOption("amatrix.arrayfire.spmm_min_nnz", Inf)
        }
        return(nnz >= min_nnz)
      }

      if (!is(x, "adgeMatrix") || !(op %in% capabilities)) {
        return(FALSE)
      }

      if (!(x@precision %in% precision_modes)) {
        return(FALSE)
      }

      if (.amatrix_arrayfire_forced_available()) {
        # Bypass size thresholds for GEMM-class ops (testing convenience),
        # but keep LAPACK availability gate for decomposition ops and
        # backend safety gates for ops known to SEGV on non-CPU AF backends.
        if (identical(op, "solve")) {
          return(.amatrix_arrayfire_solve_safe() && amatrix_arrayfire_lapack_available())
        }
        if (identical(op, "chol")) {
          return(amatrix_arrayfire_lapack_available())
        }
        if (identical(op, "qr")) {
          return(.amatrix_arrayfire_qr_safe())
        }
        return(TRUE)
      }

      if (identical(op, "matmul")) {
        rhs_width <- .amatrix_arrayfire_rhs_width(y)
        threshold <- if (!is.na(rhs_width) && rhs_width <= 1L) {
          thresholds$gemv_min_dim
        } else {
          thresholds$matmul_min_dim
        }
        return(.amatrix_arrayfire_meets_threshold(x, threshold))
      }

      if (identical(op, "crossprod")) {
        return(.amatrix_arrayfire_meets_threshold(x, thresholds$crossprod_min_dim))
      }

      if (identical(op, "tcrossprod")) {
        return(.amatrix_arrayfire_meets_threshold(x, thresholds$tcrossprod_min_dim))
      }

      if (identical(op, "ewise")) {
        return(.amatrix_arrayfire_meets_threshold(x, thresholds$ewise_min_dim))
      }

      if (identical(op, "broadcast_ewise")) {
        return(.amatrix_arrayfire_meets_threshold(x, thresholds$ewise_min_dim))
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

      if (op %in% c("rowSums", "colSums")) {
        return(.amatrix_arrayfire_meets_threshold(x, thresholds$sum_min_dim))
      }

      if (identical(op, "qr")) {
        return(.amatrix_arrayfire_qr_safe() && .amatrix_arrayfire_meets_threshold(x, thresholds$qr_min_dim))
      }

      if (identical(op, "rsvd")) {
        return(.amatrix_arrayfire_meets_threshold(x, getOption("amatrix.arrayfire.rsvd_min_dim", 400L)))
      }

      if (identical(op, "solve")) {
        return(.amatrix_arrayfire_solve_safe() &&
               amatrix_arrayfire_lapack_available() &&
               .amatrix_arrayfire_meets_threshold(x, getOption("amatrix.arrayfire.lapack_min_dim", 256L)))
      }

      if (identical(op, "chol")) {
        return(amatrix_arrayfire_lapack_available() &&
               .amatrix_arrayfire_meets_threshold(x, getOption("amatrix.arrayfire.lapack_min_dim", 256L)))
      }

      if (identical(op, "svd")) {
        return(.amatrix_arrayfire_meets_threshold(x, getOption("amatrix.arrayfire.svd_min_dim", 256L)))
      }

      if (identical(op, "covariance")) {
        return(.amatrix_arrayfire_meets_threshold(x, thresholds$crossprod_min_dim))
      }

      FALSE
    },
    matmul = function(x, y) {
      if (inherits(x, "dgCMatrix"))
        return(amatrix_arrayfire_spmm(x, y, trans_lhs = FALSE))
      amatrix_arrayfire_matmul(x, y)
    },
    crossprod = function(x, y = NULL, ...) {
      if (inherits(x, "dgCMatrix")) {
        # y != NULL guaranteed by supports() (y=NULL falls back to CPU)
        y_mat <- if (is.matrix(y)) y else as.matrix(y)
        if (!is.double(y_mat)) storage.mode(y_mat) <- "double"
        return(amatrix_arrayfire_spmm(x, y_mat, trans_lhs = TRUE))
      }
      amatrix_arrayfire_crossprod(x, y = y)
    },
    tcrossprod = function(x, y = NULL, ...) {
      if (inherits(x, "dgCMatrix")) {
        # y != NULL guaranteed by supports(); compute X %*% t(Y) as SpMM(X, t(Y))
        y_mat <- if (is.matrix(y)) y else as.matrix(y)
        if (!is.double(y_mat)) storage.mode(y_mat) <- "double"
        return(amatrix_arrayfire_spmm(x, t(y_mat), trans_lhs = FALSE))
      }
      amatrix_arrayfire_tcrossprod(x, y = y)
    },
    ewise = function(x, lhs, rhs = NULL, op, ...) {
      amatrix_arrayfire_ewise(lhs = lhs, rhs = rhs, op = op)
    },
    broadcast_ewise = function(x, lhs, v, margin, op, ...) {
      base::sweep(as.matrix(lhs), MARGIN = margin, STATS = v, FUN = op)
    },
    broadcast_ewise_resident = function(lhs_key, v, margin, op, out_key, defer = FALSE) {
      amatrix_arrayfire_broadcast_ewise_resident(lhs_key, v, margin, op, out_key,
                                                  defer = defer)
    },
    scatter_mean_resident = function(x_key, labels, K) {
      amatrix_arrayfire_scatter_mean(x_key, labels, K)
    },
    segment_sum_resident = function(x_key, labels, K, out_key) {
      amatrix_arrayfire_segment_sum(x_key, labels, K, out_key)
    },
    segment_mean_resident = function(x_key, labels, K, out_key) {
      amatrix_arrayfire_segment_mean(x_key, labels, K, out_key)
    },
    kernel_resident = function(out_key, X_mat, Y_mat = NULL, kernel,
                               sigma, degree, coef, zero_diag = FALSE) {
      amatrix_arrayfire_kernel_resident(out_key, X_mat, Y_mat, kernel,
                                        sigma, degree, coef, zero_diag)
    },
    rowargmax_resident = function(x_key) amatrix_arrayfire_argreduce(x_key, 1L, TRUE),
    rowargmin_resident = function(x_key) amatrix_arrayfire_argreduce(x_key, 1L, FALSE),
    colargmax_resident = function(x_key) amatrix_arrayfire_argreduce(x_key, 0L, TRUE),
    colargmin_resident = function(x_key) amatrix_arrayfire_argreduce(x_key, 0L, FALSE),
    rowSums = function(x, na.rm = FALSE, dims = 1L) {
      if (isTRUE(na.rm) || !identical(dims, 1L)) {
        return(base::rowSums(as.matrix(x), na.rm = na.rm, dims = dims))
      }
      amatrix_arrayfire_axis_sums(x, axis = 0L)
    },
    colSums = function(x, na.rm = FALSE, dims = 1L) {
      if (isTRUE(na.rm) || !identical(dims, 1L)) {
        return(base::colSums(as.matrix(x), na.rm = na.rm, dims = dims))
      }
      amatrix_arrayfire_axis_sums(x, axis = 1L)
    },
    qr = function(x, ...) {
      if (!.amatrix_arrayfire_qr_safe()) {
        return(base::qr(as.matrix(x), ...))
      }
      amatrix_arrayfire_qr(x)
    },
    resident_has = function(key) {
      amatrix_arrayfire_resident_has(key)
    },
    resident_store = function(key, x) {
      amatrix_arrayfire_resident_store(key, x)
    },
    resident_drop = function(key) {
      amatrix_arrayfire_resident_drop(key)
    },
    resident_materialize = function(key) {
      amatrix_arrayfire_resident_materialize(key)
    },
    matmul_resident = function(x_key, y_key, out_key, defer = FALSE) {
      amatrix_arrayfire_matmul_resident(x_key, y_key, out_key, defer = defer)
    },
    crossprod_resident = function(x_key, y_key = NULL, out_key, defer = FALSE) {
      amatrix_arrayfire_crossprod_resident(x_key, y_key = y_key, out_key = out_key,
                                            defer = defer)
    },
    tcrossprod_resident = function(x_key, y_key = NULL, out_key, defer = FALSE) {
      amatrix_arrayfire_tcrossprod_resident(x_key, y_key = y_key, out_key = out_key,
                                             defer = defer)
    },
    ewise_resident = function(lhs_key, rhs, op, out_key, defer = FALSE) {
      amatrix_arrayfire_ewise_resident(lhs_key, rhs, op, out_key, defer = defer)
    },
    rowSums_resident = function(x_key, na.rm = FALSE, dims = 1L) {
      amatrix_arrayfire_rowSums_resident(x_key, na.rm = na.rm, dims = dims)
    },
    colSums_resident = function(x_key, na.rm = FALSE, dims = 1L) {
      amatrix_arrayfire_colSums_resident(x_key, na.rm = na.rm, dims = dims)
    },
    chol = function(x, ...) {
      amatrix_arrayfire_chol(x)
    },
    solve = function(a, b = NULL, ...) {
      amatrix_arrayfire_solve(a, b = b)
    },
    covariance = function(x, center = TRUE, denom = NULL, ...) {
      amatrix_arrayfire_covariance(x, center = center, denom = denom)
    },
    solve_resident = function(a_key, b_key = NULL, out_key) {
      amatrix_arrayfire_solve_resident(a_key, b_key = b_key, out_key = out_key)
    },
    chol_resident = function(x_key, out_key) {
      amatrix_arrayfire_chol_resident(x_key, out_key)
    },
    qr_Q_resident = function(x_key, q_out_key) {
      amatrix_arrayfire_qr_Q_resident(x_key, q_out_key)
    },
    svd = function(x, nu, nv, ...) {
      amatrix_arrayfire_svd(x, nu = nu, nv = nv)
    },
    rsvd = function(x, k, n_oversamples = 10L, n_iter = 4L) {
      amatrix_arrayfire_rsvd(x, k = k, n_oversamples = n_oversamples, n_iter = n_iter)
    },
    sparse_resident_store = function(key, x_sp) {
      amatrix_arrayfire_sparse_store(key, x_sp)
    },
    sparse_resident_has = function(key) {
      amatrix_arrayfire_sparse_has(key)
    },
    sparse_resident_drop = function(key) {
      amatrix_arrayfire_sparse_drop(key)
    },
    spmm_resident = function(sp_key, B, trans_lhs = FALSE) {
      amatrix_arrayfire_spmm_resident(sp_key, B, trans_lhs = trans_lhs)
    }
  )
}

amatrix_arrayfire_register <- function(overwrite = TRUE) {
  register_backend <- getExportedValue("amatrix", "amatrix_register_backend")
  register_backend("arrayfire", amatrix_arrayfire_backend(), overwrite = overwrite)
  .amatrix_arrayfire_configure_runtime_backend(quiet = FALSE)
  invisible("arrayfire")
}
