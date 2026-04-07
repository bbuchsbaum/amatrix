amatrix_arrayfire_capabilities <- function() {
  c("matmul", "crossprod", "tcrossprod", "ewise", "rowSums", "colSums",
    "qr", "rsvd", "chol", "solve", "covariance", "svd")
}

amatrix_arrayfire_features <- function() {
  c("dense_f32", "qr", "op_resident", "rsvd", "chol_gpu", "solve_gpu")
}

amatrix_arrayfire_lapack_available <- function() {
  diag <- amatrix_arrayfire_diagnostics()
  isTRUE(diag$lapack_available)
}

amatrix_arrayfire_precision_modes <- function() {
  "fast"
}

amatrix_arrayfire_native_available <- function() {
  .Call("amatrix_arrayfire_native_available_bridge")
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

amatrix_arrayfire_matmul_resident <- function(x_key, y_key, out_key) {
  .Call("amatrix_arrayfire_matmul_resident_bridge",
        as.character(x_key), as.character(y_key), as.character(out_key))
  amatrix_arrayfire_resident_materialize(out_key)
}

amatrix_arrayfire_crossprod_resident <- function(x_key, y_key = NULL, out_key) {
  rhs_key <- if (is.null(y_key)) NULL else as.character(y_key)
  .Call("amatrix_arrayfire_crossprod_resident_bridge",
        as.character(x_key), rhs_key, as.character(out_key))
  amatrix_arrayfire_resident_materialize(out_key)
}

amatrix_arrayfire_tcrossprod_resident <- function(x_key, y_key = NULL, out_key) {
  rhs_key <- if (is.null(y_key)) NULL else as.character(y_key)
  .Call("amatrix_arrayfire_tcrossprod_resident_bridge",
        as.character(x_key), rhs_key, as.character(out_key))
  amatrix_arrayfire_resident_materialize(out_key)
}

amatrix_arrayfire_ewise_resident <- function(lhs_key, rhs, op, out_key) {
  rhs_arg <- if (is.character(rhs)) as.character(rhs)
             else if (is.numeric(rhs) && length(rhs) == 1L) as.double(rhs)
             else stop("rhs must be a resident key or numeric scalar")
  .Call("amatrix_arrayfire_ewise_resident_bridge",
        as.character(lhs_key), rhs_arg, as.character(op), as.character(out_key))
  amatrix_arrayfire_resident_materialize(out_key)
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
  if (amatrix_arrayfire_lapack_available()) {
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
# (m >> n) gives no advantage.  GPU is used for:
#   Phase 1: trailing-block updates in bidiagonalization when the block
#            has max(rows,cols) >= bdc_gemm_min (default 256).
#   Phase 3: the final Y_sub %*% V2 GEMM in each WY panel (GPU when
#            m_sub * nu_eff >= bdc_gemm_min²).
#
# Returns list(u, d, v) matching base::svd convention.
amatrix_arrayfire_bdc_svd <- function(x, nu, nv) {
  mat <- if (is.matrix(x)) x else as.matrix(x)
  storage.mode(mat) <- "double"
  m <- nrow(mat); n <- ncol(mat); k <- min(m, n)
  nu_eff <- min(as.integer(nu), k)
  nv_eff <- min(as.integer(nv), k)
  bdc_min <- as.integer(getOption("amatrix.arrayfire.bdc_gemm_min",  256L))
  nb      <- as.integer(getOption("amatrix.arrayfire.bdc_panel_width", 64L))

  # Storage:  left_U[j:m, j] = u_j   right_V[(j+1):n, j] = v_j
  left_U    <- matrix(0, m, k)
  left_tau  <- numeric(k)
  n_right   <- max(k - 1L, 0L)
  right_V   <- matrix(0, n, n_right)
  right_tau <- numeric(n_right)

  A <- mat   # working copy — destroyed during bidiagonalization

  # ── Phase 1: Golub-Kahan bidiagonalization ─────────────────────────────
  for (j in seq_len(k)) {

    # Left Householder: zero A[(j+1):m, j]
    col  <- A[j:m, j]
    nrm2 <- sum(col^2)
    if (nrm2 > .Machine$double.eps^2) {
      alpha <- -sign(col[1L] + (col[1L] == 0L)) * sqrt(nrm2)
      u     <- col; u[1L] <- u[1L] - alpha
      tau   <- 2 / sum(u^2)
      trail <- A[j:m, j:n, drop = FALSE]
      w <- if (max(m - j + 1L, n - j + 1L) >= bdc_min)
        .Call("amatrix_arrayfire_crossprod_correct_bridge",
              matrix(u, ncol = 1L), trail, PACKAGE = "amatrix.arrayfire")
      else crossprod(matrix(u, ncol = 1L), trail)
      A[j:m, j:n]    <- trail - tau * (matrix(u, ncol = 1L) %*% w)
      left_U[j:m, j] <- u
      left_tau[j]     <- tau
    }

    # Right Householder: zero A[j, (j+2):n]
    if (j < n) {
      row  <- A[j, (j + 1L):n]
      nrm2 <- sum(row^2)
      if (nrm2 > .Machine$double.eps^2) {
        alpha <- -sign(row[1L] + (row[1L] == 0L)) * sqrt(nrm2)
        v     <- row; v[1L] <- v[1L] - alpha
        tau   <- 2 / sum(v^2)
        trail <- A[j:m, (j + 1L):n, drop = FALSE]
        w <- if (max(m - j + 1L, n - j) >= bdc_min)
          .Call("amatrix_arrayfire_matmul_correct_bridge",
                trail, matrix(v, ncol = 1L), PACKAGE = "amatrix.arrayfire")
        else trail %*% matrix(v, ncol = 1L)
        A[j:m, (j + 1L):n]     <- trail - tau * (w %*% matrix(v, nrow = 1L))
        right_V[(j + 1L):n, j] <- v
        right_tau[j]            <- tau
      }
    }
  }

  # ── Phase 2: Extract bidiagonal + CPU BDC (O(k²)) ──────────────────────
  B <- matrix(0, k, k)
  diag(B) <- A[cbind(seq_len(k), seq_len(k))]
  if (k > 1L)
    B[cbind(seq_len(k - 1L), seq_len(k - 1L) + 1L)] <-
      A[cbind(seq_len(k - 1L), seq_len(k - 1L) + 1L)]
  # base::svd always returns all k singular values regardless of nu/nv.
  # Truncate $d so its length equals max(nu_eff, nv_eff) — consistent with
  # base::svd behaviour when called with truncated nu/nv on a rectangular input.
  k_out <- max(nu_eff, nv_eff, 1L)
  sv_B  <- base::svd(B, nu = nu_eff, nv = nv_eff)
  sv_B$d <- sv_B$d[seq_len(min(k_out, k))]

  # ── Phase 3: WY-blocked back-transform ─────────────────────────────────
  # U = Q_L * sv_B$u  (H_1..H_k applied last-to-first to [sv_B$u; 0_{m-k}])
  U_out <- matrix(0, m, nu_eff)
  if (nu_eff > 0L) {
    U_out[seq_len(k), ] <- sv_B$u
    for (j_start in rev(seq(1L, k, by = nb))) {
      j_end   <- min(j_start + nb - 1L, k)
      panel   <- j_start:j_end
      Y_sub   <- left_U[j_start:m, panel, drop = FALSE]
      tau_sub <- left_tau[panel]
      U_out   <- .bdc_wy_apply_left(Y_sub, tau_sub, U_out, j_start, bdc_min)
    }
  }

  # V = Q_R * sv_B$v  (G_1..G_{k-1} applied last-to-first to [sv_B$v; 0_{n-k}])
  V_out <- matrix(0, n, nv_eff)
  if (nv_eff > 0L && n_right > 0L) {
    V_out[seq_len(k), ] <- sv_B$v
    for (j_start in rev(seq(1L, n_right, by = nb))) {
      j_end   <- min(j_start + nb - 1L, n_right)
      panel   <- j_start:j_end
      # Right reflectors: v_j stored in right_V[(j+1):n, j]; panel submatrix
      # starts at global row j_start+1 with lower-triangular structure.
      Y_sub   <- right_V[(j_start + 1L):n, panel, drop = FALSE]
      tau_sub <- right_tau[panel]
      V_out   <- .bdc_wy_apply_left(Y_sub, tau_sub, V_out, j_start + 1L, bdc_min)
    }
  }

  list(u = U_out, d = sv_B$d, v = V_out)
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
  script <- paste0(
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
#   3. ts_svd fallback (QR→SVD(R)) for tall-thin or small matrices.
#
# Shape routing rationale:
#   ts_svd is optimal when m >> n (QR reduces the problem cheaply).
#   For square / near-square (aspect < 4), ts_svd offers no reduction;
#   the BDC path amortizes bidiagonalization via GPU trailing-block GEMMs.
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
    if (min(m, n) >= bdc_n && aspect < 4) {
      amatrix_arrayfire_bdc_svd(mat, nu = nu, nv = nv)
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
      if (!is(x, "adgeMatrix") || !(op %in% capabilities)) {
        return(FALSE)
      }

      if (!(x@precision %in% precision_modes)) {
        return(FALSE)
      }

      if (.amatrix_arrayfire_forced_available()) {
        # Bypass size thresholds for GEMM-class ops (testing convenience),
        # but keep LAPACK availability gate for decomposition ops.
        if (op %in% c("chol", "solve")) {
          return(amatrix_arrayfire_lapack_available())
        }
        return(TRUE)
      }

      if (identical(op, "matmul")) {
        return(.amatrix_arrayfire_meets_threshold(x, thresholds$matmul_min_dim))
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

      if (op %in% c("rowSums", "colSums")) {
        return(.amatrix_arrayfire_meets_threshold(x, thresholds$sum_min_dim))
      }

      if (identical(op, "qr")) {
        return(.amatrix_arrayfire_qr_safe() && .amatrix_arrayfire_meets_threshold(x, thresholds$qr_min_dim))
      }

      if (identical(op, "rsvd")) {
        return(.amatrix_arrayfire_meets_threshold(x, getOption("amatrix.arrayfire.rsvd_min_dim", 400L)))
      }

      if (op %in% c("chol", "solve")) {
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
      amatrix_arrayfire_matmul(x, y)
    },
    crossprod = function(x, y = NULL, ...) {
      amatrix_arrayfire_crossprod(x, y = y)
    },
    tcrossprod = function(x, y = NULL, ...) {
      amatrix_arrayfire_tcrossprod(x, y = y)
    },
    ewise = function(x, lhs, rhs = NULL, op, ...) {
      amatrix_arrayfire_ewise(lhs = lhs, rhs = rhs, op = op)
    },
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
    matmul_resident = function(x_key, y_key, out_key) {
      amatrix_arrayfire_matmul_resident(x_key, y_key, out_key)
    },
    crossprod_resident = function(x_key, y_key = NULL, out_key) {
      amatrix_arrayfire_crossprod_resident(x_key, y_key = y_key, out_key = out_key)
    },
    tcrossprod_resident = function(x_key, y_key = NULL, out_key) {
      amatrix_arrayfire_tcrossprod_resident(x_key, y_key = y_key, out_key = out_key)
    },
    ewise_resident = function(lhs_key, rhs, op, out_key) {
      amatrix_arrayfire_ewise_resident(lhs_key, rhs, op, out_key)
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
    }
  )
}

amatrix_arrayfire_register <- function(overwrite = TRUE) {
  amatrix_register_backend("arrayfire", amatrix_arrayfire_backend(), overwrite = overwrite)
  invisible("arrayfire")
}
