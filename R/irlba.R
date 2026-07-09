#' GPU-accelerated truncated SVD via irlba
#'
#' Wraps \code{irlba::irlba()} with an \code{adgeMatrix} input so that every
#' Lanczos matrix-vector product routes through the amatrix GPU dispatch path.
#' The matrix \code{A} is kept resident on device; consecutive matvecs in the
#' Lanczos loop avoid host round-trips.
#'
#' @param A A matrix, \code{adgeMatrix}, or \code{adgCMatrix}. Plain matrices
#'   are coerced via \code{adgeMatrix(A, mode=mode, backend=backend)}.
#' @param nv Number of right singular vectors.
#' @param nu Number of left singular vectors. Defaults to \code{nv}.
#' @param mode Execution mode passed to \code{adgeMatrix()} when coercing.
#'   \code{"fast"} permits float32 and enables GPU routing. Ignored if \code{A}
#'   is already an \code{adgeMatrix} or \code{adgCMatrix}.
#' @param backend Backend name (e.g. \code{"mlx"}, \code{"arrayfire"}).
#'   Ignored if \code{A} is already an amatrix object.
#' @param implementation Lanczos implementation to use. \code{"compat"}
#'   preserves the current \code{irlba::irlba()} wrapper behavior. \code{"block"}
#'   routes to \code{\link{block_lanczos}} for a GEMM-oriented approximation.
#' @param block_size Block size passed to \code{am_block_lanczos()} when
#'   \code{implementation = "block"}. Defaults to a small MLX-friendly block
#'   size derived from the requested rank.
#' @param n_steps Number of block Krylov steps passed to
#'   \code{am_block_lanczos()} when \code{implementation = "block"}.
#' @param ... Additional arguments forwarded to \code{irlba::irlba()}.
#'   \code{fastpath} is always forced to \code{FALSE} — the C fastpath bypasses
#'   S4 dispatch and cannot be GPU-accelerated.
#'
#' @return Same structure as \code{irlba::irlba()}: a list with components
#'   \code{d}, \code{u}, \code{v}, \code{iter}, \code{mprod}.
#'
#' @details
#' The hot loop in irlba is two matrix-vector products per Lanczos step:
#' \code{A \%*\% v} and \code{w \%*\% A}. Both route through \code{am_matmul()}
#' when \code{A} is an \code{adgeMatrix}, giving GPU acceleration on the dominant
#' cost. Orthogonalization, \code{svd(B)}, and convergence tests remain on CPU
#' where they belong (the subspace dimension \code{work} is always small).
#'
#' Do \strong{not} pass \code{mult=} — it is deprecated in irlba and forces a
#' non-standard dispatch path. Pass an \code{adgeMatrix} instead.
#'
#' @seealso \code{\link{adgeMatrix}}, \code{\link{svd_factor}}
#' @export
irlba <- function(A,
                     nv = 5,
                     nu = nv,
                     ...,
                     mode = "fast",
                     backend = NULL,
                     implementation = c("compat", "block"),
                     block_size = NULL,
                     n_steps = NULL) {
  implementation <- match.arg(implementation)

  if (identical(implementation, "block")) {
    dots <- list(...)
    if (length(dots) > 0L) {
      warning(
        "irlba(..., implementation = 'block') ignores irlba-specific arguments in ...; ",
        "use block_lanczos() for the block implementation surface.",
        call. = FALSE
      )
    }
    return(block_lanczos(
      A,
      nv = nv,
      nu = nu,
      block_size = block_size,
      n_steps = n_steps,
      mode = mode,
      backend = backend
    ))
  }

  if (!requireNamespace("irlba", quietly = TRUE)) {
    stop("Package 'irlba' is required: install.packages('irlba')")
  }
  if (!inherits(A, c("adgeMatrix", "adgCMatrix"))) {
    A <- adgeMatrix(A, mode = mode, backend = backend)
  }
  irlba::irlba(A, nv = nv, nu = nu, ..., fastpath = FALSE)
}

#' GPU-native truncated SVD via Lanczos bidiagonalization
#'
#' Implements Golub-Kahan Lanczos bidiagonalization directly in ArrayFire C,
#' keeping all matvecs and CGS2 reorthogonalization on the GPU. Only 2*work
#' scalars and the final basis matrices cross PCIe per restart; no per-step
#' host transfers.
#'
#' Compared to \code{am_irlba}, which routes each Lanczos matvec through S4
#' dispatch, this function:
#' \itemize{
#'   \item eliminates S4 overhead on the hot path
#'   \item replaces k sequential GEMVs for reorthogonalization with one GEMM
#'   \item uploads A once and never re-uploads it across restarts
#' }
#'
#' @param A A matrix or \code{adgeMatrix}. Coerced if necessary.
#' @param nv Number of singular values/vectors to compute.
#' @param nu Number of left singular vectors (default = \code{nv}).
#' @param tol Convergence tolerance.
#' @param maxit Maximum number of restarts.
#' @param work Size of the Lanczos subspace per restart.  Larger values
#'   converge in fewer restarts at the cost of more memory and work per
#'   restart. Default is \code{max(nv + 20L, 3L * nv)}.
#' @param v0 Optional starting vector (length \code{ncol(A)}).
#' @param mode,backend Passed to \code{adgeMatrix()} when coercing.
#'
#' @return A list with components \code{d}, \code{u}, \code{v}, \code{iter},
#'   \code{mprod}, compatible with \code{irlba::irlba()}.
#'
#' @seealso \code{\link{irlba}}, \code{\link{adgeMatrix}}
#' @export
irlba_native <- function(A,
                             nv    = 5L,
                             nu    = nv,
                             tol   = sqrt(.Machine$double.eps),
                             maxit = 100L,
                             work  = max(nv + 20L, 3L * nv),
                             v0    = NULL,
                             mode  = "fast",
                             backend = NULL) {

  if (!inherits(A, c("adgeMatrix", "adgCMatrix"))) {
    A <- adgeMatrix(A, mode = mode, backend = backend)
  }

  af_ok <- tryCatch({
    requireNamespace("amatrix.arrayfire", quietly = TRUE) &&
      amatrix.arrayfire::amatrix_arrayfire_is_available()
  }, error = function(e) FALSE)

  if (!af_ok) {
    message("irlba_native: ArrayFire unavailable, falling back to irlba")
    return(irlba(A, nv = nv, nu = nu, tol = tol, maxit = maxit))
  }

  A_mat <- amatrix_materialize_host(A)
  if (inherits(A_mat, "dgeMatrix")) A_mat <- as.matrix(A_mat)
  if (!is.matrix(A_mat) || !is.double(A_mat)) storage.mode(A_mat) <- "double"

  m    <- nrow(A_mat)
  n    <- ncol(A_mat)
  work <- min(as.integer(work), m, n)
  nv   <- min(as.integer(nv), work - 1L)
  if (nv < 1L) stop("nv must be >= 1 after clamping to work size", call. = FALSE)
  nu   <- min(as.integer(nu), nv)

  if (is.null(v0)) v0 <- rnorm(n)
  v0 <- v0 / sqrt(sum(v0^2))

  # Upload A to GPU once; it stays resident across all restarts.
  getExportedValue("amatrix.arrayfire", "amatrix_arrayfire_lanczos_upload_bridge")(A_mat)
  on.exit(
    getExportedValue("amatrix.arrayfire", "amatrix_arrayfire_lanczos_drop_bridge")(),
    add = TRUE
  )

  mprod <- 0L

  # Thick-restart state (NULL until first restart completes)
  V_warm     <- NULL   # n × nv  right Ritz vectors from last restart
  U_warm     <- NULL   # m × nv  left  Ritz vectors from last restart
  sigma_warm <- NULL   # nv      singular values from last restart
  foot       <- NULL   # nv      foot correction for next B_full col (nv+1)

  v_start <- v0
  raw     <- NULL
  svd_B   <- NULL

  for (iter in seq_len(maxit)) {

    nv_warm <- if (is.null(V_warm)) 0L else ncol(V_warm)
    k_new   <- work - nv_warm   # new Lanczos steps this restart

    # ── GPU Lanczos (warm or cold) ─────────────────────────────────
    raw <- getExportedValue("amatrix.arrayfire", "amatrix_arrayfire_lanczos_warm_bridge")(
      V_warm,               # n × nv_warm or NULL
      U_warm,               # m × nv_warm or NULL
      as.double(v_start),   # n-vector starting direction
      as.integer(k_new)
    )
    mprod <- mprod + 2L * k_new

    alpha_new <- raw$alpha   # k_new-vector (new diagonal entries)
    beta_new  <- raw$beta    # k_new-vector (last = residual norm)
    U_full    <- raw$U       # m × work
    V_full    <- raw$V       # n × (work + 1)

    # ── Assemble B_full (work × work) ─────────────────────────────
    # Top-left nv_warm × nv_warm : diag(sigma_warm)
    # Column nv_warm+1, rows 1:nv_warm : foot correction
    # Remaining lower-right block  : new upper bidiagonal
    B_full <- matrix(0.0, nrow = work, ncol = work)
    if (nv_warm > 0L) {
      diag(B_full)[seq_len(nv_warm)] <- sigma_warm
      B_full[seq_len(nv_warm), nv_warm + 1L] <- foot
    }
    for (i in seq_len(k_new)) {
      B_full[nv_warm + i, nv_warm + i] <- alpha_new[i]
      if (i < k_new)
        B_full[nv_warm + i, nv_warm + i + 1L] <- beta_new[i]
    }

    svd_B <- base::svd(B_full)

    # ── Convergence test ───────────────────────────────────────────
    # |beta_residual * Q[work, j]| <= tol * sigma_j  for j = 1..nv
    residuals <- abs(beta_new[k_new] * svd_B$v[work, seq_len(nv)])
    sigma_t   <- svd_B$d[seq_len(nv)]
    ref       <- tol * pmax(sigma_t, .Machine$double.eps)

    if (all(residuals <= ref)) break

    if (iter == maxit) {
      warning("irlba_native: did not converge after ", maxit, " restarts")
      break
    }

    # ── Thick restart: retain nv Ritz vectors ─────────────────────
    V_work <- V_full[, seq_len(work), drop = FALSE]   # n × work
    U_work <- U_full                                   # m × work

    # Ritz basis for next restart
    V_warm     <- V_work %*% svd_B$v[, seq_len(nv), drop = FALSE]  # n × nv
    U_warm     <- U_work %*% svd_B$u[, seq_len(nv), drop = FALSE]  # m × nv
    sigma_warm <- svd_B$d[seq_len(nv)]

    # Foot correction: residual norm × last row of Q restricted to nv cols.
    # This populates B_full[1:nv, nv+1] in the next restart.
    foot <- beta_new[k_new] * svd_B$v[work, seq_len(nv)]

    # New starting direction: (nv+1)-th right Ritz direction
    v_raw <- as.vector(V_work %*% svd_B$v[, nv + 1L, drop = FALSE])
    nrm   <- sqrt(sum(v_raw^2))
    v_start <- if (nrm > 1e-14) v_raw / nrm else {
      tmp <- rnorm(n); tmp / sqrt(sum(tmp^2))
    }
  }

  # ── Form final Ritz approximations ────────────────────────────
  V_work <- V_full[, seq_len(work), drop = FALSE]
  U_work <- U_full

  idx   <- seq_len(nv)
  d_out <- svd_B$d[idx]
  u_out <- U_work %*% svd_B$u[, idx, drop = FALSE]
  v_out <- V_work %*% svd_B$v[, idx, drop = FALSE]

  list(
    d     = d_out,
    u     = u_out[, seq_len(min(nu, nv)), drop = FALSE],
    v     = v_out,
    iter  = iter,
    mprod = mprod
  )
}

#' @noRd
.amatrix_block_lanczos_default_block_size <- function(k) {
  k <- as.integer(k)
  default_block <- if (k >= 16L) k else k + 1L
  max(8L, min(24L, default_block))
}

.amatrix_block_lanczos_default_steps <- function(k, block_size) {
  max(4L, ceiling(as.integer(k) / as.integer(block_size)) + 2L)
}

.amatrix_block_reorth <- function(z, basis, return_projection = FALSE) {
  if (is.null(basis) || ncol(basis) == 0L) {
    if (isTRUE(return_projection)) {
      return(list(z = z, coeff = NULL))
    }
    return(z)
  }

  z_mat <- if (is.matrix(z) && is.double(z)) z else {
    z_out <- as.matrix(z)
    if (!is.double(z_out)) {
      storage.mode(z_out) <- "double"
    }
    z_out
  }
  basis_mat <- if (is.matrix(basis) && is.double(basis)) basis else {
    basis_out <- as.matrix(basis)
    if (!is.double(basis_out)) {
      storage.mode(basis_out) <- "double"
    }
    basis_out
  }

  compiled <- tryCatch(
    .Call(
      "amatrix_block_reorth_bridge",
      z_mat,
      basis_mat,
      as.logical(return_projection)
    ),
    error = function(e) NULL
  )
  if (!is.null(compiled)) {
    return(compiled)
  }

  z_norm <- norm(z_mat, type = "F")
  coeff <- base::crossprod(basis_mat, z_mat)
  z_mat <- z_mat - basis_mat %*% coeff

  # On the current block path, a 0.5 norm-retention trigger trims late
  # redundant passes while preserving the observed singular-value envelope.
  if (is.finite(z_norm) && z_norm > 0) {
    z_reorth_norm <- norm(z_mat, type = "F")
    if (z_reorth_norm <= 0.5 * z_norm) {
      coeff2 <- base::crossprod(basis_mat, z_mat)
      z_mat <- z_mat - basis_mat %*% coeff2
      coeff <- coeff + coeff2
    }
  }

  if (isTRUE(return_projection)) {
    return(list(z = z_mat, coeff = coeff))
  }

  z_mat
}

.amatrix_block_reorth_prefix <- function(z, basis, basis_cols, return_projection = FALSE) {
  basis_cols <- as.integer(basis_cols)
  if (is.na(basis_cols) || basis_cols < 0L) {
    stop("basis_cols must be a non-negative integer", call. = FALSE)
  }

  if (basis_cols == 0L) {
    if (isTRUE(return_projection)) {
      return(list(z = z, coeff = NULL))
    }
    return(z)
  }

  z_mat <- if (is.matrix(z) && is.double(z)) z else {
    z_out <- as.matrix(z)
    if (!is.double(z_out)) {
      storage.mode(z_out) <- "double"
    }
    z_out
  }
  basis_mat <- if (is.matrix(basis) && is.double(basis)) basis else {
    basis_out <- as.matrix(basis)
    if (!is.double(basis_out)) {
      storage.mode(basis_out) <- "double"
    }
    basis_out
  }

  compiled <- tryCatch(
    .Call(
      "amatrix_block_reorth_prefix_bridge",
      z_mat,
      basis_mat,
      basis_cols,
      as.logical(return_projection)
    ),
    error = function(e) NULL
  )
  if (!is.null(compiled)) {
    return(compiled)
  }

  .amatrix_block_reorth(
    z_mat,
    basis_mat[, seq_len(basis_cols), drop = FALSE],
    return_projection = return_projection
  )
}

.amatrix_block_basis_needs_final_qr <- function(q_basis, tol = 1e-8) {
  if (is.null(q_basis) || ncol(q_basis) == 0L) {
    return(FALSE)
  }

  gram <- base::crossprod(q_basis)
  max(abs(gram - diag(ncol(q_basis)))) > tol
}

.amatrix_block_thin_qr <- function(z) {
  z_mat <- if (is.matrix(z) && is.double(z)) z else {
    z_out <- as.matrix(z)
    if (!is.double(z_out)) {
      storage.mode(z_out) <- "double"
    }
    z_out
  }

  compiled <- tryCatch(
    .Call("amatrix_block_thin_qr_bridge", z_mat),
    error = function(e) NULL
  )
  if (!is.null(compiled)) {
    return(compiled)
  }

  qr_obj <- qr(z_mat)
  list(
    q = qr.Q(qr_obj),
    r = qr.R(qr_obj, complete = FALSE)
  )
}

.amatrix_block_lanczos_source_operator <- function(A) {
  if (!inherits(A, "adgeMatrix")) {
    return(.amatrix_block_lanczos_compile_operator(A, op = "matmul"))
  }
  if (!identical(A@preferred_backend, "mlx") || !identical(A@precision, "fast")) {
    return(.amatrix_block_lanczos_compile_operator(A, op = "matmul"))
  }

  backend <- tryCatch(
    .amatrix_get_backend("mlx"),
    error = function(e) {
      if (!requireNamespace("amatrix.mlx", quietly = TRUE)) {
        return(NULL)
      }
      getExportedValue("amatrix.mlx", "amatrix_mlx_register")(overwrite = TRUE)
      tryCatch(.amatrix_get_backend("mlx"), error = function(e2) NULL)
    }
  )
  if (is.null(backend)) {
    return(NULL)
  }
  if (!.amatrix_backend_residency_capable(backend) ||
      !.amatrix_backend_supports_resident_op(backend, "matmul")) {
    return(NULL)
  }

  source_arg <- .amatrix_prepare_resident_arg(A, "mlx")
  if (is.null(source_arg)) {
    return(NULL)
  }

  list(
    backend = "mlx",
    resident_key = source_arg$key,
    temporary = isTRUE(source_arg$temporary)
  )
}

.amatrix_block_lanczos_numeric_matrix <- function(x) {
  x_mat <- if (is.matrix(x) && is.double(x)) x else {
    x_out <- as.matrix(x)
    if (!is.double(x_out)) {
      storage.mode(x_out) <- "double"
    }
    x_out
  }
  x_mat
}

.amatrix_block_lanczos_compile_operator <- function(A, op) {
  if (!inherits(A, "aMatrix")) {
    return(NULL)
  }

  backend_name <- tryCatch(
    amatrix_resident_backend_for(A, op = op),
    error = function(e) NULL
  )
  if (is.null(backend_name) || identical(backend_name, "cpu")) {
    return(NULL)
  }

  tryCatch(
    amatrix_compile_product(A, op = op, backend = backend_name),
    error = function(e) NULL
  )
}

.amatrix_block_lanczos_right_operator <- function(A, source_operator = NULL) {
  if (!inherits(A, "adgeMatrix")) {
    return(.amatrix_block_lanczos_compile_operator(A, op = "crossprod"))
  }
  if (!identical(A@preferred_backend, "mlx") || !identical(A@precision, "fast")) {
    return(.amatrix_block_lanczos_compile_operator(A, op = "crossprod"))
  }

  backend <- tryCatch(
    .amatrix_get_backend("mlx"),
    error = function(e) {
      if (!requireNamespace("amatrix.mlx", quietly = TRUE)) {
        return(NULL)
      }
      getExportedValue("amatrix.mlx", "amatrix_mlx_register")(overwrite = TRUE)
      tryCatch(.amatrix_get_backend("mlx"), error = function(e2) NULL)
    }
  )
  if (is.null(backend)) {
    return(NULL)
  }
  if (!.amatrix_backend_residency_capable(backend) ||
      !.amatrix_backend_supports_resident_op(backend, "transpose") ||
      !.amatrix_backend_supports_resident_op(backend, "matmul")) {
    return(NULL)
  }

  source_key <- if (is.list(source_operator)) source_operator$resident_key else NULL
  if (is.null(source_key) || !nzchar(source_key)) {
    source_arg <- .amatrix_prepare_resident_arg(A, "mlx")
    if (is.null(source_arg)) {
      return(NULL)
    }
    source_key <- source_arg$key
  }

  out_key <- .amatrix_next_resident_key("mlx")
  success <- tryCatch(
    {
      backend$transpose_resident(source_key, out_key)
      TRUE
    },
    error = function(e) FALSE
  )
  if (!success) {
    if (isTRUE(backend$resident_has(out_key))) {
      backend$resident_drop(out_key)
    }
    return(NULL)
  }

  list(
    backend = "mlx",
    resident_key = out_key
  )
}

.amatrix_block_lanczos_drop_source_operator <- function(operator) {
  if (inherits(operator, "am_product_plan")) {
    return(.amatrix_release_product_plan(operator))
  }

  resident_key <- if (is.list(operator)) operator$resident_key else NULL
  if (is.null(operator) || !is.list(operator) || is.null(resident_key) || !nzchar(resident_key) ||
      !isTRUE(operator$temporary)) {
    return(invisible(FALSE))
  }

  backend <- tryCatch(
    .amatrix_get_backend(operator$backend),
    error = function(e) NULL
  )
  if (is.null(backend) || !.amatrix_backend_residency_capable(backend)) {
    return(invisible(FALSE))
  }

  if (isTRUE(backend$resident_has(resident_key))) {
    backend$resident_drop(resident_key)
    return(invisible(TRUE))
  }

  invisible(FALSE)
}

.amatrix_block_lanczos_drop_right_operator <- function(operator) {
  if (inherits(operator, "am_product_plan")) {
    return(.amatrix_release_product_plan(operator))
  }

  resident_key <- if (is.list(operator)) operator$resident_key else NULL
  if (is.null(operator) || !is.list(operator) || is.null(resident_key) || !nzchar(resident_key)) {
    return(invisible(FALSE))
  }

  backend <- tryCatch(
    .amatrix_get_backend(operator$backend),
    error = function(e) NULL
  )
  if (is.null(backend) || !.amatrix_backend_residency_capable(backend)) {
    return(invisible(FALSE))
  }

  if (isTRUE(backend$resident_has(operator$resident_key))) {
    backend$resident_drop(operator$resident_key)
    return(invisible(TRUE))
  }

  invisible(FALSE)
}

.amatrix_block_lanczos_right_product <- function(A, q_block, operator = NULL) {
  if (inherits(operator, "am_product_plan")) {
    return(operator(q_block, materialize = "matrix"))
  }
  if (is.null(operator)) {
    return(as.matrix(crossprod(A, q_block)))
  }
  if (inherits(operator, "adgeMatrix")) {
    return(as.matrix(operator %*% q_block))
  }
  if (!is.list(operator) || is.null(operator$backend) || is.null(operator$resident_key)) {
    return(as.matrix(crossprod(A, q_block)))
  }

  backend_name <- operator$backend
  backend <- tryCatch(
    .amatrix_get_backend(backend_name),
    error = function(e) NULL
  )
  if (is.null(backend)) {
    return(as.matrix(crossprod(A, q_block)))
  }
  if (!.amatrix_backend_supports_resident_op(backend, "matmul")) {
    return(as.matrix(crossprod(A, q_block)))
  }
  q_mat <- .amatrix_block_lanczos_numeric_matrix(q_block)
  if (is.function(backend$matmul_resident_host)) {
    value <- tryCatch(
      backend$matmul_resident_host(operator$resident_key, q_mat),
      error = function(e) NULL
    )
    if (!is.null(value)) {
      return(as.matrix(value))
    }
  }

  rhs <- .amatrix_prepare_resident_arg(q_block, backend_name)
  if (is.null(rhs)) {
    return(as.matrix(crossprod(A, q_mat)))
  }
  on.exit(.amatrix_cleanup_temp_resident(list(rhs), backend_name), add = TRUE)

  out_key <- .amatrix_next_resident_key(backend_name)
  on.exit(
    {
      if (isTRUE(backend$resident_has(out_key))) {
        backend$resident_drop(out_key)
      }
    },
    add = TRUE
  )

  value <- tryCatch(
    backend$matmul_resident(operator$resident_key, rhs$key, out_key),
    error = function(e) NULL
  )
  if (is.null(value)) {
    return(as.matrix(crossprod(A, q_mat)))
  }

  as.matrix(.amatrix_host_arg(value))
}

.amatrix_block_lanczos_left_product <- function(A, q_block, operator = NULL, backend_name = NULL) {
  if (inherits(operator, "am_product_plan")) {
    return(operator(q_block, materialize = "matrix"))
  }
  if (is.list(operator) && !is.null(operator$backend) && !is.null(operator$resident_key)) {
    backend_name <- operator$backend
  }
  if (is.null(backend_name)) {
    if (!inherits(A, "adgeMatrix") ||
        !identical(A@preferred_backend, "mlx") ||
        !identical(A@precision, "fast")) {
      return(as.matrix(A %*% q_block))
    }
    backend_name <- "mlx"
  }

  backend <- tryCatch(
    .amatrix_get_backend(backend_name),
    error = function(e) NULL
  )
  if (is.null(backend) || !.amatrix_backend_supports_resident_op(backend, "matmul")) {
    return(as.matrix(A %*% q_block))
  }
  q_mat <- .amatrix_block_lanczos_numeric_matrix(q_block)
  if (is.list(operator) && !is.null(operator$resident_key) && nzchar(operator$resident_key) &&
      is.function(backend$matmul_resident_host)) {
    value <- tryCatch(
      backend$matmul_resident_host(operator$resident_key, q_mat),
      error = function(e) NULL
    )
    if (!is.null(value)) {
      return(as.matrix(value))
    }
  }

  lhs <- if (is.list(operator) && !is.null(operator$resident_key) && nzchar(operator$resident_key)) {
    list(key = operator$resident_key, temporary = FALSE, tracked = !isTRUE(operator$temporary))
  } else {
    .amatrix_prepare_resident_arg(A, backend_name)
  }
  rhs <- .amatrix_prepare_resident_arg(q_mat, backend_name)
  if (is.null(lhs) || is.null(rhs)) {
    .amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name)
    return(as.matrix(A %*% q_mat))
  }
  on.exit(.amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name), add = TRUE)

  out_key <- .amatrix_next_resident_key(backend_name)
  on.exit(
    {
      if (isTRUE(backend$resident_has(out_key))) {
        backend$resident_drop(out_key)
      }
    },
    add = TRUE
  )

  value <- tryCatch(
    backend$matmul_resident(lhs$key, rhs$key, out_key),
    error = function(e) NULL
  )
  if (is.null(value)) {
    return(as.matrix(A %*% q_mat))
  }

  as.matrix(.amatrix_host_arg(value))
}

#' Block Lanczos SVD via block Krylov iteration
#'
#' Computes a truncated SVD using a block Lanczos bidiagonalization. Each
#' Krylov step issues one GPU GEMM per block rather than sequential GEMVs,
#' significantly reducing kernel-launch overhead on accelerated backends.
#'
#' @param A Numeric matrix, \code{adgeMatrix}, or \code{adgCMatrix}. Plain
#'   matrices are coerced to \code{adgeMatrix} using \code{mode} and
#'   \code{backend}.
#' @param nv Number of right singular vectors to return.
#' @param nu Number of left singular vectors to return. Defaults to \code{nv}.
#' @param block_size Integer block width for the Krylov iteration. When
#'   \code{NULL} (default), a size is chosen automatically based on \code{nv}
#'   and \code{nu}.
#' @param n_steps Number of Krylov steps. When \code{NULL} (default), chosen
#'   automatically.
#' @param mode Execution mode passed to \code{adgeMatrix()} when coercing
#'   plain matrices.
#' @param backend Backend name passed to \code{adgeMatrix()} when coercing
#'   plain matrices. Ignored when \code{A} is already an \code{aMatrix}.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{u}{Numeric matrix \code{[m, nu]}: left singular vectors.}
#'     \item{d}{Numeric vector of length \code{min(nu, nv)}: singular values
#'       in decreasing order.}
#'     \item{v}{Numeric matrix \code{[n, nv]}: right singular vectors.}
#'     \item{iter}{Integer: number of Krylov steps performed.}
#'     \item{mprod}{Integer: total matrix-vector products issued.}
#'   }
#'
#' @examples
#' A <- matrix(rnorm(200), nrow = 20)
#' res <- block_lanczos(A, nv = 3L)
#' length(res$d)
#'
#' @seealso \code{\link{rsvd}}, \code{\link{block_svd}}
#' @export
block_lanczos <- function(A,
                             nv = 5L,
                             nu = nv,
                             block_size = NULL,
                             n_steps = NULL,
                             mode = "fast",
                             backend = NULL) {
  if (!inherits(A, c("adgeMatrix", "adgCMatrix"))) {
    A <- adgeMatrix(A, mode = mode, backend = backend)
  }

  nv <- as.integer(nv)
  nu <- as.integer(nu)
  if (is.na(nv) || nv < 1L) {
    stop("nv must be a positive integer", call. = FALSE)
  }
  if (is.na(nu) || nu < 1L) {
    stop("nu must be a positive integer", call. = FALSE)
  }
  k <- max(nv, nu)

  if (is.null(block_size)) {
    block_size <- .amatrix_block_lanczos_default_block_size(k)
  }
  block_size <- as.integer(block_size)
  if (is.na(block_size) || block_size < 1L) {
    stop("block_size must be a positive integer", call. = FALSE)
  }

  b <- min(block_size, NROW(A), NCOL(A))
  if (b < 1L) {
    stop("block_size must not exceed matrix dimensions", call. = FALSE)
  }
  J <- if (is.null(n_steps)) .amatrix_block_lanczos_default_steps(k, b) else as.integer(n_steps)
  if (is.na(J) || J < 1L) {
    stop("n_steps must be a positive integer", call. = FALSE)
  }

  # Block Krylov iteration: GEMM-per-block instead of sequential GEMVs
  m <- NROW(A)
  n <- NCOL(A)
  A_left <- .amatrix_block_lanczos_source_operator(A)
  # Register the source-operator drop immediately after its (possibly resident)
  # upload, BEFORE building the right operator. If the right-operator build
  # throws, this on.exit still fires and releases A_left's device buffer.
  on.exit(.amatrix_block_lanczos_drop_source_operator(A_left), add = TRUE)
  A_right <- .amatrix_block_lanczos_right_operator(A, source_operator = A_left)
  on.exit(.amatrix_block_lanczos_drop_right_operator(A_right), add = TRUE)

  Q_left_basis <- matrix(0, nrow = m, ncol = J * b)
  Q_right_basis <- matrix(0, nrow = n, ncol = J * b)
  B_proj <- matrix(0, nrow = J * b, ncol = J * b)

  # Starting right block: random n×b, orthonormalized
  Q_cur <- .amatrix_block_thin_qr(matrix(rnorm(n * b), n, b))$q
  storage.mode(Q_cur) <- "double"

  for (j in seq_len(J)) {
    prev_cols <- (j - 1L) * b
    cur_cols <- ((j - 1L) * b + 1L):(j * b)

    # Z = A %*% Q_cur  — GPU GEMM (b columns, not b sequential GEMVs)
    Z_left_raw <- .amatrix_block_lanczos_left_product(A, Q_cur, operator = A_left)  # m × b
    left_reorth <- .amatrix_block_reorth_prefix(
      Z_left_raw,
      Q_left_basis,
      prev_cols,
      return_projection = TRUE
    )
    Z_left <- left_reorth$z
    left_factor <- .amatrix_block_thin_qr(Z_left)
    QL_j    <- left_factor$q               # m × b, CPU thin factor (b << m)
    storage.mode(QL_j) <- "double"
    R_left_j <- left_factor$r
    storage.mode(R_left_j) <- "double"
    Q_left_basis[, cur_cols] <- QL_j
    left_coeff <- left_reorth$coeff

    # A %*% Q_{R,j-1} lives entirely in span(Q_{L,1:j}) after reorthogonalization,
    # so we can assemble that projected block from the same coefficients instead of
    # issuing a second full projection pass later.
    if (j > 1L) {
      col_idx <- (prev_cols - b + 1L):prev_cols
      if (!is.null(left_coeff) && nrow(left_coeff) > 0L) {
        B_proj[seq_len(nrow(left_coeff)), col_idx] <- left_coeff
      }
      row_idx <- cur_cols
      B_proj[row_idx, col_idx] <- R_left_j[seq_len(b), seq_len(b), drop = FALSE]
    }

    # W = t(A) %*% Q_j  — GPU GEMM
    Z_right_raw <- .amatrix_block_lanczos_right_product(A, QL_j, operator = A_right)
    Z_right <- .amatrix_block_reorth_prefix(Z_right_raw, Q_right_basis, prev_cols)
    QR_j    <- .amatrix_block_thin_qr(Z_right)$q
    storage.mode(QR_j) <- "double"
    Q_right_basis[, cur_cols] <- QR_j

    Q_cur <- QR_j
  }

  # Collect and re-orthogonalize bases (CPU; J*b << m,n for typical settings)
  total_cols <- J * b
  Q_L <- Q_left_basis[, seq_len(total_cols), drop = FALSE]    # m × (J*b)
  Q_R <- Q_right_basis[, seq_len(total_cols), drop = FALSE]   # n × (J*b)
  needs_final_left_qr <- J > 1L && .amatrix_block_basis_needs_final_qr(Q_L)
  if (needs_final_left_qr) {
    Q_L <- qr.Q(qr(Q_L))
    storage.mode(Q_L) <- "double"
  }
  needs_final_right_qr <- J > 1L && .amatrix_block_basis_needs_final_qr(Q_R)
  if (needs_final_right_qr) {
    Q_R <- qr.Q(qr(Q_R))
    storage.mode(Q_R) <- "double"
  }

  if (needs_final_left_qr || needs_final_right_qr) {
    # If a final global QR changes the basis, fall back to the explicit
    # projection so the small matrix still matches the re-orthogonalized basis.
    AQ_R <- as.matrix(A %*% Q_R)           # m × (J*b), GPU GEMM
    B    <- base::crossprod(Q_L, AQ_R)           # (J*b) × (J*b), CPU
  } else {
    last_col_idx <- ((J - 1L) * b + 1L):(J * b)
    AQ_last <- .amatrix_block_lanczos_left_product(A, Q_cur, operator = A_left)  # A %*% Q_{R,J}
    B <- B_proj
    B[, last_col_idx] <- base::crossprod(Q_L, AQ_last)
  }

  k_out <- min(k, ncol(Q_L), ncol(Q_R), nrow(B), ncol(B))
  svd_B <- base::svd(B, nu = k_out, nv = k_out)

  # Lift singular vectors back to original space
  U <- Q_L %*% svd_B$u[, seq_len(k_out), drop = FALSE]  # m × k_out
  V <- Q_R %*% svd_B$v[, seq_len(k_out), drop = FALSE]  # n × k_out

  list(
    u     = U[, seq_len(min(nu, k_out)), drop = FALSE],
    d     = svd_B$d[seq_len(k_out)],
    v     = V[, seq_len(min(nv, k_out)), drop = FALSE],
    iter  = J,
    mprod = 2L * J * as.integer(b)
  )
}

#' @rdname block_lanczos
#' @param k Number of singular values/vectors. Alias for \code{nv = nu = k}.
#' @export
block_svd <- function(A,
                         k,
                         block_size = NULL,
                         n_steps = NULL,
                         mode = "fast",
                         backend = NULL) {
  block_lanczos(
    A,
    nv = as.integer(k),
    nu = as.integer(k),
    block_size = block_size,
    n_steps = n_steps,
    mode = mode,
    backend = backend
  )
}
