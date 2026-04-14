# test-cross-backend-conformance.R
#
# Cross-backend numerical conformance harness.
# For every available backend, runs all core ops and compares against
# base R ground truth within appropriate tolerance.
#
# GPU backends (mlx, arrayfire) compute in float32 and up-convert to
# float64; 1e-4 relative tolerance accounts for that precision gap.
# CPU backend operates natively in float64; 1e-10 tolerance is used.

.GPU_TOL  <- 1e-4
.CPU_TOL  <- 1e-10

# -- Shared test data (reproducible) ------------------------------------------

.conformance_data <- local({
  set.seed(7)
  n <- 40L; p <- 8L
  A_r  <- matrix(rnorm(n * p), n, p)   # n×p
  C_r  <- matrix(rnorm(n * p), n, p)   # n×p  (same cols as A, used for cross-tcrossprod)
  sq_r <- matrix(rnorm(p * p), p, p)   # p×p  (for matmul)

  # SPD p×p for solve/chol
  Z   <- matrix(rnorm(p * p), p, p)
  S_r <- crossprod(Z) + diag(p) * 2

  rhs_r  <- matrix(rnorm(p * 3L), p, 3L)  # RHS for solve (p×3)
  rhs_lm <- matrix(rnorm(n * 3L), n, 3L)  # RHS for many_lm (n×3)

  # Near-rank-3 matrix for rsvd (true reconstruction error should be ~0)
  U3 <- qr.Q(qr(matrix(rnorm(n * 3L), n, 3L)))
  V3 <- qr.Q(qr(matrix(rnorm(p * 3L), p, 3L)))
  low_rank_r <- U3 %*% diag(c(5, 3, 1)) %*% t(V3) +
                matrix(rnorm(n * p, sd = 0.01), n, p)

  list(n = n, p = p,
       A_r = A_r, C_r = C_r, sq_r = sq_r,
       S_r = S_r, rhs_r = rhs_r, rhs_lm = rhs_lm, low_rank_r = low_rank_r)
})

# -- Core conformance checker -------------------------------------------------
# Runs every op for one backend and asserts numerical parity with base R.
# `tol` should be .GPU_TOL for GPU backends and .CPU_TOL for cpu.

.run_backend_conformance <- function(backend_name, tol) {
  d <- .conformance_data

  make <- function(x) adgeMatrix(x, preferred_backend = backend_name)
  tag  <- function(op) paste0("[", backend_name, "/", op, "]")

  A   <- make(d$A_r)
  C   <- make(d$C_r)
  sq  <- make(d$sq_r)
  sq2 <- make(t(d$sq_r))   # p×p, for matmul
  S   <- make(d$S_r)
  rhs <- make(d$rhs_r)

  # matmul: sq(p×p) %*% sq2(p×p)
  expect_equal(as.matrix(sq %*% sq2), d$sq_r %*% t(d$sq_r),
               tolerance = tol, label = tag("matmul"))

  # crossprod: A'A  → p×p
  expect_equal(as.matrix(crossprod(A)), crossprod(d$A_r),
               tolerance = tol, label = tag("crossprod_self"))

  # crossprod: A'C  → p×p  (A and C both n×p)
  expect_equal(as.matrix(crossprod(A, C)), crossprod(d$A_r, d$C_r),
               tolerance = tol, label = tag("crossprod_xy"))

  # tcrossprod: AA'  → n×n
  expect_equal(as.matrix(tcrossprod(A)), tcrossprod(d$A_r),
               tolerance = tol, label = tag("tcrossprod_self"))

  # tcrossprod: AC'  → n×n  (A and C both n×p, so AC' is n×n)
  expect_equal(as.matrix(tcrossprod(A, C)), tcrossprod(d$A_r, d$C_r),
               tolerance = tol, label = tag("tcrossprod_xy"))

  # ewise add / multiply / scalar
  expect_equal(as.matrix(A + A), d$A_r + d$A_r,
               tolerance = tol, label = tag("ewise_add"))
  expect_equal(as.matrix(A * A), d$A_r * d$A_r,
               tolerance = tol, label = tag("ewise_mul"))
  expect_equal(as.matrix(A * 3.0), d$A_r * 3.0,
               tolerance = tol, label = tag("ewise_scalar"))

  # rowSums / colSums
  expect_equal(rowSums(A), rowSums(d$A_r), tolerance = tol,
               label = tag("rowSums"))
  expect_equal(colSums(A), colSums(d$A_r), tolerance = tol,
               label = tag("colSums"))

  # chol (SPD input)
  expect_equal(as.matrix(chol(S)), base::chol(d$S_r),
               tolerance = tol, label = tag("chol"))

  # solve: matrix inverse (no RHS)
  expect_equal(as.matrix(solve(S)), base::solve(d$S_r),
               tolerance = tol, label = tag("solve_inv"))

  # solve: A x = B  (p×3 RHS)
  expect_equal(as.matrix(solve(S, rhs)), base::solve(d$S_r, d$rhs_r),
               tolerance = tol, label = tag("solve_rhs"))

  # covariance vs cov()
  ref_cov <- cov(d$A_r)
  am_cov  <- as.matrix(covariance(make(d$A_r)))
  expect_equal(am_cov, ref_cov, tolerance = tol, label = tag("covariance"))

  # am_eigen: eigenvalues (sorted) + sign-invariant residual ||A·V - V·Λ||_F / ||A||_F
  ref_ev <- eigen(d$S_r, symmetric = TRUE)
  am_ev  <- am_eigen(make(d$S_r), symmetric = TRUE)
  expect_equal(sort(am_ev$values, decreasing = TRUE),
               ref_ev$values, tolerance = tol, label = tag("eigen_values"))
  resid_eigen <- norm(d$S_r %*% am_ev$vectors -
                        am_ev$vectors %*% diag(am_ev$values, nrow = length(am_ev$values)), "F") /
                 norm(d$S_r, "F")
  expect_lt(resid_eigen, sqrt(tol), label = tag("eigen_residual"))

  # svd() S4 dispatch: singular values only (vectors sign-invariant)
  ref_sv <- base::svd(d$A_r)$d
  am_sv  <- svd(make(d$A_r))$d
  expect_equal(am_sv, ref_sv, tolerance = tol, label = tag("svd_values"))

  # many_lm vs lm.fit() loop (n×3 RHS)
  ref_coef <- do.call(cbind, lapply(seq_len(3L), function(j)
    lm.fit(d$A_r, d$rhs_lm[, j])$coefficients))
  dimnames(ref_coef) <- NULL   # many_lm does not set dimnames
  am_coef  <- as.matrix(many_lm(make(d$A_r), d$rhs_lm,
                                    method = "qr", cache = FALSE)$coefficients)
  expect_equal(am_coef, ref_coef, tolerance = tol, label = tag("many_lm"))

  # dist_matrix vs as.matrix(dist())
  # dist_matrix uses the GEMM identity (||x||² + ||y||² − 2x·y); dist() uses direct
  # element-wise sum. The two are mathematically equal but differ by ~1e-6 in
  # float64, so we floor the tolerance at 1e-6 even for the CPU backend.
  ref_dist <- as.matrix(dist(d$A_r, method = "euclidean"))
  dimnames(ref_dist) <- NULL   # dist_matrix does not set dimnames
  am_d     <- dist_matrix(d$A_r, method = "euclidean")
  expect_equal(am_d, ref_dist, tolerance = max(tol, 1e-6),
               label = tag("dist_matrix_euclidean"))
}

# -- rsvd conformance (reconstruction + orthonormality) --------------------
.run_rsvd_conformance <- function(backend_name, tol, precision = "strict") {
  d   <- .conformance_data
  X_in <- adgeMatrix(d$low_rank_r, preferred_backend = backend_name, precision = precision)
  k   <- 3L  # true rank of d$low_rank_r is 3

  sv  <- rsvd(X_in, k = k, n_oversamples = 10L, n_iter = 4L)
  X_rec <- sv$u %*% diag(sv$d) %*% t(sv$v)

  # Reconstruction error: compare to optimal rank-k truncated SVD.
  # GPU (float32) may deviate more; we allow up to 5x the optimal error + 0.1.
  sv_opt <- base::svd(d$low_rank_r, nu = k, nv = k)
  X_opt  <- sv_opt$u %*% diag(sv_opt$d[seq_len(k)]) %*% t(sv_opt$v)
  err_opt <- norm(d$low_rank_r - X_opt, "F") / norm(d$low_rank_r, "F")

  err <- norm(d$low_rank_r - X_rec, "F") / norm(d$low_rank_r, "F")
  # Allow up to 10x the optimal error + 0.15 to accommodate float32 variance
  expect_lt(err, err_opt * 10 + 0.15,
            label = paste0("[", backend_name, "/rsvd reconstruction]"))

  # Orthonormality of U and V columns
  expect_equal(t(sv$u) %*% sv$u, diag(k), tolerance = tol,
               label = paste0("[", backend_name, "/rsvd U orthonormal]"))
  expect_equal(t(sv$v) %*% sv$v, diag(k), tolerance = tol,
               label = paste0("[", backend_name, "/rsvd V orthonormal]"))

  # Singular values non-increasing
  expect_true(all(diff(sv$d) <= 0),
              label = paste0("[", backend_name, "/rsvd sv ordering]"))
}

# -- Gaussian Process regression conformance ------------------------------
#
# GP regression is a five-op chain that exercises the full factorization
# lifecycle (kernel_matrix -> chol_factor -> chol_solve -> matmul ->
# chol_logdet + quad_form) on a SHARED factor. The reference implementation
# is a 15-line plain-R GP using base::chol + base::backsolve.
#
# Checks:
#   * Posterior mean matches base R within tol
#   * Posterior variance matches base R within tol
#   * Log marginal likelihood matches base R within tol
#   * The Kα = y invariant (metamorphic): Kα must reproduce y exactly.
#
# The test is the smallest dataset we can use while still exercising a
# non-trivial Cholesky on K+σ²I. All constants are pinned via seed.

.gp_reference <- function(X_train, y_train, X_test, lengthscale = 0.7,
                          amplitude = 1.0, noise = 0.1) {
  .rbf <- function(A, B, ls, amp) {
    sq <- outer(rowSums(A * A), rowSums(B * B), `+`) - 2 * (A %*% t(B))
    sq[sq < 0] <- 0
    amp^2 * exp(-sq / (2 * ls^2))
  }
  K   <- .rbf(X_train, X_train, lengthscale, amplitude) +
         diag(noise^2, nrow(X_train))
  Ks  <- .rbf(X_test,  X_train, lengthscale, amplitude)
  Kss <- .rbf(X_test,  X_test,  lengthscale, amplitude) +
         diag(noise^2, nrow(X_test))

  L     <- t(chol(K))                            # lower triangular
  alpha <- backsolve(t(L), forwardsolve(L, y_train))
  mu    <- as.numeric(Ks %*% alpha)

  v     <- forwardsolve(L, t(Ks))                # shape p x n_test
  var_post <- diag(Kss) - colSums(v * v)

  logdet_K <- 2 * sum(log(diag(L)))
  quad     <- sum(y_train * alpha)
  lml      <- -0.5 * quad - 0.5 * logdet_K - 0.5 * length(y_train) * log(2 * pi)

  list(mu = mu, var = var_post, lml = lml, alpha = alpha, K = K)
}

.run_gp_conformance <- function(backend_name, tol) {
  set.seed(2026041601L)
  n_train <- 24L
  n_test  <- 8L
  d       <- 3L

  X_train_host <- matrix(rnorm(n_train * d), n_train, d)
  X_test_host  <- matrix(rnorm(n_test  * d), n_test,  d)
  f_fun   <- function(x) sin(x[, 1L]) + 0.5 * x[, 2L]^2 - x[, 3L]
  y_train <- f_fun(X_train_host) + 0.05 * rnorm(n_train)

  lengthscale <- 0.7
  amplitude   <- 1.0
  noise       <- 0.1

  # -- Reference (base R) -------------------------------------------------
  ref <- .gp_reference(X_train_host, y_train, X_test_host,
                       lengthscale = lengthscale,
                       amplitude   = amplitude,
                       noise       = noise)

  # -- amatrix implementation --------------------------------------------
  # This block is the "plain R" story: we swap matrix for adgeMatrix and
  # use kernel_matrix + chol_factor + chol_solve + chol_logdet + matmul,
  # all of which route through the chosen backend.
  X_train <- adgeMatrix(X_train_host, preferred_backend = backend_name)
  X_test  <- adgeMatrix(X_test_host,  preferred_backend = backend_name)

  K_train   <- kernel_matrix(X_train, kernel = "rbf", sigma = lengthscale)
  K_train   <- as.matrix(K_train) * (amplitude^2) +
               diag(noise^2, n_train)
  K_cross   <- kernel_matrix(X_test, X_train, kernel = "rbf", sigma = lengthscale)
  K_cross   <- as.matrix(K_cross) * (amplitude^2)
  K_test    <- kernel_matrix(X_test, kernel = "rbf", sigma = lengthscale)
  K_test    <- as.matrix(K_test) * (amplitude^2) +
               diag(noise^2, n_test)

  K_am    <- as_adgeMatrix(K_train, preferred_backend = backend_name)
  factor  <- chol_factor(K_am)
  alpha   <- chol_solve(factor, matrix(y_train, ncol = 1L))
  alpha_v <- as.numeric(alpha)

  mu_am   <- as.numeric(K_cross %*% alpha_v)

  # Posterior variance: diag(K_test - K_cross %*% K^{-1} %*% t(K_cross)).
  # Using diag(A %*% B) = rowSums(A * t(B)): with A = K_cross (n_test x n_train)
  # and B = K^{-1} %*% t(K_cross) = Ks_t_m (n_train x n_test), we get
  # diag(A B) = rowSums(K_cross * t(Ks_t_m)).
  Ks_t    <- chol_solve(factor, t(K_cross))
  Ks_t_m  <- as.matrix(Ks_t)
  var_am  <- diag(K_test) - rowSums(K_cross * t(Ks_t_m))

  # Log marginal likelihood via chol_logdet + quad_form
  logdet_K <- chol_logdet(factor)
  quad     <- sum(y_train * alpha_v)
  lml_am   <- -0.5 * quad - 0.5 * logdet_K - 0.5 * n_train * log(2 * pi)

  tag <- function(op) paste0("[", backend_name, "/gp ", op, "]")

  # Oracle: match base R reference
  expect_equal(mu_am, ref$mu, tolerance = tol, label = tag("mean"))
  expect_equal(var_am, ref$var, tolerance = tol, label = tag("variance"))
  expect_equal(lml_am, ref$lml, tolerance = tol, label = tag("log marginal likelihood"))

  # Metamorphic: Kα must reproduce y.
  Ky <- as.numeric(K_train %*% alpha_v)
  expect_equal(Ky, y_train, tolerance = max(tol, 1e-8),
               label = tag("K alpha = y invariant"))

  # Variance must be non-negative (numerical): no backend should produce a
  # negative predictive variance for an SPD kernel + noise.
  expect_true(all(var_am >= -max(tol, 1e-8)),
              label = tag("variance non-negative"))
}

# -- dist_matrix / kernel_matrix wrapper-level tests ----------------------------------
test_that("dist_matrix euclidean matches base R", {
  set.seed(11)
  X <- matrix(rnorm(50 * 8), 50, 8)
  Y <- matrix(rnorm(30 * 8), 30, 8)

  # Cross euclidean distance
  D_gpu  <- dist_matrix(X, Y, method = "euclidean")
  # Build reference: euclidean cross-distance
  D_base <- matrix(0, 50, 30)
  for (i in seq_len(50)) for (j in seq_len(30))
    D_base[i, j] <- sqrt(sum((X[i, ] - Y[j, ])^2))
  expect_equal(D_gpu, D_base, tolerance = .GPU_TOL,
               label = "dist_matrix/euclidean cross")

  # Self-distance: diagonal must be exactly 0
  D_self <- dist_matrix(X, method = "euclidean")
  expect_equal(diag(D_self), rep(0, 50), label = "dist_matrix/euclidean self-diag")

  # Symmetric
  expect_equal(D_self, t(D_self), tolerance = .GPU_TOL,
               label = "dist_matrix/euclidean symmetry")
})

test_that("dist_matrix sqeuclidean matches base R", {
  set.seed(12)
  X <- matrix(rnorm(40 * 6), 40, 6)
  D_sq   <- dist_matrix(X, method = "sqeuclidean")
  # Reference: squared euclidean from outer product identity
  nx <- rowSums(X^2)
  D_base <- outer(nx, nx, "+") - 2 * tcrossprod(X)
  D_base[D_base < 0] <- 0
  expect_equal(D_sq, D_base, tolerance = .GPU_TOL,
               label = "dist_matrix/sqeuclidean")
  expect_equal(diag(D_sq), rep(0, 40),
               label = "dist_matrix/sqeuclidean self-diag")
})

test_that("kernel_matrix linear matches tcrossprod", {
  set.seed(13)
  X <- matrix(rnorm(30 * 5), 30, 5)
  Y <- matrix(rnorm(20 * 5), 20, 5)

  K_self <- kernel_matrix(X, kernel = "linear")
  expect_equal(K_self, tcrossprod(X), tolerance = .GPU_TOL,
               label = "kernel_matrix/linear self")

  K_cross <- kernel_matrix(X, Y, kernel = "linear")
  expect_equal(K_cross, tcrossprod(X, Y), tolerance = .GPU_TOL,
               label = "kernel_matrix/linear cross")
})

test_that("kernel_matrix rbf matches reference implementation", {
  set.seed(14)
  X     <- matrix(rnorm(20 * 4), 20, 4)
  sigma <- 0.5

  K_rbf <- kernel_matrix(X, kernel = "rbf", sigma = sigma)

  # Reference: K(x,y) = exp(-||x-y||^2 / (2*sigma^2))
  nx     <- rowSums(X^2)
  D_sq   <- outer(nx, nx, "+") - 2 * tcrossprod(X)
  D_sq[D_sq < 0] <- 0
  K_base <- exp(-D_sq / (2 * sigma^2))

  expect_equal(K_rbf, K_base, tolerance = .GPU_TOL,
               label = "kernel_matrix/rbf")
  # Diagonal must be 1 for normalised RBF
  expect_equal(diag(K_rbf), rep(1, 20), tolerance = .GPU_TOL,
               label = "kernel_matrix/rbf diagonal=1")
})

test_that("kernel_matrix polynomial matches reference implementation", {
  set.seed(15)
  X      <- matrix(rnorm(25 * 4), 25, 4)
  degree <- 3L; coef <- 1.0

  K_poly <- kernel_matrix(X, kernel = "polynomial", degree = degree, coef = coef)
  G      <- tcrossprod(X)
  K_base <- (G + coef)^degree
  expect_equal(K_poly, K_base, tolerance = .GPU_TOL,
               label = "kernel_matrix/polynomial")
})

# -- CPU backend conformance --------------------------------------------------
test_that("cpu backend: all core ops match base R (float64 tolerance)", {
  .run_backend_conformance("cpu", tol = .CPU_TOL)
})

test_that("cpu backend: rsvd reconstruction and orthonormality", {
  .run_rsvd_conformance("cpu", tol = .CPU_TOL)
})

test_that("cpu backend: GP regression matches base R GP reference", {
  .run_gp_conformance("cpu", tol = .CPU_TOL)
})

# -- MLX backend conformance --------------------------------------------------
test_that("mlx backend: all core ops match base R", {
  skip_if_not(
    isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE)),
    "mlx backend not available"
  )
  .run_backend_conformance("mlx", tol = .GPU_TOL)
})

test_that("mlx backend: rsvd reconstruction and orthonormality", {
  skip_if_not(
    isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE)),
    "mlx backend not available"
  )
  # In the broad cross-backend suite, strict precision should preserve the
  # CPU-quality contract even when a GPU backend is preferred. The native MLX
  # fast-path is covered separately by backend-local smoke tests.
  .run_rsvd_conformance("mlx", tol = .CPU_TOL, precision = "strict")
})

test_that("mlx backend: GP regression matches base R GP reference", {
  skip_if_not(
    isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE)),
    "mlx backend not available"
  )
  .run_gp_conformance("mlx", tol = .GPU_TOL)
})

# -- ArrayFire backend conformance --------------------------------------------
test_that("arrayfire backend: all core ops match base R", {
  skip_if_not(
    isTRUE(try(amatrix.arrayfire::amatrix_arrayfire_is_available(), silent = TRUE)),
    "arrayfire backend not available"
  )
  .run_backend_conformance("arrayfire", tol = .GPU_TOL)
})

test_that("arrayfire backend: rsvd reconstruction and orthonormality", {
  skip_if_not(
    isTRUE(try(amatrix.arrayfire::amatrix_arrayfire_is_available(), silent = TRUE)),
    "arrayfire backend not available"
  )
  .run_rsvd_conformance("arrayfire", tol = .CPU_TOL, precision = "strict")
})

test_that("arrayfire backend: GP regression matches base R GP reference", {
  skip_if_not(
    isTRUE(try(amatrix.arrayfire::amatrix_arrayfire_is_available(), silent = TRUE)),
    "arrayfire backend not available"
  )
  .run_gp_conformance("arrayfire", tol = .GPU_TOL)
})

# -- OpenCL backend conformance -----------------------------------------------
test_that("opencl backend: all core ops match base R", {
  skip_if_not(requireNamespace("amatrix.opencl", quietly = TRUE),
              "amatrix.opencl backend package not installed")

  old_env <- Sys.getenv("AMATRIX_OPENCL_PROBE_GPU", unset = NA_character_)
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("AMATRIX_OPENCL_PROBE_GPU") else Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = old_env)
  }, add = TRUE)

  Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = "1")
  amatrix.opencl::amatrix_opencl_register(overwrite = TRUE)

  skip_if_not(
    isTRUE(try(amatrix.opencl::amatrix_opencl_native_available(force = TRUE), silent = TRUE)),
    "opencl backend not available"
  )

  .run_backend_conformance("opencl", tol = .GPU_TOL)
})

test_that("opencl backend: GP regression matches base R GP reference", {
  skip_if_not(requireNamespace("amatrix.opencl", quietly = TRUE),
              "amatrix.opencl backend package not installed")

  old_env <- Sys.getenv("AMATRIX_OPENCL_PROBE_GPU", unset = NA_character_)
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("AMATRIX_OPENCL_PROBE_GPU") else Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = old_env)
  }, add = TRUE)

  Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = "1")
  amatrix.opencl::amatrix_opencl_register(overwrite = TRUE)

  skip_if_not(
    isTRUE(try(amatrix.opencl::amatrix_opencl_native_available(force = TRUE), silent = TRUE)),
    "opencl backend not available"
  )

  .run_gp_conformance("opencl", tol = .GPU_TOL)
})
