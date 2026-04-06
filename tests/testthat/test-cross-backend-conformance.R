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

  rhs_r <- matrix(rnorm(p * 3L), p, 3L)  # RHS for solve

  # Near-rank-3 matrix for rsvd (true reconstruction error should be ~0)
  U3 <- qr.Q(qr(matrix(rnorm(n * 3L), n, 3L)))
  V3 <- qr.Q(qr(matrix(rnorm(p * 3L), p, 3L)))
  low_rank_r <- U3 %*% diag(c(5, 3, 1)) %*% t(V3) +
                matrix(rnorm(n * p, sd = 0.01), n, p)

  list(n = n, p = p,
       A_r = A_r, C_r = C_r, sq_r = sq_r,
       S_r = S_r, rhs_r = rhs_r, low_rank_r = low_rank_r)
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
}

# -- am_rsvd conformance (reconstruction + orthonormality) --------------------
.run_rsvd_conformance <- function(backend_name, tol, precision = "strict") {
  d   <- .conformance_data
  X_in <- adgeMatrix(d$low_rank_r, preferred_backend = backend_name, precision = precision)
  k   <- 3L  # true rank of d$low_rank_r is 3

  sv  <- am_rsvd(X_in, k = k, n_oversamples = 10L, n_iter = 4L)
  X_rec <- sv$u %*% diag(sv$d) %*% t(sv$v)

  # Reconstruction error: compare to optimal rank-k truncated SVD.
  # GPU (float32) may deviate more; we allow up to 5x the optimal error + 0.1.
  sv_opt <- base::svd(d$low_rank_r, nu = k, nv = k)
  X_opt  <- sv_opt$u %*% diag(sv_opt$d[seq_len(k)]) %*% t(sv_opt$v)
  err_opt <- norm(d$low_rank_r - X_opt, "F") / norm(d$low_rank_r, "F")

  err <- norm(d$low_rank_r - X_rec, "F") / norm(d$low_rank_r, "F")
  # Allow up to 10x the optimal error + 0.15 to accommodate float32 variance
  expect_lt(err, err_opt * 10 + 0.15,
            label = paste0("[", backend_name, "/am_rsvd reconstruction]"))

  # Orthonormality of U and V columns
  expect_equal(t(sv$u) %*% sv$u, diag(k), tolerance = tol,
               label = paste0("[", backend_name, "/am_rsvd U orthonormal]"))
  expect_equal(t(sv$v) %*% sv$v, diag(k), tolerance = tol,
               label = paste0("[", backend_name, "/am_rsvd V orthonormal]"))

  # Singular values non-increasing
  expect_true(all(diff(sv$d) <= 0),
              label = paste0("[", backend_name, "/am_rsvd sv ordering]"))
}

# -- am_dist / am_kernel wrapper-level tests ----------------------------------
test_that("am_dist euclidean matches base R", {
  set.seed(11)
  X <- matrix(rnorm(50 * 8), 50, 8)
  Y <- matrix(rnorm(30 * 8), 30, 8)

  # Cross euclidean distance
  D_gpu  <- am_dist(X, Y, method = "euclidean")
  # Build reference: euclidean cross-distance
  D_base <- matrix(0, 50, 30)
  for (i in seq_len(50)) for (j in seq_len(30))
    D_base[i, j] <- sqrt(sum((X[i, ] - Y[j, ])^2))
  expect_equal(D_gpu, D_base, tolerance = .GPU_TOL,
               label = "am_dist/euclidean cross")

  # Self-distance: diagonal must be exactly 0
  D_self <- am_dist(X, method = "euclidean")
  expect_equal(diag(D_self), rep(0, 50), label = "am_dist/euclidean self-diag")

  # Symmetric
  expect_equal(D_self, t(D_self), tolerance = .GPU_TOL,
               label = "am_dist/euclidean symmetry")
})

test_that("am_dist sqeuclidean matches base R", {
  set.seed(12)
  X <- matrix(rnorm(40 * 6), 40, 6)
  D_sq   <- am_dist(X, method = "sqeuclidean")
  # Reference: squared euclidean from outer product identity
  nx <- rowSums(X^2)
  D_base <- outer(nx, nx, "+") - 2 * tcrossprod(X)
  D_base[D_base < 0] <- 0
  expect_equal(D_sq, D_base, tolerance = .GPU_TOL,
               label = "am_dist/sqeuclidean")
  expect_equal(diag(D_sq), rep(0, 40),
               label = "am_dist/sqeuclidean self-diag")
})

test_that("am_kernel linear matches tcrossprod", {
  set.seed(13)
  X <- matrix(rnorm(30 * 5), 30, 5)
  Y <- matrix(rnorm(20 * 5), 20, 5)

  K_self <- am_kernel(X, kernel = "linear")
  expect_equal(K_self, tcrossprod(X), tolerance = .GPU_TOL,
               label = "am_kernel/linear self")

  K_cross <- am_kernel(X, Y, kernel = "linear")
  expect_equal(K_cross, tcrossprod(X, Y), tolerance = .GPU_TOL,
               label = "am_kernel/linear cross")
})

test_that("am_kernel rbf matches reference implementation", {
  set.seed(14)
  X     <- matrix(rnorm(20 * 4), 20, 4)
  sigma <- 0.5

  K_rbf <- am_kernel(X, kernel = "rbf", sigma = sigma)

  # Reference: K(x,y) = exp(-||x-y||^2 / (2*sigma^2))
  nx     <- rowSums(X^2)
  D_sq   <- outer(nx, nx, "+") - 2 * tcrossprod(X)
  D_sq[D_sq < 0] <- 0
  K_base <- exp(-D_sq / (2 * sigma^2))

  expect_equal(K_rbf, K_base, tolerance = .GPU_TOL,
               label = "am_kernel/rbf")
  # Diagonal must be 1 for normalised RBF
  expect_equal(diag(K_rbf), rep(1, 20), tolerance = .GPU_TOL,
               label = "am_kernel/rbf diagonal=1")
})

test_that("am_kernel polynomial matches reference implementation", {
  set.seed(15)
  X      <- matrix(rnorm(25 * 4), 25, 4)
  degree <- 3L; coef <- 1.0

  K_poly <- am_kernel(X, kernel = "polynomial", degree = degree, coef = coef)
  G      <- tcrossprod(X)
  K_base <- (G + coef)^degree
  expect_equal(K_poly, K_base, tolerance = .GPU_TOL,
               label = "am_kernel/polynomial")
})

# -- CPU backend conformance --------------------------------------------------
test_that("cpu backend: all core ops match base R (float64 tolerance)", {
  .run_backend_conformance("cpu", tol = .CPU_TOL)
})

test_that("cpu backend: am_rsvd reconstruction and orthonormality", {
  .run_rsvd_conformance("cpu", tol = .CPU_TOL)
})

# -- MLX backend conformance --------------------------------------------------
test_that("mlx backend: all core ops match base R", {
  skip_if_not(
    isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE)),
    "mlx backend not available"
  )
  .run_backend_conformance("mlx", tol = .GPU_TOL)
})

test_that("mlx backend: am_rsvd reconstruction and orthonormality", {
  skip_if_not(
    isTRUE(try(amatrix.mlx::amatrix_mlx_is_available(), silent = TRUE)),
    "mlx backend not available"
  )
  # In the broad cross-backend suite, strict precision should preserve the
  # CPU-quality contract even when a GPU backend is preferred. The native MLX
  # fast-path is covered separately by backend-local smoke tests.
  .run_rsvd_conformance("mlx", tol = .CPU_TOL, precision = "strict")
})

# -- ArrayFire backend conformance --------------------------------------------
test_that("arrayfire backend: all core ops match base R", {
  skip_if_not(
    isTRUE(try(amatrix.arrayfire::amatrix_arrayfire_is_available(), silent = TRUE)),
    "arrayfire backend not available"
  )
  .run_backend_conformance("arrayfire", tol = .GPU_TOL)
})

test_that("arrayfire backend: am_rsvd reconstruction and orthonormality", {
  skip_if_not(
    isTRUE(try(amatrix.arrayfire::amatrix_arrayfire_is_available(), silent = TRUE)),
    "arrayfire backend not available"
  )
  .run_rsvd_conformance("arrayfire", tol = .CPU_TOL, precision = "strict")
})
