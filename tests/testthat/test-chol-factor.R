make_spd <- function(n, seed = 1L) {
  set.seed(seed)
  A <- matrix(rnorm(n * n), n, n)
  crossprod(A) + diag(n)
}

make_ridge_spd <- function(n_obs, p, lambda = 0.5, seed = 1L) {
  set.seed(seed)
  X <- matrix(rnorm(n_obs * p), nrow = n_obs, ncol = p)
  crossprod(X) + diag(lambda, p)
}

make_kernel_spd <- function(n_obs, p, sigma = 1.0, jitter = 0.1, seed = 1L) {
  set.seed(seed)
  X <- matrix(rnorm(n_obs * p), nrow = n_obs, ncol = p)
  kernel_matrix(X, kernel = "rbf", sigma = sigma) + diag(jitter, n_obs)
}

.mlx_test_spec <- function() {
  specs <- optional_backend_specs()
  specs[[match("mlx", vapply(specs, `[[`, character(1), "backend"))]]
}

.frob_norm <- function(x) {
  sqrt(sum(x * x))
}

test_that("chol_factor returns correct upper triangular factor", {
  M <- make_spd(8L, seed = 42L)
  X <- as_adgeMatrix(M)

  fac <- chol_factor(X)
  expect_s4_class(fac, "amChol")

  R <- as.matrix(fac)
  R_ref <- chol(as.matrix(M))

  expect_equal(R, R_ref, tolerance = 1e-10)
  expect_equal(t(R) %*% R, as.matrix(M), tolerance = 1e-10)
  # upper triangular
  expect_true(all(abs(R[lower.tri(R)]) < 1e-12))
})

test_that("chol_solve matches base::solve for single RHS vector", {
  M <- make_spd(10L, seed = 7L)
  X <- as_adgeMatrix(M)
  fac <- chol_factor(X)

  b <- rnorm(10L)
  sol <- chol_solve(fac, b)
  ref <- base::solve(M, b)

  expect_equal(as.numeric(sol), as.numeric(ref), tolerance = 1e-10)
})

test_that("chol_solve matches base::solve for multi-column B (k=1, 10, 100)", {
  n <- 12L
  M <- make_spd(n, seed = 11L)
  X <- as_adgeMatrix(M)
  fac <- chol_factor(X)

  for (k in c(1L, 10L, 100L)) {
    set.seed(k)
    B <- matrix(rnorm(n * k), n, k)
    sol <- chol_solve(fac, B)
    ref <- base::solve(M, B)
    expect_equal(dim(sol), c(n, k))
    expect_lt(max(abs(sol - ref)), 1e-10)
  }
})

test_that("chol_factor reuses cache on second call", {
  M <- make_spd(6L, seed = 3L)
  X <- as_adgeMatrix(M)

  fac1 <- chol_factor(X)
  fac2 <- chol_factor(X)

  expect_identical(fac1@factor, fac2@factor)
  expect_identical(fac1, fac2)
})

test_that("amChol show method runs without error", {
  M <- make_spd(4L, seed = 5L)
  X <- as_adgeMatrix(M)
  fac <- chol_factor(X)

  expect_output(show(fac), "amChol")
})

test_that("fast MLX Cholesky path matches CPU on ridge-like SPD many-RHS solves", {
  spec <- .mlx_test_spec()
  skip_if_backend_package_missing(spec)

  mlx_ns <- optional_backend_namespace(spec$package)
  native_available <- get("amatrix_mlx_native_available", envir = mlx_ns, inherits = FALSE)
  skip_if_not(isTRUE(native_available()), "mlx native backend not available")

  with_optional_backend_available(spec, {
    A <- make_ridge_spd(n_obs = 320L, p = 96L, lambda = 0.75, seed = 20260406L)
    B <- matrix(rnorm(96L * 24L), nrow = 96L, ncol = 24L)
    X_mlx <- adgeMatrix(A, preferred_backend = "mlx", precision = "fast")

    expect_identical(amatrix_backend_plan(X_mlx, "chol")$chosen, "mlx")

    fac <- chol_factor(X_mlx)
    sol <- chol_solve(fac, B)
    ref_sol <- solve(A, B)
    by_col <- vapply(
      seq_len(ncol(B)),
      function(j) as.numeric(chol_solve(fac, B[, j])),
      numeric(nrow(B))
    )
    recon_rel <- .frob_norm(crossprod(fac@factor) - A) / .frob_norm(A)
    solve_ref_rel <- .frob_norm(sol - ref_sol) / .frob_norm(ref_sol)
    solve_resid_rel <- .frob_norm(A %*% sol - B) / .frob_norm(B)
    batched_rel <- .frob_norm(sol - by_col) / .frob_norm(by_col)

    expect_s4_class(fac, "amChol")
    expect_identical(fac@precision, "fast")
    expect_identical(fac@backend, "mlx")
    expect_true(all(is.finite(sol)))
    expect_lt(recon_rel, 5e-6)
    expect_lt(solve_ref_rel, 5e-6)
    expect_lt(solve_resid_rel, 5e-6)
    expect_lt(batched_rel, 5e-6)
  })
})

test_that("fast MLX Cholesky path stays aligned on kernel-like SPD systems", {
  spec <- .mlx_test_spec()
  skip_if_backend_package_missing(spec)

  mlx_ns <- optional_backend_namespace(spec$package)
  native_available <- get("amatrix_mlx_native_available", envir = mlx_ns, inherits = FALSE)
  skip_if_not(isTRUE(native_available()), "mlx native backend not available")

  with_optional_backend_available(spec, {
    A <- make_kernel_spd(n_obs = 96L, p = 8L, sigma = 1.1, jitter = 0.2, seed = 20260407L)
    B <- matrix(rnorm(96L * 12L), nrow = 96L, ncol = 12L)
    X_mlx <- adgeMatrix(A, preferred_backend = "mlx", precision = "fast")

    expect_identical(amatrix_backend_plan(X_mlx, "chol")$chosen, "mlx")

    fac <- chol_factor(X_mlx)
    sol <- chol_solve(fac, B)
    ref_sol <- solve(A, B)
    recon_rel <- .frob_norm(crossprod(fac@factor) - A) / .frob_norm(A)
    solve_ref_rel <- .frob_norm(sol - ref_sol) / .frob_norm(ref_sol)
    solve_resid_rel <- .frob_norm(A %*% sol - B) / .frob_norm(B)

    expect_s4_class(fac, "amChol")
    expect_true(all(is.finite(fac@factor)))
    expect_true(all(is.finite(sol)))
    expect_equal(dim(sol), dim(B))
    expect_lt(recon_rel, 5e-6)
    expect_lt(solve_ref_rel, 5e-6)
    expect_lt(solve_resid_rel, 5e-6)
  })
})
