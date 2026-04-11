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

.opencl_test_spec <- function() {
  specs <- optional_backend_specs()
  specs[[match("opencl", vapply(specs, `[[`, character(1), "backend"))]]
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

test_that("chol_solve_batches matches repeated chol_solve calls", {
  n <- 10L
  M <- make_spd(n, seed = 13L)
  X <- as_adgeMatrix(M)
  fac <- chol_factor(X)

  rhs <- list(
    rnorm(n),
    matrix(rnorm(n * 3L), nrow = n, ncol = 3L),
    matrix(rnorm(n), nrow = n, ncol = 1L)
  )

  out <- chol_solve_batches(fac, rhs)
  ref <- lapply(rhs, function(b) chol_solve(fac, b))

  expect_length(out, length(rhs))
  for (idx in seq_along(rhs)) {
    expect_equal(out[[idx]], ref[[idx]], tolerance = 1e-10)
  }

  rhs_arr <- array(rnorm(n * 2L * 4L), dim = c(n, 2L, 4L))
  arr_out <- chol_solve_batches(fac, rhs_arr)
  expect_length(arr_out, 4L)
  for (idx in seq_along(arr_out)) {
    expect_equal(arr_out[[idx]], solve(M, rhs_arr[, , idx]), tolerance = 1e-10)
  }
})

test_that("chol_solve_batches validates RHS dimensions", {
  M <- make_spd(5L, seed = 17L)
  fac <- chol_factor(as_adgeMatrix(M))
  expect_error(chol_solve_batches(fac, list(matrix(rnorm(6L), nrow = 6L))), "nrow")
  expect_error(chol_solve_batches(fac, matrix(rnorm(5L), nrow = 5L)), "list of RHS")
})

test_that("chol_factor reuses cache on second call", {
  M <- make_spd(6L, seed = 3L)
  X <- as_adgeMatrix(M)

  fac1 <- chol_factor(X)
  fac2 <- chol_factor(X)

  expect_equal(as.matrix(fac1), as.matrix(fac2), tolerance = 1e-10)
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
    recon_rel <- .frob_norm(crossprod(as.matrix(fac)) - A) / .frob_norm(A)
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
    recon_rel <- .frob_norm(crossprod(as.matrix(fac)) - A) / .frob_norm(A)
    solve_ref_rel <- .frob_norm(sol - ref_sol) / .frob_norm(ref_sol)
    solve_resid_rel <- .frob_norm(A %*% sol - B) / .frob_norm(B)

    expect_s4_class(fac, "amChol")
    expect_true(all(is.finite(as.matrix(fac))))
    expect_true(all(is.finite(sol)))
    expect_equal(dim(sol), dim(B))
    expect_lt(recon_rel, 5e-6)
    expect_lt(solve_ref_rel, 5e-6)
    expect_lt(solve_resid_rel, 5e-6)
  })
})

test_that("fast OpenCL Cholesky path matches CPU on dense SPD solves", {
  spec <- .opencl_test_spec()
  skip_if_backend_package_missing(spec)

  ns <- optional_backend_namespace(spec$package)
  register_backend <- get(spec$register_fun, envir = ns, inherits = FALSE)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    A <- make_ridge_spd(n_obs = 256L, p = 64L, lambda = 0.5, seed = 20260408L)
    B <- matrix(rnorm(64L * 12L), nrow = 64L, ncol = 12L)
    X_opencl <- adgeMatrix(A, preferred_backend = "opencl", precision = "fast")

    expect_identical(amatrix_backend_plan(X_opencl, "chol")$chosen, "opencl")
    expect_identical(amatrix_backend_plan(X_opencl, "solve", y = adgeMatrix(B, preferred_backend = "opencl", precision = "fast"))$chosen, "opencl")

    fac <- chol_factor(X_opencl)
    sol <- chol_solve(fac, B)
    ref_sol <- solve(A, B)
    recon_rel <- .frob_norm(crossprod(as.matrix(fac)) - A) / .frob_norm(A)
    solve_ref_rel <- .frob_norm(sol - ref_sol) / .frob_norm(ref_sol)
    solve_resid_rel <- .frob_norm(A %*% sol - B) / .frob_norm(B)

    expect_s4_class(fac, "amChol")
    expect_identical(fac@backend, "opencl")
    expect_true(all(is.finite(as.matrix(fac))))
    expect_true(all(is.finite(sol)))
    expect_lt(recon_rel, 5e-6)
    expect_lt(solve_ref_rel, 5e-6)
    expect_lt(solve_resid_rel, 5e-6)
  })
})

test_that("fast OpenCL chol_factor retains resident factor for repeated solves", {
  spec <- .opencl_test_spec()
  skip_if_backend_package_missing(spec)

  ns <- optional_backend_namespace(spec$package)
  register_backend <- get(spec$register_fun, envir = ns, inherits = FALSE)
  native_available <- get("amatrix_opencl_native_available", envir = ns, inherits = FALSE)
  bridge_info <- get("amatrix_opencl_bridge_info", envir = ns, inherits = FALSE)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    skip_if_not(isTRUE(native_available(force = TRUE)))
    skip_if_not(isTRUE(bridge_info()$clblast))

    old <- options(
      amatrix.opencl.factor_gpu = TRUE,
      amatrix.opencl.factor_min_dim = 1L,
      amatrix.opencl.trsm_min_work = 1
    )
    on.exit(options(old), add = TRUE)

    A <- make_ridge_spd(n_obs = 192L, p = 48L, lambda = 0.75, seed = 20260409L)
    B <- matrix(rnorm(48L * 8L), nrow = 48L, ncol = 8L)
    X_opencl <- adgeMatrix(A, preferred_backend = "opencl", precision = "fast")
    B_opencl <- adgeMatrix(B, preferred_backend = "opencl", precision = "fast")

    fac <- chol_factor(X_opencl)
    tri <- solve_triangular(fac, B)
    sol <- chol_solve(fac, B)
    tri_opencl <- solve_triangular(fac, B_opencl)
    sol_opencl <- chol_solve(fac, B_opencl)
    batch_opencl <- chol_solve_batches(fac, list(B_opencl, B_opencl))

    expect_true(inherits(fac@factor_obj, "aMatrix"))
    expect_identical(amatrix:::.amatrix_live_resident_backend(fac@factor_obj), "opencl")
    expect_identical(length(fac@factor), 0L)
    expect_equal(tri, backsolve(chol(A), B), tolerance = 5e-6)
    expect_equal(sol, solve(A, B), tolerance = 5e-6)
    expect_equal(tri_opencl, backsolve(chol(A), B), tolerance = 5e-6)
    expect_s4_class(sol_opencl, "adgeMatrix")
    expect_identical(amatrix:::.amatrix_live_resident_backend(sol_opencl), "opencl")
    expect_equal(as.matrix(sol_opencl), solve(A, B), tolerance = 5e-6)
    expect_length(batch_opencl, 2L)
    expect_equal(batch_opencl[[1L]], solve(A, B), tolerance = 5e-6)
    expect_equal(batch_opencl[[2L]], solve(A, B), tolerance = 5e-6)
    expect_identical(length(fac@factor), 0L)
  })
})
