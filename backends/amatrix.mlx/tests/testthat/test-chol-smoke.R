test_that("mlx chol bridge and factor path smoke test on native backend", {
  skip_if_not(amatrix_mlx_native_available(), "mlx native backend not available")

  old <- getOption("amatrix.mlx.available")
  options(amatrix.mlx.available = TRUE)
  on.exit(options(amatrix.mlx.available = old), add = TRUE)

  set.seed(20260406L)
  z <- matrix(rnorm(160L * 160L), nrow = 160L, ncol = 160L)
  spd <- crossprod(z) + diag(0.5, 160L)
  rhs <- matrix(rnorm(160L * 32L), nrow = 160L, ncol = 32L)

  bridge_fac <- .Call("amatrix_mlx_chol_factor_bridge", spd, PACKAGE = "amatrix.mlx")
  bridge_fac[lower.tri(bridge_fac)] <- 0
  bridge_sol <- amatrix.mlx:::amatrix_mlx_chol_solve_factor(bridge_fac, rhs)

  factor_fit <- amatrix::chol_factor(
    amatrix::adgeMatrix(spd, preferred_backend = "mlx", precision = "fast")
  )
  factor_sol <- amatrix::chol_solve(factor_fit, rhs)

  ref_sol <- base::solve(spd, rhs)
  by_col <- vapply(
    seq_len(ncol(rhs)),
    function(j) {
      drop(amatrix.mlx:::amatrix_mlx_chol_solve_factor(bridge_fac, rhs[, j, drop = FALSE]))
    },
    numeric(nrow(rhs))
  )
  frob_norm <- function(x) sqrt(sum(x * x))
  bridge_recon_rel <- frob_norm(crossprod(bridge_fac) - spd) / frob_norm(spd)
  bridge_ref_rel <- frob_norm(bridge_sol - ref_sol) / frob_norm(ref_sol)
  bridge_resid_rel <- frob_norm(spd %*% bridge_sol - rhs) / frob_norm(rhs)
  bridge_batched_rel <- frob_norm(bridge_sol - by_col) / frob_norm(by_col)
  # @factor may be empty when the factor lives resident (factor_obj); the
  # public accessor as.matrix() materializes the dense upper factor.
  factor_recon_rel <- frob_norm(crossprod(as.matrix(factor_fit)) - spd) / frob_norm(spd)
  factor_ref_rel <- frob_norm(factor_sol - ref_sol) / frob_norm(ref_sol)
  factor_resid_rel <- frob_norm(spd %*% factor_sol - rhs) / frob_norm(rhs)

  expect_true(all(is.finite(bridge_sol)))
  expect_lt(bridge_recon_rel, 5e-5)
  expect_lt(bridge_ref_rel, 5e-5)
  expect_lt(bridge_resid_rel, 5e-5)
  expect_lt(bridge_batched_rel, 5e-5)

  expect_s4_class(factor_fit, "amChol")
  expect_identical(factor_fit@backend, "mlx")
  expect_true(all(is.finite(factor_fit@factor)))
  expect_true(all(is.finite(factor_sol)))
  expect_lt(factor_recon_rel, 5e-5)
  expect_lt(factor_ref_rel, 5e-5)
  expect_lt(factor_resid_rel, 5e-5)
})
