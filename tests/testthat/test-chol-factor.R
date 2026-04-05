make_spd <- function(n, seed = 1L) {
  set.seed(seed)
  A <- matrix(rnorm(n * n), n, n)
  crossprod(A) + diag(n)
}

test_that("am_chol_factor returns correct upper triangular factor", {
  M <- make_spd(8L, seed = 42L)
  X <- as_adgeMatrix(M)

  fac <- am_chol_factor(X)
  expect_s4_class(fac, "amChol")

  R <- as.matrix(fac)
  R_ref <- chol(as.matrix(M))

  expect_equal(R, R_ref, tolerance = 1e-10)
  expect_equal(t(R) %*% R, as.matrix(M), tolerance = 1e-10)
  # upper triangular
  expect_true(all(abs(R[lower.tri(R)]) < 1e-12))
})

test_that("am_chol_solve matches base::solve for single RHS vector", {
  M <- make_spd(10L, seed = 7L)
  X <- as_adgeMatrix(M)
  fac <- am_chol_factor(X)

  b <- rnorm(10L)
  sol <- am_chol_solve(fac, b)
  ref <- base::solve(M, b)

  expect_equal(as.numeric(sol), as.numeric(ref), tolerance = 1e-10)
})

test_that("am_chol_solve matches base::solve for multi-column B (k=1, 10, 100)", {
  n <- 12L
  M <- make_spd(n, seed = 11L)
  X <- as_adgeMatrix(M)
  fac <- am_chol_factor(X)

  for (k in c(1L, 10L, 100L)) {
    set.seed(k)
    B <- matrix(rnorm(n * k), n, k)
    sol <- am_chol_solve(fac, B)
    ref <- base::solve(M, B)
    expect_equal(dim(sol), c(n, k))
    expect_lt(max(abs(sol - ref)), 1e-10)
  }
})

test_that("am_chol_factor reuses cache on second call", {
  M <- make_spd(6L, seed = 3L)
  X <- as_adgeMatrix(M)

  fac1 <- am_chol_factor(X)
  fac2 <- am_chol_factor(X)

  expect_identical(fac1@factor, fac2@factor)
  expect_identical(fac1, fac2)
})

test_that("amChol show method runs without error", {
  M <- make_spd(4L, seed = 5L)
  X <- as_adgeMatrix(M)
  fac <- am_chol_factor(X)

  expect_output(show(fac), "amChol")
})
