make_spd_regression_case <- function(n, seed = 1L) {
  set.seed(seed)
  a <- matrix(rnorm(n * n), nrow = n, ncol = n)
  crossprod(a) + diag(n)
}

test_that("chol_factor regression: basic SPD input produces an amChol factor", {
  m <- make_spd_regression_case(8L, seed = 42L)
  x <- as_adgeMatrix(m)

  fac <- chol_factor(x)

  expect_s4_class(fac, "amChol")
  expect_equal(as.matrix(fac), chol(m), tolerance = 1e-10)
})
