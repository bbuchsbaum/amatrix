# Extracted from test-chol-factor.R:14

# prequel ----------------------------------------------------------------------
make_spd <- function(n, seed = 1L) {
  set.seed(seed)
  A <- matrix(rnorm(n * n), n, n)
  crossprod(A) + diag(n)
}

# test -------------------------------------------------------------------------
M <- make_spd(8L, seed = 42L)
X <- as_adgeMatrix(M)
fac <- am_chol_factor(X)
expect_s4_class(fac, "amChol")
R <- as.matrix(fac)
