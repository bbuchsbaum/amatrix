# Matrix's unary minus / round / signif / zapsmall methods update `x@x` in
# place on a copy and so keep every amatrix slot, including `object_id`. The
# result then aliased the identity-keyed caches of its source matrix and
# am_chol_factor()/am_lm_fit() silently returned factors of the *original* data.

test_that("value-changing unary and rounding ops mint a new object identity", {
  m <- matrix(sin(1:16) * 3, 4, dimnames = list(letters[1:4], LETTERS[1:4]))
  A <- adgeMatrix(m)
  sp <- Matrix::sparseMatrix(
    i = c(1, 3, 2, 4, 1), j = c(1, 1, 2, 3, 4), x = c(1.26, -2.51, 3.14159, 0.5, -7.77),
    dims = c(4, 4), dimnames = list(letters[1:4], LETTERS[1:4])
  )
  B <- adgCMatrix(sp)

  cases <- list(
    list(function(x) -x, "neg"),
    list(function(x) round(x, 1), "round1"),
    list(function(x) round(x), "round0"),
    list(function(x) signif(x, 2), "signif"),
    list(function(x) zapsmall(x), "zapsmall")
  )
  for (case in cases) {
    f <- case[[1L]]
    for (pair in list(list(A, m, "adgeMatrix"), list(B, sp, "adgCMatrix"))) {
      got <- f(pair[[1L]])
      ref <- f(pair[[2L]])
      expect_s4_class(got, pair[[3L]])
      expect_false(identical(got@object_id, pair[[1L]]@object_id), info = case[[2L]])
      expect_equal(as.matrix(got), as.matrix(ref), info = case[[2L]])
      expect_identical(dimnames(got), dimnames(ref), info = case[[2L]])
    }
  }
})

test_that("chol_factor and lm_fit do not reuse a cached factor after round()", {
  set.seed(3)
  S <- adgeMatrix(crossprod(matrix(rnorm(40), 10)) + diag(4))
  am_chol_factor(S)
  S2 <- round(S, 0)
  expect_equal(
    as.matrix(am_chol_factor(S2)),
    chol(as.matrix(S2)),
    ignore_attr = TRUE
  )

  X <- adgeMatrix(cbind(1, matrix(rnorm(60), 20)))
  y <- rnorm(20)
  am_lm_fit(X, y)
  X2 <- -X
  expect_equal(
    as.numeric(as.matrix(coef(am_lm_fit(X2, y)))),
    as.numeric(qr.coef(qr(as.matrix(X2)), y)),
    tolerance = 1e-8
  )
})

test_that("model-cache entries are released when their source matrix is collected", {
  cache <- amatrix:::.amatrix_state$model_cache
  gc()
  before <- length(ls(cache))
  set.seed(11)
  for (i in 1:8) {
    X <- matrix(rnorm(120), 30)
    am_lm_fit(X, rnorm(30))
    am_chol_factor(adgeMatrix(crossprod(X) + diag(4)))
  }
  rm(X)
  gc()
  gc()
  expect_lte(length(ls(cache)), before + 2L)
})
