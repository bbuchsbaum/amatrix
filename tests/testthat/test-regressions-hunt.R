# Pinned regression tests for already-fixed silent-wrong-answer bugs surfaced
# by the bug-hunt rounds. Each test_that asserts the CORRECT post-fix behavior
# so a future regression fails loudly. CPU-only, exact / 1e-10 tolerance.
#
# Only bugs WITHOUT an existing dedicated pin are covered here; already-pinned
# bugs (e.g. amatrix-1ha, -fya, -x5o, -adu, -e97, -6lg's NA sibling -8li, -inz,
# -n5x, -p24, -03p, -3gv, ...) live in their own regression / bughunt files.

test_that("regression amatrix-6lg: dist/kernel never silently zero non-finite rows", {
  # Bug: a row containing NaN/Inf produced exact 0 distances (and an rbf kernel
  # value of 1) for every pair involving that row, silently contaminating
  # downstream ML. Correct post-fix behavior is to fail fast on non-finite
  # input. We assert the anti-regression rather than the exact policy: either
  # the call errors, or -- if a value is returned -- the contaminated row must
  # contain a non-finite entry (never a clean, all-finite 0 / 1 row).
  assert_not_silently_finite <- function(fun, X, bad_row) {
    res <- tryCatch(as.matrix(fun(adgeMatrix(X, preferred_backend = "cpu"))),
                    error = function(e) e)
    if (inherits(res, "error")) {
      succeed() # fail-fast on non-finite input is the current, acceptable policy
    } else {
      expect_false(
        all(is.finite(res[bad_row, ])),
        label = "non-finite row collapsed to finite distances/kernel values"
      )
    }
  }
  X_nan <- matrix(c(1, 2, NaN, 4, 5, 6), 3, 2)
  X_inf <- matrix(c(1, 2, Inf, 4, 5, 6), 3, 2)
  assert_not_silently_finite(function(A) dist_matrix(A, A), X_nan, 3L)
  assert_not_silently_finite(function(A) dist_matrix(A, A), X_inf, 3L)
  assert_not_silently_finite(
    function(A) kernel_matrix(A, kernel = "rbf", sigma = 1), X_nan, 3L
  )
})

test_that("regression amatrix-47w: Compare ops keep amatrix class, not bare Matrix logical", {
  # Bug: ==, !=, <, >, <=, >= on adgeMatrix/adgCMatrix demoted to the plain
  # Matrix-pkg logical classes (lgeMatrix / lgCMatrix), dropping all amatrix
  # slots and forcing a CPU round-trip in any downstream any()/all() pipeline.
  ops <- c("==", "!=", "<", ">", "<=", ">=")

  A_r <- matrix(c(1, 2, 3, 4), 2, 2)
  A <- adgeMatrix(A_r, preferred_backend = "cpu")
  for (op in ops) {
    r <- get(op)(A, 2)
    expect_false(
      identical(class(r)[[1L]], "lgeMatrix"),
      label = paste0("adgeMatrix ", op, " demoted to bare lgeMatrix")
    )
    expect_true(
      grepl("^adl", class(r)[[1L]]),
      label = paste0("adgeMatrix ", op, " lost amatrix logical class")
    )
    expect_equal(as.matrix(r), get(op)(A_r, 2), label = paste0("dense ", op))
  }

  # non-triangular sparse pattern -> dgCMatrix (avoids dtCMatrix coercion path)
  Sp <- Matrix::Matrix(matrix(c(0, 2, 3, 4), 2, 2), sparse = TRUE)
  As <- adgCMatrix(Sp, preferred_backend = "cpu")
  Sp_r <- as.matrix(Sp)
  for (op in ops) {
    r <- get(op)(As, 1)
    expect_false(
      identical(class(r)[[1L]], "lgCMatrix"),
      label = paste0("adgCMatrix ", op, " demoted to bare lgCMatrix")
    )
    expect_true(
      grepl("^adl", class(r)[[1L]]),
      label = paste0("adgCMatrix ", op, " lost amatrix logical class")
    )
    expect_equal(as.matrix(r), get(op)(Sp_r, 1), label = paste0("sparse ", op))
  }
})

test_that("regression amatrix-8oy: rbf self-kernel diagonal is exactly 1 on cpu", {
  # Bug: kernel_matrix(Y, Y, 'rbf') left f32 self-distance drift in the
  # diagonal (values in [0.999999, 1]) instead of exactly 1. On the cpu
  # (float64) path the diagonal must be exactly 1.
  set.seed(8)
  Y <- matrix(rnorm(20), 5, 4)
  K <- kernel_matrix(adgeMatrix(Y, preferred_backend = "cpu"), kernel = "rbf", sigma = 1)
  expect_equal(diag(as.matrix(K)), rep(1, 5), tolerance = 0)
})

test_that("regression amatrix-yb1: sparse rbind/cbind preserve adgCMatrix class and values", {
  # Bug: rbind/cbind on adgCMatrix could crash or demote the class. Once
  # constructed from a dgCMatrix, binding must return an adgCMatrix whose
  # materialized values match the base-R bind.
  Sp <- Matrix::Matrix(matrix(c(0, 2, 3, 4), 2, 2), sparse = TRUE)
  Xs <- adgCMatrix(Sp, preferred_backend = "cpu")
  Sp_r <- as.matrix(Sp)

  rb <- rbind(Xs, Xs)
  expect_s4_class(rb, "adgCMatrix")
  expect_equal(as.matrix(rb), rbind(Sp_r, Sp_r))

  cb <- cbind(Xs, Xs)
  expect_s4_class(cb, "adgCMatrix")
  expect_equal(as.matrix(cb), cbind(Sp_r, Sp_r))
})

test_that("regression amatrix-juq: row/col Sums/Means give correct values on adgeMatrix", {
  # amatrix-juq (duplicate root cause of amatrix-1ha): base reduction primitives
  # bypassed S4 dispatch. amatrix-1ha pins the fresh-attach dispatch path in a
  # child process; here we pin exact in-process numeric parity of the S4
  # reduction methods against base R.
  M <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3)
  A <- adgeMatrix(M, preferred_backend = "cpu")
  expect_equal(rowSums(A), rowSums(M), tolerance = 1e-10)
  expect_equal(colSums(A), colSums(M), tolerance = 1e-10)
  expect_equal(rowMeans(A), rowMeans(M), tolerance = 1e-10)
  expect_equal(colMeans(A), colMeans(M), tolerance = 1e-10)
})
