# test-dispatch-primitives.R
#
# Dispatch gap regression tests. Each test verifies that the expected S4 method
# fires (i.e., the GPU path is taken) rather than falling through to base R
# which would coerce adgeMatrix to a plain matrix.
#
# Tests are deterministic and run against the cpu backend only; GPU backends
# are tested by the cross-backend conformance harness.

.disp_tol <- 1e-10

test_that("crossprod(matrix, adgeMatrix) routes to GPU path, not base coercion", {
  set.seed(31)
  A_r <- matrix(rnorm(12 * 5), 12, 5)   # 12×5 plain matrix
  B_r <- matrix(rnorm(12 * 7), 12, 7)   # 12×7 → B adgeMatrix: m=12 cols=7
  B   <- adgeMatrix(B_r, preferred_backend = "cpu")

  # crossprod(A, B) = t(A) %*% B  →  5×7
  result <- crossprod(A_r, B)

  expect_true(inherits(result, "adgeMatrix"),
              info = "crossprod(matrix, adgeMatrix) must return an adgeMatrix, not a plain matrix")
  expect_equal(as.matrix(result), base::crossprod(A_r, B_r),
               tolerance = .disp_tol,
               label = "crossprod(matrix, adgeMatrix) numerics")
})

test_that("tcrossprod(matrix, adgeMatrix) routes to GPU path, not base coercion", {
  set.seed(32)
  A_r <- matrix(rnorm(6 * 9), 6, 9)    # 6×9 plain matrix
  B_r <- matrix(rnorm(4 * 9), 4, 9)    # 4×9 → B adgeMatrix
  B   <- adgeMatrix(B_r, preferred_backend = "cpu")

  # tcrossprod(A, B) = A %*% t(B)  →  6×4
  result <- tcrossprod(A_r, B)

  expect_true(inherits(result, "adgeMatrix"),
              info = "tcrossprod(matrix, adgeMatrix) must return an adgeMatrix, not a plain matrix")
  expect_equal(as.matrix(result), base::tcrossprod(A_r, B_r),
               tolerance = .disp_tol,
               label = "tcrossprod(matrix, adgeMatrix) numerics")
})

test_that("crossprod(numeric, adgeMatrix) routes to GPU path", {
  set.seed(33)
  v   <- rnorm(10)                       # length-10 numeric vector
  B_r <- matrix(rnorm(10 * 6), 10, 6)   # 10×6 adgeMatrix
  B   <- adgeMatrix(B_r, preferred_backend = "cpu")

  # crossprod(v, B) = t(as.matrix(v)) %*% B  →  1×6
  result <- crossprod(v, B)

  expect_true(inherits(result, "adgeMatrix"),
              info = "crossprod(numeric, adgeMatrix) must return an adgeMatrix")
  expect_equal(as.matrix(result), base::crossprod(matrix(v, ncol = 1L), B_r),
               tolerance = .disp_tol,
               label = "crossprod(numeric, adgeMatrix) numerics")
})

test_that("matrix %*% adgeMatrix routes to GPU path (pre-existing method)", {
  set.seed(34)
  A_r <- matrix(rnorm(3 * 8), 3, 8)    # 3×8 plain matrix
  B_r <- matrix(rnorm(8 * 5), 8, 5)    # 8×5 adgeMatrix
  B   <- adgeMatrix(B_r, preferred_backend = "cpu")

  result <- A_r %*% B

  expect_true(inherits(result, "adgeMatrix"),
              info = "matrix %*% adgeMatrix must return adgeMatrix")
  expect_equal(as.matrix(result), A_r %*% B_r,
               tolerance = .disp_tol,
               label = "matrix %*% adgeMatrix numerics")
})

test_that("numeric %*% adgeMatrix routes to GPU path (pre-existing method)", {
  set.seed(35)
  v   <- rnorm(7)
  B_r <- matrix(rnorm(7 * 4), 7, 4)
  B   <- adgeMatrix(B_r, preferred_backend = "cpu")

  # numeric (length 7) %*% adgeMatrix (7×4) → 1×4 (drop to vector by base)
  result <- v %*% B

  # The numeric %*% adgeMatrix method returns crossprod(B, matrix(v)) = t(B) %*% v  (4×1),
  # and base::`%*%`(v, B_r) = 1×4. Both give the same absolute values.
  ref <- drop(t(B_r) %*% v)
  expect_equal(drop(as.matrix(result)), ref,
               tolerance = .disp_tol,
               label = "numeric %*% adgeMatrix numerics")
})

test_that("adgeMatrix crossprod self/cross shape and type are correct", {
  set.seed(36)
  A_r <- matrix(rnorm(20 * 4), 20, 4)
  C_r <- matrix(rnorm(20 * 6), 20, 6)
  A   <- adgeMatrix(A_r, preferred_backend = "cpu")
  C   <- adgeMatrix(C_r, preferred_backend = "cpu")

  cp_self <- crossprod(A)
  expect_equal(dim(cp_self), c(4L, 4L))
  expect_true(inherits(cp_self, "adgeMatrix"))
  expect_equal(as.matrix(cp_self), base::crossprod(A_r), tolerance = .disp_tol)

  cp_cross <- crossprod(A, C)
  expect_equal(dim(cp_cross), c(4L, 6L))
  expect_equal(as.matrix(cp_cross), base::crossprod(A_r, C_r), tolerance = .disp_tol)
})

test_that("adgeMatrix tcrossprod self/cross shape and type are correct", {
  set.seed(37)
  A_r <- matrix(rnorm(15 * 5), 15, 5)
  C_r <- matrix(rnorm(12 * 5), 12, 5)
  A   <- adgeMatrix(A_r, preferred_backend = "cpu")
  C   <- adgeMatrix(C_r, preferred_backend = "cpu")

  tc_self <- tcrossprod(A)
  expect_equal(dim(tc_self), c(15L, 15L))
  expect_true(inherits(tc_self, "adgeMatrix"))
  expect_equal(as.matrix(tc_self), base::tcrossprod(A_r), tolerance = .disp_tol)

  tc_cross <- tcrossprod(A, C)
  expect_equal(dim(tc_cross), c(15L, 12L))
  expect_equal(as.matrix(tc_cross), base::tcrossprod(A_r, C_r), tolerance = .disp_tol)
})

# --- zero-copy t() fast-path tests -------------------------------------------

test_that("t() of a non-resident adgeMatrix has empty src_id", {
  set.seed(41)
  A_r <- matrix(rnorm(8 * 3), 8, 3)
  A   <- adgeMatrix(A_r, preferred_backend = "cpu")

  tA <- t(A)

  expect_true(inherits(tA, "adgeMatrix"))
  expect_equal(dim(tA), c(3L, 8L))
  expect_true(!nzchar(tA@src_id),
              info = "non-resident t() must not set src_id")
  expect_equal(as.matrix(tA), t(A_r), tolerance = .disp_tol)
})

test_that("t() of a resident adgeMatrix sets src_id and swaps dims", {
  counter <- new.env(parent = emptyenv())
  backend <- make_recording_backend(counter)

  with_registered_backend("mock-t-src", backend, {
    set.seed(40)
    A_r <- matrix(rnorm(12 * 5), 12, 5)
    A   <- adgeMatrix(A_r, preferred_backend = "mock-t-src")

    # trigger residency
    invisible(A %*% matrix(rnorm(5 * 2), 5, 2))

    tA <- t(A)

    expect_true(inherits(tA, "adgeMatrix"))
    expect_equal(dim(tA), c(5L, 12L))
    expect_equal(tA@src_id, A@object_id,
                 info = "transposed view must carry original object_id in src_id")
    expect_equal(as.matrix(tA), t(A_r), tolerance = .disp_tol,
                 label = "t() host data must be correct")
  })
})

test_that("t(A) %*% B routes to crossprod_resident (not matmul_resident)", {
  counter <- new.env(parent = emptyenv())
  backend <- make_recording_backend(counter,
    resident_supported_ops = c("matmul", "crossprod", "tcrossprod", "ewise"))

  with_registered_backend("mock-t-cp", backend, {
    set.seed(42)
    A_r <- matrix(rnorm(12 * 5), 12, 5)   # 12x5
    B_r <- matrix(rnorm(12 * 7), 12, 7)   # 12x7
    A   <- adgeMatrix(A_r, preferred_backend = "mock-t-cp")
    B   <- adgeMatrix(B_r, preferred_backend = "mock-t-cp")

    # trigger residency for both
    invisible(A %*% matrix(rnorm(5), 5, 1))
    invisible(B %*% matrix(rnorm(7), 7, 1))

    counter$crossprod_resident <- 0L
    counter$matmul_resident    <- 0L

    result <- t(A) %*% B   # 5x7 via crossprod fast path

    expect_true(inherits(result, "adgeMatrix"))
    expect_equal(dim(result), c(5L, 7L))
    expect_equal(as.matrix(result), base::crossprod(A_r, B_r),
                 tolerance = .disp_tol,
                 label = "t(A) %*% B numerics via src_id fast path")
    expect_equal(counter$crossprod_resident, 1L,
                 info = "fast path must call crossprod_resident exactly once")
    expect_equal(counter$matmul_resident %||% 0L, 0L,
                 info = "fast path must not call matmul_resident")
  })
})

test_that("A %*% t(B) routes to tcrossprod_resident (not matmul_resident)", {
  counter <- new.env(parent = emptyenv())
  backend <- make_recording_backend(counter,
    resident_supported_ops = c("matmul", "crossprod", "tcrossprod", "ewise"))

  with_registered_backend("mock-t-tcp", backend, {
    set.seed(43)
    A_r <- matrix(rnorm(6 * 9), 6, 9)    # 6x9
    B_r <- matrix(rnorm(4 * 9), 4, 9)    # 4x9
    A   <- adgeMatrix(A_r, preferred_backend = "mock-t-tcp")
    B   <- adgeMatrix(B_r, preferred_backend = "mock-t-tcp")

    invisible(A %*% matrix(rnorm(9), 9, 1))
    invisible(B %*% matrix(rnorm(9), 9, 1))

    counter$tcrossprod_resident <- 0L
    counter$matmul_resident     <- 0L

    result <- A %*% t(B)   # 6x4 via tcrossprod fast path

    expect_true(inherits(result, "adgeMatrix"))
    expect_equal(dim(result), c(6L, 4L))
    expect_equal(as.matrix(result), base::tcrossprod(A_r, B_r),
                 tolerance = .disp_tol,
                 label = "A %*% t(B) numerics via src_id fast path")
    expect_equal(counter$tcrossprod_resident, 1L,
                 info = "fast path must call tcrossprod_resident exactly once")
    expect_equal(counter$matmul_resident %||% 0L, 0L,
                 info = "fast path must not call matmul_resident")
  })
})

test_that("t(A) %*% B falls back correctly when source is no longer resident", {
  counter <- new.env(parent = emptyenv())
  backend <- make_recording_backend(counter)

  with_registered_backend("mock-t-evict", backend, {
    set.seed(44)
    A_r <- matrix(rnorm(8 * 4), 8, 4)
    B_r <- matrix(rnorm(8 * 3), 8, 3)
    A   <- adgeMatrix(A_r, preferred_backend = "mock-t-evict")
    B   <- adgeMatrix(B_r, preferred_backend = "mock-t-evict")

    invisible(A %*% matrix(rnorm(4), 4, 1))

    tA <- t(A)   # src_id set to A@object_id

    # Simulate eviction by dropping A's resident binding
    amatrix:::.amatrix_drop_resident_binding(A)

    result <- tA %*% B   # src no longer live; falls back to normal matmul path

    expect_true(inherits(result, "adgeMatrix"))
    expect_equal(as.matrix(result), base::crossprod(A_r, B_r),
                 tolerance = .disp_tol,
                 label = "fallback when src evicted must still be numerically correct")
  })
})
