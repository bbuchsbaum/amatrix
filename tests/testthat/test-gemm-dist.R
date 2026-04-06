# test-gemm-dist.R — tests for .am_gemm() and am_dist(tile_size)

# ── .am_gemm ──────────────────────────────────────────────────────────────────

test_that(".am_gemm(A, B) matches A %*% B", {
  set.seed(51)
  A_r <- matrix(rnorm(8 * 5), 8, 5)
  B_r <- matrix(rnorm(5 * 7), 5, 7)
  A <- adgeMatrix(A_r, preferred_backend = "cpu")
  B <- adgeMatrix(B_r, preferred_backend = "cpu")

  res  <- amatrix:::.am_gemm(A, B)
  ref  <- A_r %*% B_r
  expect_true(inherits(res, "adgeMatrix"))
  expect_equal(dim(res), c(8L, 7L))
  expect_equal(as.matrix(res), ref, tolerance = 1e-10)
})

test_that(".am_gemm(A, B, alpha=2) matches 2 * A %*% B", {
  set.seed(52)
  A_r <- matrix(rnorm(6 * 4), 6, 4)
  B_r <- matrix(rnorm(4 * 3), 4, 3)
  A <- adgeMatrix(A_r, preferred_backend = "cpu")
  B <- adgeMatrix(B_r, preferred_backend = "cpu")

  res <- amatrix:::.am_gemm(A, B, alpha = 2.0)
  expect_equal(as.matrix(res), 2 * A_r %*% B_r, tolerance = 1e-10)
})

test_that(".am_gemm(A, B, C, alpha=2, beta=0.5) matches 2*A%*%B + 0.5*C", {
  set.seed(53)
  A_r <- matrix(rnorm(5 * 4), 5, 4)
  B_r <- matrix(rnorm(4 * 6), 4, 6)
  C_r <- matrix(rnorm(5 * 6), 5, 6)
  A <- adgeMatrix(A_r, preferred_backend = "cpu")
  B <- adgeMatrix(B_r, preferred_backend = "cpu")
  C <- adgeMatrix(C_r, preferred_backend = "cpu")

  res <- amatrix:::.am_gemm(A, B, C, alpha = 2.0, beta = 0.5)
  ref <- 2 * A_r %*% B_r + 0.5 * C_r
  expect_equal(as.matrix(res), ref, tolerance = 1e-10)
})

test_that(".am_gemm transA=TRUE routes through am_crossprod", {
  set.seed(54)
  A_r <- matrix(rnorm(10 * 4), 10, 4)   # 10×4; transA gives 4×10
  B_r <- matrix(rnorm(10 * 3), 10, 3)   # 10×3
  A <- adgeMatrix(A_r, preferred_backend = "cpu")
  B <- adgeMatrix(B_r, preferred_backend = "cpu")

  res <- amatrix:::.am_gemm(A, B, transA = TRUE)    # t(A) %*% B → 4×3
  ref <- t(A_r) %*% B_r
  expect_equal(dim(res), c(4L, 3L))
  expect_equal(as.matrix(res), ref, tolerance = 1e-10)
})

test_that(".am_gemm transB=TRUE routes through am_tcrossprod", {
  set.seed(55)
  A_r <- matrix(rnorm(6 * 9), 6, 9)   # 6×9
  B_r <- matrix(rnorm(4 * 9), 4, 9)   # 4×9; transB gives 9×4
  A <- adgeMatrix(A_r, preferred_backend = "cpu")
  B <- adgeMatrix(B_r, preferred_backend = "cpu")

  res <- amatrix:::.am_gemm(A, B, transB = TRUE)    # A %*% t(B) → 6×4
  ref <- A_r %*% t(B_r)
  expect_equal(dim(res), c(6L, 4L))
  expect_equal(as.matrix(res), ref, tolerance = 1e-10)
})

# ── am_dist tiled ─────────────────────────────────────────────────────────────

test_that("am_dist tiled self-distance matches untiled on small matrix", {
  set.seed(61)
  X_r <- matrix(rnorm(50 * 8), 50, 8)

  ref   <- am_dist(X_r, method = "euclidean")
  tiled <- am_dist(X_r, method = "euclidean", tile_size = 15L)

  expect_equal(dim(tiled), c(50L, 50L))
  expect_equal(tiled, ref, tolerance = 1e-6)
})

test_that("am_dist tiled sqeuclidean matches untiled", {
  set.seed(62)
  X_r <- matrix(rnorm(40 * 5), 40, 5)

  ref   <- am_dist(X_r, method = "sqeuclidean")
  tiled <- am_dist(X_r, method = "sqeuclidean", tile_size = 12L)

  expect_equal(tiled, ref, tolerance = 1e-6)
})

test_that("am_dist tiled diagonal is zero", {
  set.seed(63)
  X_r <- matrix(rnorm(30 * 6), 30, 6)
  tiled <- am_dist(X_r, method = "euclidean", tile_size = 11L)
  expect_equal(diag(tiled), rep(0, 30))
})

test_that("am_dist tiled result is symmetric", {
  set.seed(64)
  X_r <- matrix(rnorm(25 * 4), 25, 4)
  tiled <- am_dist(X_r, method = "euclidean", tile_size = 8L)
  expect_equal(tiled, t(tiled), tolerance = 1e-10)
})

test_that("am_dist tiled X,Y matches untiled", {
  set.seed(65)
  X_r <- matrix(rnorm(20 * 5), 20, 5)
  Y_r <- matrix(rnorm(15 * 5), 15, 5)

  ref   <- am_dist(X_r, Y_r, method = "euclidean")
  tiled <- am_dist(X_r, Y_r, method = "euclidean", tile_size = 7L)

  expect_equal(dim(tiled), c(20L, 15L))
  expect_equal(tiled, ref, tolerance = 1e-6)
})
