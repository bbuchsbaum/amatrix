# test-bughunt-coverage-gap.R
# amatrix-lcn, amatrix-jbs, amatrix-doj: exported functions with ZERO test coverage.
# These tests FAIL intentionally (via expect_true(FALSE, ...)) to force attention;
# replace with real assertions when implementations are verified.

# ── amatrix-lcn: batch 1 ─────────────────────────────────────────────────────

test_that("addmm: A + alpha*(B %*% C) returns correct result [amatrix-lcn]", {
  A <- matrix(1:4, 2, 2) * 1.0
  B <- matrix(c(1, 0, 0, 1), 2, 2)
  C <- matrix(c(2, 0, 0, 3), 2, 2)
  result <- addmm(adgeMatrix(A), B, C, alpha = 1.0, beta = 1.0)
  expected <- A + 1.0 * (B %*% C)
  expect_equal(as.matrix(result), expected, tolerance = 1e-10)
})

test_that("am_colargmax returns per-column argmax row index [amatrix-lcn]", {
  # am_colargmax(x)[j] = which.max(x[, j]) for each column j
  m <- matrix(c(1, 5, 3, 2, 4, 6), nrow = 2)
  x <- adgeMatrix(m)
  result <- am_colargmax(x)
  expected <- apply(m, 2, which.max)
  expect_equal(as.integer(result), as.integer(expected))
})

test_that("am_colargmin returns per-column argmin row index [amatrix-lcn]", {
  # am_colargmin(x)[j] = which.min(x[, j]) for each column j
  m <- matrix(c(5, 1, 3, 2, 4, 6), nrow = 2)
  x <- adgeMatrix(m)
  result <- am_colargmin(x)
  expected <- apply(m, 2, which.min)
  expect_equal(as.integer(result), as.integer(expected))
})

test_that("am_ewise_inplace modifies resident handle in place [amatrix-lcn]", {
  skip_if_not_installed("amatrix.mlx")
  m <- matrix(1:4 * 1.0, 2, 2)
  x <- adgeMatrix(m, preferred_backend = "mlx")
  h <- tryCatch(resident_handle(x), error = function(e) NULL)
  skip_if(is.null(h), "No GPU backend available for residency")
  am_ewise_inplace(h, 2.0, "+")
  result <- as.matrix(h)
  expect_equal(result, m + 2.0, tolerance = 1e-10)
})

test_that("am_rowargmax returns per-row argmax col index [amatrix-lcn]", {
  # am_rowargmax(x)[i] = which.max(x[i, ]) for each row i
  m <- matrix(c(1, 5, 3, 2, 4, 6), nrow = 2)
  # row1=(1,3,4)->col3, row2=(5,2,6)->col3
  x <- adgeMatrix(m)
  result <- am_rowargmax(x)
  expected <- apply(m, 1, which.max)
  expect_equal(as.integer(result), as.integer(expected))
})

test_that("am_rowargmin returns per-row argmin col index [amatrix-lcn]", {
  # am_rowargmin(x)[i] = which.min(x[i, ]) for each row i
  m <- matrix(c(5, 1, 3, 2, 6, 4), nrow = 2)
  # row1=(5,3,6)->col2, row2=(1,2,4)->col1
  x <- adgeMatrix(m)
  result <- am_rowargmin(x)
  expected <- apply(m, 1, which.min)
  expect_equal(as.integer(result), as.integer(expected))
})

test_that("am_scatter_mean computes group means correctly [amatrix-lcn]", {
  m <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)
  x <- adgeMatrix(m)
  labels <- c(1L, 1L, 2L)
  result <- am_scatter_mean(x, labels, K = 2L)
  expected <- rbind(colMeans(m[1:2, , drop = FALSE]), m[3, ])
  expect_equal(as.matrix(result), expected, tolerance = 1e-10)
})

test_that("am_sweep_inplace applies sweep in place on resident handle [amatrix-lcn]", {
  skip_if_not_installed("amatrix.mlx")
  m <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2) * 1.0
  x <- adgeMatrix(m, preferred_backend = "mlx")
  h <- tryCatch(resident_handle(x), error = function(e) NULL)
  skip_if(is.null(h), "No GPU backend available for residency")
  am_sweep_inplace(h, MARGIN = 1L, STATS = c(10.0, 20.0), FUN = "+")
  result <- as.matrix(h)
  expected <- sweep(m, 1, c(10.0, 20.0), "+")
  expect_equal(result, expected, tolerance = 1e-10)
})

# ── amatrix-jbs: batch 2 ─────────────────────────────────────────────────────

test_that("amatrix_backend_precision_modes returns character vector for cpu [amatrix-jbs]", {
  modes <- amatrix_backend_precision_modes("cpu")
  expect_true(is.character(modes))
  expect_true(length(modes) >= 1L)
  # cpu backend uses "strict" and "fast" precision labels (not "double")
  expect_true("strict" %in% modes)
})

test_that("amatrix_cache_max_size returns numeric [amatrix-jbs]", {
  sz <- amatrix_cache_max_size()
  expect_true(is.numeric(sz) && length(sz) == 1L)
  expect_true(sz >= 0)
})

test_that("amatrix_calibration_info runs without error (NULL before calibration) [amatrix-jbs]", {
  info <- amatrix_calibration_info()
  expect_true(is.null(info) || is.list(info))
})

test_that("amatrix_explain returns invisibly without error for cpu matrix [amatrix-jbs]", {
  m <- adgeMatrix(matrix(1:4 * 1.0, 2, 2))
  expect_no_error(amatrix_explain(m, "matmul"))
})

test_that("amatrix_gc runs without error and returns invisibly [amatrix-jbs]", {
  expect_no_error(amatrix_gc())
})

test_that("amatrix_memory_stats returns a list with residency and model_cache fields [amatrix-jbs]", {
  stats <- amatrix_memory_stats()
  expect_true(is.list(stats))
  expect_true("residency" %in% names(stats))
  expect_true("model_cache" %in% names(stats))
  expect_true("resident_objects" %in% names(stats$residency))
})

test_that("amatrix_set_cache_max_size updates and amatrix_cache_max_size reflects it [amatrix-jbs]", {
  old <- amatrix_cache_max_size()
  on.exit(amatrix_set_cache_max_size(old), add = TRUE)
  amatrix_set_cache_max_size(5L)
  expect_equal(amatrix_cache_max_size(), 5L)
})

test_that("amatrix_warm runs without error (no GPU available path) [amatrix-jbs]", {
  expect_no_error(amatrix_warm(quiet = TRUE))
})

# ── amatrix-doj: batch 3 ─────────────────────────────────────────────────────

test_that("as.matrix.KronMatrix round-trips to dense matrix [amatrix-doj]", {
  A <- matrix(1:4 * 1.0, 2, 2)
  B <- matrix(1:9 * 1.0, 3, 3)
  K <- kron_matrix(A, B)
  result <- as.matrix(K)
  expected <- kronecker(A, B)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("as.matrix.resident_handle materializes to host matrix [amatrix-doj]", {
  skip_if_not_installed("amatrix.mlx")
  m <- matrix(1:6 * 1.0, 2, 3)
  x <- adgeMatrix(m, preferred_backend = "mlx")
  h <- tryCatch(resident_handle(x), error = function(e) NULL)
  skip_if(is.null(h), "No GPU backend available for residency")
  result <- as.matrix(h)
  expect_equal(result, m, tolerance = 1e-10)
})

test_that("colsums returns column sums matching base R [amatrix-doj]", {
  m <- matrix(1:6 * 1.0, 2, 3)
  x <- adgeMatrix(m)
  result <- colsums(x)
  expect_equal(result, colSums(m), tolerance = 1e-10)
})

test_that("irlba_native returns nv singular values [amatrix-doj]", {
  skip_if_not_installed("irlba")
  set.seed(42)
  m <- matrix(rnorm(50), 10, 5)
  result <- irlba_native(m, nv = 3L)
  expect_true(!is.null(result$d))
  expect_equal(length(result$d), 3L)
  ref <- base::svd(m)$d[1:3]
  expect_equal(result$d, ref, tolerance = 1e-6)
})

test_that("ncol.resident_handle returns number of columns [amatrix-doj]", {
  skip_if_not_installed("amatrix.mlx")
  m <- matrix(1:6 * 1.0, 2, 3)
  x <- adgeMatrix(m, preferred_backend = "mlx")
  h <- tryCatch(resident_handle(x), error = function(e) NULL)
  skip_if(is.null(h), "No GPU backend available for residency")
  expect_equal(ncol(h), 3L)
})

test_that("nrow.resident_handle returns number of rows [amatrix-doj]", {
  skip_if_not_installed("amatrix.mlx")
  m <- matrix(1:6 * 1.0, 2, 3)
  x <- adgeMatrix(m, preferred_backend = "mlx")
  h <- tryCatch(resident_handle(x), error = function(e) NULL)
  skip_if(is.null(h), "No GPU backend available for residency")
  expect_equal(nrow(h), 2L)
})

test_that("pairwise_sqdist_argmin returns nearest centroid index [amatrix-doj]", {
  # X: n×p matrix. Rows are points, cols are dims.
  # Ct: p×K matrix. Rows are dims, cols are centroids.
  # C1=(0,0), C2=(10,0): Ct = cbind(c(0,0), c(10,0))
  # P1=(1,0) nearest C1; P2=(9,0) nearest C2 -> expected c(1L, 2L)
  # BUG: CPU path uses rep(c_norms, each=nrow) instead of rep(c_norms, nrow),
  # causing incorrect centroid distance broadcast -> wrong assignments.
  X <- matrix(c(1, 9, 0, 0), nrow = 2, ncol = 2, byrow = FALSE)  # P1=(1,0),P2=(9,0)
  Ct <- cbind(c(0, 0), c(10, 0))  # C1=(0,0), C2=(10,0)
  result <- pairwise_sqdist_argmin(X, Ct)
  expect_equal(as.integer(result), c(1L, 2L))
})

test_that("rowsums returns row sums matching base R [amatrix-doj]", {
  m <- matrix(1:6 * 1.0, 2, 3)
  x <- adgeMatrix(m)
  result <- rowsums(x)
  expect_equal(result, rowSums(m), tolerance = 1e-10)
})
