# Regressions for the fourth batch of pre-CRAN review fixes.

test_that("dist_matrix cosine is a distance, not a similarity", {
  set.seed(1)
  X <- matrix(rnorm(20), 5)
  nx <- sqrt(rowSums(X^2))
  ref <- 1 - tcrossprod(X) / outer(nx, nx)
  diag(ref) <- 0

  D <- as.matrix(dist_matrix(X, method = "cosine"))
  expect_equal(D, ref, ignore_attr = TRUE, tolerance = 1e-10)
  expect_true(all(diag(D) == 0))
  expect_true(all(D >= 0 & D <= 2))

  Dxy <- as.matrix(dist_matrix(adgeMatrix(X), X[1:2, ], method = "cosine"))
  expect_equal(Dxy, 1 - tcrossprod(X, X[1:2, ]) / outer(nx, nx[1:2]),
               ignore_attr = TRUE, tolerance = 1e-10)
  # The similarity remains available as a kernel.
  expect_equal(as.matrix(kernel_matrix(X, kernel = "cosine")), 1 - ref,
               ignore_attr = TRUE, tolerance = 1e-10)
})

test_that("products with a plain matrix operand match base above the in-place threshold", {
  set.seed(2)
  n <- 520L
  M <- matrix(rnorm(n * n), n, dimnames = list(paste0("r", seq_len(n)), NULL))
  A <- adgeMatrix(M)
  expect_s4_class(A %*% M, "adgeMatrix")
  expect_equal(as.matrix(A %*% M), M %*% M)
  expect_identical(dimnames(A %*% M), dimnames(M %*% M))
  expect_equal(as.matrix(M %*% A), M %*% M)
  expect_equal(as.matrix(crossprod(A, M)), crossprod(M), ignore_attr = TRUE)
})

test_that("the fallback log is bounded", {
  state <- amatrix:::.amatrix_state
  old_log <- state$fallback_log
  on.exit(state$fallback_log <- old_log, add = TRUE)
  amatrix_fallback_log_reset()
  for (i in seq_len(1100L)) {
    amatrix:::.amatrix_log_fallback(op = "matmul", backend = "fake", reason = as.character(i))
  }
  log <- amatrix_fallback_log()
  expect_lte(nrow(log), 1000L)
  expect_identical(utils::tail(log$reason, 1L), "1100")
})

test_that("amatrix_compile_product leaves the defer_host option as it found it", {
  withr::local_options(amatrix.defer_host = NULL)
  x <- adgeMatrix(matrix(rnorm(20), 5, 4))
  plan <- amatrix_compile_product(x, op = "matmul")
  plan(matrix(rnorm(8), 4, 2))
  expect_null(getOption("amatrix.defer_host"))
})
