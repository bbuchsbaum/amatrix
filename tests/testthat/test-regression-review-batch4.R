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

test_that("sweep() dispatches on amatrix objects and keeps sparse scaling sparse", {
  m <- matrix(c(1, 0, 3, 0, 5, 6), 2, 3)
  A <- adgeMatrix(m)
  S <- adgCMatrix(Matrix::Matrix(m, sparse = TRUE))
  hv <- function(x) if (methods::is(x, "Matrix")) as.matrix(x) else x
  for (margin in 1:2) {
    stats <- if (margin == 1L) c(2, 4) else c(1, 2, 4)
    for (fun in c("-", "+", "*", "/")) {
      ref <- base::sweep(m, margin, stats, fun)
      expect_equal(hv(sweep(A, margin, stats, fun)), ref, ignore_attr = TRUE, info = paste(margin, fun))
      expect_equal(hv(sweep(S, margin, stats, fun)), ref, ignore_attr = TRUE, info = paste(margin, fun))
    }
  }
  # Default operation is "-", as in base, for both sweep() and am_sweep().
  expect_equal(hv(sweep(A, 2, c(1, 2, 4))), base::sweep(m, 2, c(1, 2, 4)), ignore_attr = TRUE)
  expect_equal(hv(am_sweep(A, 2, c(1, 2, 4))), base::sweep(m, 2, c(1, 2, 4)), ignore_attr = TRUE)
  expect_s4_class(sweep(A, 2, c(1, 2, 4)), "adgeMatrix")
  expect_s4_class(sweep(S, 2, c(1, 2, 4), "*"), "adgCMatrix")
  expect_s4_class(colscale(S, c(1, 2, 4)), "adgCMatrix")
  expect_s4_class(rowscale(S, c(2, 4)), "adgCMatrix")
  expect_equal(hv(colscale(S, c(1, 2, 4))), t(t(m) * c(1, 2, 4)), ignore_attr = TRUE)
  expect_equal(hv(rowscale(S, c(2, 4))), m * c(2, 4), ignore_attr = TRUE)
  # Non-amatrix inputs are untouched.
  expect_identical(sweep(m, 2, c(1, 2, 4)), base::sweep(m, 2, c(1, 2, 4)))
})

test_that("mean, which and isSymmetric work on amatrix objects", {
  m <- matrix(c(1, 0, 3, 0, 5, 6), 2, 3)
  A <- adgeMatrix(m)
  S <- adgCMatrix(Matrix::Matrix(m, sparse = TRUE))
  expect_identical(mean(A), mean(m))
  expect_equal(mean(S), mean(m))
  expect_equal(mean(A > 2), mean(m > 2))
  expect_identical(as.integer(which(A > 2)), which(m > 2))
  expect_identical(as.integer(which(S > 2)), which(m > 2))
  sq <- matrix(c(2, 1, 1, 3), 2)
  expect_true(isSymmetric(adgeMatrix(sq)))
  expect_true(isSymmetric(adgCMatrix(Matrix::Matrix(sq, sparse = TRUE))))
  expect_false(isSymmetric(adgeMatrix(matrix(1:4 + 0, 2))))
})

test_that("mean() and which() work with only amatrix attached", {
  skip_if_not_installed("callr")
  out <- callr::r(function() {
    library(amatrix)
    A <- adgeMatrix(matrix(c(1, 0, 3, 0, 5, 6), 2, 3))
    list(
      matrix_attached = "package:Matrix" %in% search(),
      mean = mean(A),
      which = as.integer(which(A > 2)),
      sym = isSymmetric(adgeMatrix(diag(2)))
    )
  })
  expect_false(out$matrix_attached)
  expect_identical(out$mean, 2.5)
  expect_identical(out$which, c(3L, 5L, 6L))
  expect_true(out$sym)
})

test_that("dim<- reshapes to a valid object with a new identity", {
  m <- matrix(c(1, 0, 3, 0, 5, 6), 2, 3, dimnames = list(c("a", "b"), NULL))
  ref <- m
  dim(ref) <- c(3L, 2L)
  A <- adgeMatrix(m)
  d <- A
  dim(d) <- c(3L, 2L)
  expect_true(methods::validObject(d))
  expect_equal(as.matrix(d), ref, ignore_attr = TRUE)
  expect_null(rownames(d))
  expect_false(identical(d@object_id, A@object_id))

  S <- adgCMatrix(Matrix::Matrix(m, sparse = TRUE))
  ds <- S
  dim(ds) <- c(3L, 2L)
  expect_s4_class(ds, "adgCMatrix")
  expect_equal(as.matrix(ds), ref, ignore_attr = TRUE)
  expect_error({ z <- A; dim(z) <- c(4L, 2L) }, "do not match")
})
