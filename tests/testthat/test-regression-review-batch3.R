# Regressions for the third batch of pre-CRAN review fixes.

test_that("block_lanczos converges by default and warns when it cannot", {
  set.seed(7)
  G <- matrix(rnorm(300 * 120), 300)
  sg <- svd(G, nu = 0, nv = 0)$d
  set.seed(1)
  res <- block_lanczos(G, nv = 4)
  expect_true(res$converged)
  expect_lt(max(abs(res$d - sg[1:4]) / sg[1:4]), 1e-6)
  expect_length(res$residuals, length(res$d))
  expect_warning(block_lanczos(G, nv = 4, n_steps = 1L), class = "amatrix_not_converged")
  expect_error(block_lanczos(matrix(c(1, NA, 3, 4), 2), nv = 1), "NA")
})

test_that("svd_factor subspace is scale invariant and tolerates rank deficiency", {
  set.seed(2)
  W <- matrix(rnorm(2000), 100)
  sw <- svd(W, nu = 0, nv = 0)$d
  for (scale in c(1, 1e-9, 1e-12, 1e9)) {
    f <- svd_factor(adgeMatrix(W * scale), k = 10, method = "subspace")
    expect_lt(max(abs(f@d / scale - sw[1:10]) / sw[1:10]), 1e-8)
  }
  low <- matrix(rnorm(60 * 3), 60) %*% matrix(rnorm(3 * 30), 3)
  f3 <- svd_factor(adgeMatrix(low), k = 5, method = "subspace")
  expect_length(f3@d, 5L)
  # Null singular values must not blow up the coefficients.
  Y <- matrix(rnorm(120), 60)
  expect_lt(max(abs(pca_coef(f3, Y))), 1e3)
})

test_that("normal-equation fits signal rank deficiency and ill-conditioning", {
  set.seed(3)
  y <- rnorm(50)
  Xd <- cbind(1, matrix(rnorm(150), 50))
  expect_error(lm_fit(cbind(Xd, Xd[, 2]), y), class = "amatrix_rank_deficient")
  expect_error(wls_fit(cbind(Xd, Xd[, 2]), y, runif(50)), class = "amatrix_rank_deficient")
  expect_no_error(ridge_fit(cbind(Xd, Xd[, 2]), y, 1))

  U <- qr.Q(qr(matrix(rnorm(200), 50)))
  V <- qr.Q(qr(matrix(rnorm(16), 4)))
  ill <- U %*% diag(c(1, 1, 1, 1e-7)) %*% t(V)
  well <- U %*% diag(c(1, 1, 1, 1e-2)) %*% t(V)
  expect_warning(lm_fit(ill, y), class = "amatrix_ill_conditioned")
  expect_no_warning(lm_fit(well, y))
  expect_no_warning(lm_fit(ill, y, method = "qr"))
})

test_that("irlba wrapper accepts sparse input and a user fastpath argument", {
  skip_if_not_installed("irlba")
  set.seed(4)
  As <- Matrix::rsparsematrix(200, 100, 0.1)
  ref <- svd(as.matrix(As), nu = 0, nv = 0)$d[1:3]
  expect_equal(amatrix::irlba(As, nv = 3)$d, ref, tolerance = 1e-6)
  dense <- as.matrix(As)
  expect_equal(amatrix::irlba(dense, nv = 3, fastpath = TRUE)$d, ref, tolerance = 1e-6)
  expect_equal(amatrix::irlba(adgeMatrix(dense), nv = 3)$d, ref, tolerance = 1e-6)
})

test_that("rowsum and segment reductions handle names and missing groups", {
  m <- matrix(1:12 + 0, 4, dimnames = list(NULL, c("a", "b", "c")))
  A <- adgeMatrix(m)
  S <- adgCMatrix(Matrix::Matrix(m, sparse = TRUE))
  g <- c("x", "y", "x", "y")
  expect_identical(dimnames(rowsum(A, g)), dimnames(rowsum(m, g)))
  expect_equal(as.matrix(rowsum(A, g)), rowsum(m, g))

  gn <- c("x", NA, "x", "y")
  ref <- suppressWarnings(rowsum(m, gn))
  expect_equal(as.matrix(suppressWarnings(rowsum(A, gn))), ref)
  expect_warning(got <- rowsum(S, gn), "missing values")
  expect_equal(unname(as.matrix(got)), unname(ref))

  labels <- c(1L, NA, 1L, 2L)
  expect_equal(unname(as.matrix(segment_sum(m, labels, 2L))),
               unname(rbind(m[1, ] + m[3, ], m[4, ])))
  expect_equal(unname(as.matrix(segment_mean(m, labels, 2L))),
               unname(rbind((m[1, ] + m[3, ]) / 2, m[4, ])))
})

test_that("assorted wrappers validate inputs and match base", {
  Q <- matrix(c(2, 1, 1, 3), 2, dimnames = list(NULL, c("p", "q")))
  expect_identical(names(solve(adgeMatrix(Q), c(1, 2))), names(solve(Q, c(1, 2))))

  set.seed(5)
  X <- Matrix::rsparsematrix(20, 5, 0.5)
  y <- rnorm(20)
  qs <- am_qr(X)
  expect_type(qr_info(qs), "list")
  expect_equal(as.numeric(qr.solve(qs, y)),
               as.numeric(qr.coef(qr(as.matrix(X)), y)), tolerance = 1e-8)

  expect_error(kron_matrix(diag(2), diag(3)) %*% 1:5, "non-conformable")
  m <- matrix(1:12 + 0, 4)
  expect_error(pairwise_sqdist_argmin(m, t(m[1:2, ]), x_norms = 1:2), "x_norms")
  expect_error(sinkhorn(matrix(runif(9), 3), max_iter = 1e10), "positive integer")

  A0 <- diag(3) * 2
  f <- chol_factor(adgeMatrix(A0))
  U <- matrix(c(1, 0, 0), 3)
  expect_equal(woodbury_logdet(f, U, V = t(U), C_inv = matrix(1)),
               as.numeric(determinant(A0 + U %*% t(U))$modulus))
  expect_warning(neg <- woodbury_logdet(f, U, V = t(U), C_inv = matrix(-1 / 5)), "negative")
  expect_true(is.nan(neg))

  Xs <- Matrix::Matrix(cbind(rnorm(200, 1e8, 1), rnorm(200)), sparse = TRUE)
  expect_warning(covariance(Xs), class = "amatrix_precision_loss")
})

test_that("large dense products and scalar extraction give base results", {
  set.seed(6)
  n <- 520L   # above the in-place product threshold (n * n > 250000)
  M <- matrix(rnorm(n * n), n, dimnames = list(NULL, paste0("c", seq_len(n))))
  A <- adgeMatrix(M)
  B <- adgeMatrix(M)
  prod_ab <- A %*% B
  expect_s4_class(prod_ab, "adgeMatrix")
  expect_equal(as.matrix(prod_ab), M %*% M)
  expect_identical(dimnames(prod_ab), dimnames(M %*% M))
  cp <- crossprod(A)
  expect_equal(as.matrix(cp), crossprod(M))
  expect_true(isSymmetric(as.matrix(cp)))
  expect_equal(as.matrix(tcrossprod(A, B)), tcrossprod(M, M), ignore_attr = TRUE)
  expect_identical(A[[3, 4]], M[[3, 4]])
  expect_identical(A[[7]], M[[7]])
  expect_error(A[[n + 1L, 1L]])
})

test_that("amatrix_calibrate does not write to disk unless asked", {
  expect_false(isTRUE(formals(amatrix_calibrate)$persist))
})
