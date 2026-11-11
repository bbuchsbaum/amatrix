# Drop-in semantics regressions from the pre-CRAN review. Each block compares
# an amatrix operation against the same operation on the plain base / Matrix
# object.

host_value <- function(x) if (methods::is(x, "Matrix")) as.matrix(x) else x

expect_same_as_host <- function(got, ref, info = NULL) {
  got <- host_value(got)
  ref <- host_value(ref)
  expect_identical(dim(got), dim(ref), info = info)
  expect_equal(got, ref, ignore_attr = TRUE, info = info)
}

test_that("one-subscript indexing is linear, not row indexing", {
  m <- matrix(c(1, 0, 3, 0, 5, 6), 2, 3)
  sp <- Matrix::Matrix(m, sparse = TRUE)
  for (pair in list(list(adgeMatrix(m), m), list(adgCMatrix(sp), sp))) {
    X <- pair[[1L]]
    R <- pair[[2L]]
    expect_same_as_host(X[2], R[2])
    expect_same_as_host(X[c(1, 4)], R[c(1, 4)])
    expect_same_as_host(X[5], R[5])
    expect_same_as_host(X[-1], R[-1])
    expect_same_as_host(X[cbind(1:2, 2:3)], R[cbind(1:2, 2:3)])
    # Two-subscript forms are unchanged.
    expect_same_as_host(X[2, ], R[2, ])
    expect_same_as_host(X[, 2], R[, 2])
    expect_same_as_host(X[2, , drop = FALSE], R[2, , drop = FALSE])
    expect_same_as_host(X[1, 2], R[1, 2])

    assigners <- list(
      function(z) { z[2] <- 99; z },
      function(z) { z[c(1, 6)] <- c(7, 8); z },
      function(z) { z[2, ] <- 9; z },
      function(z) { z[, 3] <- 4; z },
      function(z) { z[cbind(1:2, 2:3)] <- -1; z }
    )
    for (f in assigners) {
      expect_same_as_host(f(X), f(R), info = deparse(body(f))[2L])
    }
  }
  A <- adgeMatrix(m)
  expect_same_as_host(t(A)[2], t(m)[2])
  expect_equal(rev(A), rev(m))
  expect_equal(median(A), median(m))
})

test_that("log() honours the base argument", {
  p <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3)
  P <- adgeMatrix(p)
  SP <- adgCMatrix(Matrix::Matrix(p, sparse = TRUE))
  expect_same_as_host(log(P, 2), log(p, 2))
  expect_same_as_host(log(P, base = 10), log10(p))
  expect_same_as_host(log(P), log(p))
  expect_same_as_host(log(SP, 2), log(p, 2))
  expect_same_as_host(exp(P), exp(p))
})

test_that("norm() matches base for every type and leaves plain matrices alone", {
  m <- matrix(c(1, -2, 3, 4, 5, -6), 2, 3)
  A <- adgeMatrix(m)
  S <- adgCMatrix(Matrix::Matrix(m, sparse = TRUE))

  for (type in c("O", "o", "1", "I", "i", "F", "f", "E", "M", "m", "2")) {
    expect_equal(norm(A, type), base::norm(m, type), info = type)
    if (!identical(type, "2")) {
      expect_equal(norm(S, type), base::norm(m, type), info = type)
    }
  }
  # Missing type defaults to the one norm, as in base.
  expect_equal(norm(A), base::norm(m))
  expect_equal(norm(S), base::norm(m))
  expect_error(norm(A, "zz"))

  # Plain matrices are untouched by amatrix being loaded.
  expect_equal(norm(m), base::norm(m))
  expect_equal(norm(m, "O"), base::norm(m, "O"))
  expect_equal(norm(m, "E"), base::norm(m, "E"))
  expect_equal(Matrix::norm(m), base::norm(m))

  # The spectral norm is exact, not a randomized estimate.
  set.seed(1)
  X <- matrix(rnorm(200 * 150), 200)
  expect_equal(norm(adgeMatrix(X), "2"), base::norm(X, "2"), tolerance = 1e-10)
  expect_equal(norm(X, "2"), base::norm(X, "2"))
})

test_that("svd() still works for non-amatrix inputs", {
  m <- matrix(c(1, -2, 3, 4, 5, -6), 2, 3)
  expect_equal(svd(1:3)$d, base::svd(1:3)$d)
  expect_equal(svd(m)$d, base::svd(m)$d)
  expect_equal(svd(Matrix::Matrix(m))$d, base::svd(m)$d)
  expect_equal(svd(adgeMatrix(m))$d, base::svd(m)$d)
})

test_that("dense elementwise ops survive sparse or diagonal Matrix results", {
  m <- matrix(c(1, 0, 3, 0, 5, 6), 2, 3)
  A <- adgeMatrix(m)
  B <- adgeMatrix(m + 1)
  S <- adgCMatrix(Matrix::Matrix(m, sparse = TRUE))
  zeros <- matrix(c(0, 0, 0, 0, 0, 1), 2, 3)

  expect_s4_class(A * (m > 4), "adgeMatrix")
  expect_same_as_host(A * (m > 4), m * (m > 4))
  expect_same_as_host(A * zeros, m * zeros)
  expect_same_as_host(adgeMatrix(diag(2) + 1) * diag(2), (diag(2) + 1) * diag(2))
  expect_same_as_host(A == matrix(0, 2, 3), m == matrix(0, 2, 3))
  expect_same_as_host(A & B, m & (m + 1))
  expect_same_as_host(A * S, m * m)
  expect_same_as_host(A == S, m == m)
  expect_same_as_host(adgeMatrix(matrix(2)) * matrix(3), matrix(6))
})

test_that("sparse results are rewrapped without allocating a dense copy", {
  n <- 60000L
  i <- seq(1L, n, length.out = 3000L)
  big <- adgCMatrix(Matrix::sparseMatrix(
    i = as.integer(i), j = as.integer(rev(i)), x = seq_along(i) / 10, dims = c(n, n)
  ))
  # Each of these would need ~27 GB if the result were densified.
  expect_s4_class(big * 2, "adgCMatrix")
  expect_s4_class(big + big, "adgCMatrix")
  expect_s4_class(abs(big), "adgCMatrix")
  expect_s4_class(big %*% big, "adgCMatrix")
  expect_s4_class(cbind(big, big), "adgCMatrix")
  expect_s4_class(big > 1, "adlgCMatrix")
  expect_length(diag(big), n)
  z <- big
  z[1, 1] <- 5
  expect_equal(z[1, 1], 5)
})

test_that("cbind/rbind work on logical amatrix objects", {
  m <- matrix(c(1, 0, 3, 0, 5, 6), 2, 3)
  sp <- Matrix::Matrix(m, sparse = TRUE)
  A <- adgeMatrix(m)
  L <- A > 2
  LS <- adgCMatrix(sp) > 2
  expect_same_as_host(cbind(L, L), cbind(m > 2, m > 2))
  expect_same_as_host(rbind(L, A > 3), rbind(m > 2, m > 3))
  expect_same_as_host(cbind(A, L), cbind(m, m > 2))
  expect_same_as_host(cbind(LS, LS), cbind(sp > 2, sp > 2))
})

test_that("numeric %*% adgeMatrix has base R's shape", {
  m <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3, dimnames = list(NULL, c("a", "b", "c")))
  A <- adgeMatrix(m)
  got <- 1:2 %*% A
  expect_s4_class(got, "adgeMatrix")
  expect_same_as_host(got, 1:2 %*% m)
  expect_identical(colnames(got), c("a", "b", "c"))

  one_row <- matrix(c(1, 2, 3), 1)
  expect_same_as_host(1:4 %*% adgeMatrix(one_row), 1:4 %*% one_row)
  expect_error(1:5 %*% A, "non-conformable")
})

test_that("matrix functions are correct for 1 x 1 input and keep dimnames", {
  expect_equal(as.numeric(as.matrix(mat_sqrt(matrix(2.25)))), 1.5)
  expect_equal(as.numeric(as.matrix(mat_sqrt(matrix(4)))), 2)
  expect_equal(as.numeric(as.matrix(mat_log(matrix(exp(2))))), 2)

  S <- crossprod(matrix(sin(1:20), 5)) + diag(4)
  dimnames(S) <- list(letters[1:4], letters[1:4])
  root <- mat_sqrt(S)
  expect_equal(as.matrix(root) %*% as.matrix(root), S, ignore_attr = TRUE)
  expect_identical(dimnames(root), dimnames(S))
  expect_same_as_host(mat_pow(adgeMatrix(S), 2), S %*% S)
})

test_that("wrappers accept plain matrices as documented", {
  m <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3)
  m2 <- matrix(1:9 + 0, 3)
  A <- adgeMatrix(m)
  S <- crossprod(matrix(sin(1:20), 5)) + diag(4)

  expect_s4_class(matmul(m, m2), "adgeMatrix")
  expect_same_as_host(matmul(m, m2), m %*% m2)
  expect_same_as_host(matmul(Matrix::Matrix(m, sparse = TRUE), m2), m %*% m2)
  expect_same_as_host(gemm(m, m2), m %*% m2)
  expect_same_as_host(gemm(A, adgeMatrix(m2), C = m, beta = 2), m %*% m2 + 2 * m)
  expect_equal(rowsums(m), rowSums(m))
  expect_equal(colsums(m), colSums(m))
  expect_equal(sort(eigh(S)$values), sort(eigen(S)$values))
  expect_same_as_host(ewise("+", m, m), m + m)
})

test_that("diag() extraction matches base without copying the matrix", {
  S <- crossprod(matrix(sin(1:20), 5)) + diag(4)
  dimnames(S) <- list(letters[1:4], letters[1:4])
  expect_identical(diag(adgeMatrix(S)), diag(S))
  expect_equal(diag(adgCMatrix(Matrix::Matrix(S, sparse = TRUE))), diag(S))
  expect_equal(diag(adgeMatrix(S[, 1:2])), diag(S[, 1:2]))
})

test_that("t() views behave like ordinary matrices", {
  m <- matrix(c(1, 2, 3, 4, 5, 7), 2, 3)
  A <- adgeMatrix(m)
  tA <- t(A)
  tm <- t(m)
  expect_identical(length(tA), length(tm))
  expect_equal(sum(tA), sum(tm))
  expect_equal(range(tA), range(tm))
  expect_equal(rowSums(tA), rowSums(tm))
  expect_equal(colMeans(tA), colMeans(tm))
  expect_equal(as.vector(tA), as.vector(tm))
  expect_equal(diag(tA), diag(tm))
  expect_same_as_host(-tA, -tm)
  expect_same_as_host(exp(tA), exp(tm))
  expect_same_as_host(crossprod(tA), crossprod(tm))
  expect_same_as_host(tcrossprod(tA), tcrossprod(tm))
  expect_same_as_host(1:3 %*% tA, 1:3 %*% tm)
  expect_equal(norm(tA, "F"), norm(tm, "F"))
  sq <- matrix(c(2, 1, 1, 3), 2)
  expect_same_as_host(solve(t(adgeMatrix(sq))), solve(t(sq)))
  expect_true(isSymmetric(adgeMatrix(sq)))
  expect_false(isSymmetric(adgeMatrix(matrix(1:4 + 0, 2))))
})

test_that("results keep the amatrix class where Matrix methods used to win dispatch", {
  m <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3)
  A <- adgeMatrix(m)
  S <- adgCMatrix(Matrix::Matrix(m, sparse = TRUE))
  dg <- Matrix::Matrix(m)

  z <- A; z[1, ] <- 0
  expect_s4_class(z, "adgeMatrix")
  z <- A; z[, 1] <- c(1, 2)
  expect_s4_class(z, "adgeMatrix")
  z <- A; z[2] <- 0
  expect_s4_class(z, "adgeMatrix")
  z <- A; z[A > 2] <- 0
  expect_s4_class(z, "adgeMatrix")
  expect_same_as_host(z, {r <- m; r[m > 2] <- 0; r})
  z <- S; z[S > 2] <- 0
  expect_s4_class(z, "adgCMatrix")

  expect_s4_class(A %*% 1:3, "adgeMatrix")
  expect_same_as_host(A %*% 1:3, m %*% 1:3)
  expect_same_as_host(S %*% 1:3, m %*% 1:3)
  expect_true(inherits(S %*% 1:3, "aMatrix"))
  expect_s4_class(crossprod(A, m), "adgeMatrix")
  expect_s4_class(crossprod(A, dg), "adgeMatrix")
  expect_s4_class(tcrossprod(A, m), "adgeMatrix")
  expect_same_as_host(crossprod(A, 1:2), crossprod(m, 1:2))
  expect_s4_class(kronecker(A, dg), "adgeMatrix")
  expect_same_as_host(kronecker(A, dg), kronecker(m, m))
  expect_s4_class((S > 2) & (S > 0), "adlgCMatrix")
})
