# Invariant-driven wrapper contract checks.
#
# Focus:
# - reductions with na.rm must match base / Matrix semantics
# - rowsum grouping semantics must match base
# - distance / kernel wrappers must reject non-finite inputs consistently

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(amatrix))

test_that("reduction wrappers match host references across matrix classes", {
  dense <- .invariant_na_dense_fixture(backend = "cpu", policy = "opencl", precision = "strict")
  sparse_host <- methods::as(dense$host, "dgCMatrix")
  sparse <- as_adgCMatrix(sparse_host, preferred_backend = "cpu", policy = "opencl", precision = "strict")

  cases <- list(
    list(name = "dense", host = dense$host, x = dense$am),
    list(name = "sparse", host = sparse_host, x = sparse)
  )
  reducers <- list(
    list(name = "rowSums", fn = rowSums, ref = function(x, na.rm) Matrix::rowSums(x, na.rm = na.rm)),
    list(name = "colSums", fn = colSums, ref = function(x, na.rm) Matrix::colSums(x, na.rm = na.rm)),
    list(name = "rowMeans", fn = rowmeans, ref = function(x, na.rm) Matrix::rowMeans(x, na.rm = na.rm)),
    list(name = "colMeans", fn = colmeans, ref = function(x, na.rm) Matrix::colMeans(x, na.rm = na.rm))
  )

  for (case in cases) {
    for (reducer in reducers) {
      for (na_rm in c(FALSE, TRUE)) {
        info <- sprintf("case=%s reducer=%s na.rm=%s", case$name, reducer$name, na_rm)
        expect_equal(
          reducer$fn(case$x, na.rm = na_rm),
          reducer$ref(case$host, na.rm = na_rm),
          tolerance = 1e-12,
          info = info
        )
      }
    }
  }
})

test_that("rowsum wrappers match base semantics across reorder and na.rm", {
  host <- matrix(
    c(1, 2, NA, 4,
      5, NA, 7, 8,
      9, 10, 11, NA),
    nrow = 4L,
    ncol = 3L
  )
  groups <- c("b", "a", "b", "a")
  dense <- as_adgeMatrix(host, preferred_backend = "cpu", policy = "opencl", precision = "strict")
  sparse <- as_adgCMatrix(methods::as(host, "dgCMatrix"), preferred_backend = "cpu", policy = "opencl", precision = "strict")

  cases <- list(
    list(name = "dense", fn = rowsum.adgeMatrix, x = dense),
    list(name = "sparse", fn = rowsum.adgCMatrix, x = sparse)
  )

  for (case in cases) {
    for (reorder in c(FALSE, TRUE)) {
      for (na_rm in c(FALSE, TRUE)) {
        info <- sprintf("case=%s reorder=%s na.rm=%s", case$name, reorder, na_rm)
        expect_equal(
          case$fn(case$x, groups, reorder = reorder, na.rm = na_rm),
          base::rowsum(host, groups, reorder = reorder, na.rm = na_rm),
          tolerance = 1e-12,
          info = info
        )
      }
    }
  }
})

test_that("distance and kernel wrappers reject non-finite contamination by type", {
  bad_inputs <- list(
    list(name = "NA", value = NA_real_),
    list(name = "Inf", value = Inf),
    list(name = "-Inf", value = -Inf),
    list(name = "NaN", value = NaN)
  )

  for (bad in bad_inputs) {
    x <- matrix(rnorm(12L), nrow = 3L, ncol = 4L)
    x[2L, 3L] <- bad$value
    info <- paste("bad=", bad$name)

    expect_error(dist_matrix(x), regexp = NULL, info = info)
    expect_error(kernel_matrix(x, kernel = "rbf", sigma = 1), regexp = NULL, info = info)
  }
})
