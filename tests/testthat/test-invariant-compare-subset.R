# Invariant checks for compare and subset semantics.
#
# Focus:
# - comparisons stay logical, match host references, and preserve the
#   amatrix logical wrapper classes (adlgeMatrix/adlgCMatrix) with their
#   backend metadata — the amatrix-47w contract (compare must NOT demote
#   to plain lgeMatrix/lgCMatrix, which loses residency)
# - subset extraction matches host semantics for drop=TRUE/FALSE
# - matrix-like subset results preserve amatrix metadata when still matrix-like

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(amatrix))

test_that("compare invariants match host semantics and keep amatrix logical classes", {
  dense <- .invariant_dense_fixture(backend = "cpu", policy = "opencl", precision = "strict")
  sparse <- .invariant_sparse_fixture(backend = "cpu", policy = "opencl", precision = "strict")
  ops <- c("==", "!=", "<", ">", "<=", ">=")

  cases <- list(
    list(name = "dense_scalar", lhs = dense$am, rhs = 0),
    list(name = "dense_matrix", lhs = dense$am, rhs = dense$host),
    list(name = "sparse_scalar", lhs = sparse$am, rhs = 0),
    list(name = "sparse_matrix", lhs = sparse$am, rhs = sparse$matrix)
  )

  for (case in cases) {
    for (op in ops) {
      info <- sprintf("case=%s op=%s", case$name, op)
      result <- do.call(op, list(case$lhs, case$rhs))
      reference <- do.call(op, list(amatrix_materialize_host(case$lhs), case$rhs))

      expect_true(inherits(result, "aMatrix"), info = info)
      expect_true(inherits(result, "adlgeMatrix") || inherits(result, "adlgCMatrix"),
                  info = info)
      expect_identical(result@preferred_backend, case$lhs@preferred_backend, info = info)
      expect_identical(as.matrix(result), as.matrix(reference), info = info)
    }
  }
})

test_that("subset invariants match host semantics across drop modes", {
  dense <- .invariant_dense_fixture(backend = "cpu", policy = "opencl", precision = "strict")
  sparse <- .invariant_sparse_fixture(backend = "cpu", policy = "opencl", precision = "strict")

  cases <- list(
    list(name = "dense", host = dense$host, x = dense$am),
    list(name = "sparse", host = sparse$host, x = sparse$am)
  )
  specs <- list(
    list(name = "row_drop", eval = function(obj) obj[1L, ]),
    list(name = "col_drop", eval = function(obj) obj[, 2L]),
    list(name = "row_keep", eval = function(obj) obj[1L, , drop = FALSE]),
    list(name = "col_keep", eval = function(obj) obj[, 2L, drop = FALSE]),
    list(name = "block_default", eval = function(obj) obj[, 1:2]),
    list(name = "block_keep", eval = function(obj) obj[, 1:2, drop = FALSE])
  )

  for (case in cases) {
    for (spec in specs) {
      info <- sprintf("case=%s subset=%s", case$name, spec$name)
      result <- spec$eval(case$x)
      reference <- spec$eval(case$host)

      if (is.matrix(reference) || inherits(reference, "Matrix")) {
        expect_true(inherits(result, "aMatrix"), info = info)
        expect_equal(as.matrix(result), as.matrix(reference), tolerance = 0, info = info)
        expect_identical(result@preferred_backend, case$x@preferred_backend, info = info)
        expect_identical(result@policy, case$x@policy, info = info)
        expect_identical(result@precision, case$x@precision, info = info)
      } else {
        expect_equal(result, reference, tolerance = 0, info = info)
      }
    }
  }
})
