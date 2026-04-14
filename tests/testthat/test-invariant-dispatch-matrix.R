# Invariant-driven dispatch sweep.
#
# Purpose:
# - generated mixed-class `Ops` checks instead of one-off bug repros
# - generated `[<-` checks for missing-i / missing-j signatures
# - assert numeric correctness, class preservation, metadata preservation,
#   and dimname preservation against host references

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(amatrix))

.dispatch_binary_cases <- function(fixtures) {
  list(
    list(name = "dense_am_matrix", lhs = fixtures$dense$am, rhs = fixtures$dense$host),
    list(name = "dense_am_dge", lhs = fixtures$dense$am, rhs = fixtures$dense$dge),
    list(name = "dense_matrix_am", lhs = fixtures$dense$host, rhs = fixtures$dense$am),
    list(name = "dense_dge_am", lhs = fixtures$dense$dge, rhs = fixtures$dense$am),
    list(name = "sparse_am_sparse", lhs = fixtures$sparse$am, rhs = fixtures$sparse$dgC),
    list(name = "sparse_am_dense", lhs = fixtures$sparse$am, rhs = fixtures$sparse$matrix),
    list(name = "sparse_sparse_am", lhs = fixtures$sparse$dgC, rhs = fixtures$sparse$am),
    list(name = "sparse_dense_am", lhs = fixtures$sparse$matrix, rhs = fixtures$sparse$am)
  )
}

test_that("generated mixed Ops sweep preserves template semantics", {
  fixtures <- list(
    dense = .invariant_dense_fixture(backend = "cpu", policy = "opencl", precision = "strict"),
    sparse = .invariant_sparse_fixture(backend = "cpu", policy = "opencl", precision = "strict")
  )
  ops <- c("+", "-", "*")

  for (case in .dispatch_binary_cases(fixtures)) {
    for (op in ops) {
      info <- sprintf("case=%s op=%s", case$name, op)
      host_reference <- do.call(op, list(amatrix_materialize_host(case$lhs), amatrix_materialize_host(case$rhs)))
      result <- do.call(op, list(case$lhs, case$rhs))
      template <- .invariant_template(case$lhs, case$rhs)
      .invariant_expect_wrapped_result(result, template, host_reference, tolerance = 1e-12, info = info)
    }
  }
})

test_that("generated replacement sweep preserves class, metadata, and values", {
  dense <- .invariant_dense_fixture(backend = "cpu", policy = "opencl", precision = "strict")$am
  sparse <- .invariant_sparse_fixture(backend = "cpu", policy = "opencl", precision = "strict")$am

  cases <- list(
    list(
      name = "dense_row_replace",
      template = dense,
      host = amatrix_materialize_host(dense),
      apply = function(x) {
        x[1L, ] <- c(10, 11, 12, 13)
        x
      }
    ),
    list(
      name = "dense_col_replace",
      template = dense,
      host = amatrix_materialize_host(dense),
      apply = function(x) {
        x[, 2L] <- c(-4, -5, -6)
        x
      }
    ),
    list(
      name = "sparse_row_replace",
      template = sparse,
      host = amatrix_materialize_host(sparse),
      apply = function(x) {
        x[2L, ] <- c(7, 0, 8, 0)
        x
      }
    ),
    list(
      name = "sparse_col_replace",
      template = sparse,
      host = amatrix_materialize_host(sparse),
      apply = function(x) {
        x[, 3L] <- c(0, 9, 10)
        x
      }
    )
  )

  for (case in cases) {
    info <- case$name
    host_reference <- case$apply(case$host)
    result <- case$apply(case$template)
    .invariant_expect_wrapped_result(result, case$template, host_reference, tolerance = 1e-12, info = info)
  }
})

test_that("generated coercion round-trips preserve metadata across class grid", {
  dense <- .invariant_dense_fixture(backend = "cpu", policy = "opencl", precision = "strict")$am
  sparse <- .invariant_sparse_fixture(backend = "cpu", policy = "opencl", precision = "strict")$am

  roundtrips <- list(
    list(
      name = "dense_via_dge",
      source = dense,
      mid_class = "dgeMatrix",
      out_class = "adgeMatrix"
    ),
    list(
      name = "sparse_via_dgC",
      source = sparse,
      mid_class = "dgCMatrix",
      out_class = "adgCMatrix"
    ),
    list(
      name = "sparse_to_dense",
      source = sparse,
      mid_class = "dgeMatrix",
      out_class = "adgeMatrix"
    ),
    list(
      name = "dense_to_sparse",
      source = dense,
      mid_class = "dgCMatrix",
      out_class = "adgCMatrix"
    )
  )

  for (case in roundtrips) {
    info <- case$name
    mid <- methods::as(case$source, case$mid_class)
    out <- methods::as(mid, case$out_class)
    expect_equal(as.matrix(out), as.matrix(case$source), tolerance = 1e-12, info = info)
    .invariant_expect_template_metadata(out, case$source, info = info)
    expect_equal(base::dimnames(out), base::dimnames(case$source), info = info)
  }
})
