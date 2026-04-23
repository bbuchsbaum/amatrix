# Regression repro metadata
# Seed: 20260423
# Dimensions: 3x3 dense positive-definite-like input
# Backend / precision / dispatch: cpu, strict, fresh-process attach path
# R version / platform: captured by child sessionInfo() on failure
# Issue: amatrix-1ha

.matrix_generic_repo_dir <- function() {
  candidates <- unique(c(
    tryCatch(getNamespaceInfo(asNamespace("amatrix"), "path"), error = function(e) NULL),
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", "..")
  ))
  candidates <- Filter(Negate(is.null), candidates)
  matches <- candidates[file.exists(file.path(candidates, "DESCRIPTION"))]
  if (length(matches) == 0L) {
    return(NULL)
  }
  normalizePath(matches[[1L]], winslash = "/", mustWork = TRUE)
}

test_that("Matrix-style generics work after plain amatrix attach without Matrix attached [amatrix-1ha]", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("pkgload")
  repo_dir <- .matrix_generic_repo_dir()
  skip_if(is.null(repo_dir), "source tree not reachable (installed-pkg context)")

  result <- callr::r(
    function(repo_dir) {
      pkgload::load_all(repo_dir, quiet = TRUE)
      set.seed(20260423L)
      x <- adgeMatrix(matrix(runif(9L) + 1, 3L, 3L) + diag(3L) * 5, backend = "cpu")

      list(
        matrix_attached = "package:Matrix" %in% search(),
        rowSums_symbol = find("rowSums")[[1L]],
        colSums_symbol = find("colSums")[[1L]],
        t_class = class(t(x)),
        chol_class = class(chol(x)),
        rowSums_value = rowSums(x),
        colSums_value = colSums(x),
        diag_value = diag(x),
        solve_class = class(solve(x)),
        session = capture.output(sessionInfo())
      )
    },
    args = list(repo_dir = repo_dir)
  )

  expect_false(result$matrix_attached, info = paste(result$session, collapse = "\n"))
  expect_identical(result$rowSums_symbol, "package:amatrix")
  expect_identical(result$colSums_symbol, "package:amatrix")
  expect_s3_class(structure(list(), class = result$t_class), "aTransposeView")
  expect_identical(result$chol_class[[1L]], "adgeMatrix")
  expect_type(result$rowSums_value, "double")
  expect_length(result$rowSums_value, 3L)
  expect_type(result$colSums_value, "double")
  expect_length(result$colSums_value, 3L)
  expect_type(result$diag_value, "double")
  expect_length(result$diag_value, 3L)
  expect_identical(result$solve_class[[1L]], "adgeMatrix")
})
