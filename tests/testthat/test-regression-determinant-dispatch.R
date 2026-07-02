# Regression repro metadata
# Seed: 20260702
# Dimensions: 3x3 dense double input with a positive diagonal shift
# Backend / precision / dispatch: cpu, strict, fresh-process base::det() path
# R version / platform: captured by child sessionInfo() on failure
# Issue: bd-01KWHDRQV6BX162NSXQP5DWAAC

test_that("base det() dispatch works on adgeMatrix without Matrix attached [bd-01KWHDRQV6BX162NSXQP5DWAAC]", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("pkgload")
  repo_dir <- .amatrix_source_tree_dir()
  skip_if(is.null(repo_dir), "source tree not reachable (installed-pkg context)")

  result <- callr::r(
    function(repo_dir) {
      pkgload::load_all(repo_dir, quiet = TRUE)
      set.seed(20260702L)
      host <- matrix(runif(9L) + 1, 3L, 3L) + diag(3L) * 5
      x <- adgeMatrix(host, preferred_backend = "cpu")

      list(
        matrix_attached = "package:Matrix" %in% search(),
        value = det(x),
        expected = det(host),
        session = capture.output(sessionInfo())
      )
    },
    args = list(repo_dir = repo_dir)
  )

  expect_false(result$matrix_attached, info = paste(result$session, collapse = "\n"))
  expect_equal(result$value, result$expected, tolerance = 1e-10)
})
