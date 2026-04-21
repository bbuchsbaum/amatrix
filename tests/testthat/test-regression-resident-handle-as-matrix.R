# Regression repro metadata
# Seed: none (deterministic 2 x 2 dense fixture)
# Dimensions: 2 x 2 dense matrix with row/column names
# Backend / precision / dispatch: mock recording backend / strict / fresh-process
#   source-tree pkgload::load_all() resident_handle materialization
# R version / platform: captured by CI sessionInfo() on failure
# Link: amatrix-ze2

test_that("resident_handle materializes through generic as.matrix in a fresh load_all process [regression]", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("pkgload")
  skip_if_not(file.exists("DESCRIPTION"),
              "source tree not reachable (installed-pkg context)")

  result <- callr::r(
    function(path) {
      pkgload::load_all(path, quiet = TRUE)

      store <- new.env(parent = emptyenv())
      backend <- list(
        capabilities = function() c("matmul", "crossprod", "tcrossprod", "ewise", "rowSums", "colSums"),
        features = function() character(0),
        precision_modes = function() "strict",
        available = function() TRUE,
        supports = function(op, x, y = NULL) TRUE,
        matmul = function(x, y) x %*% y,
        crossprod = function(x, y = NULL, ...) if (is.null(y)) base::crossprod(x) else base::crossprod(x, y),
        tcrossprod = function(x, y = NULL, ...) if (is.null(y)) base::tcrossprod(x) else base::tcrossprod(x, y),
        ewise = function(x, lhs, rhs = NULL, op, ...) if (is.null(rhs)) do.call(op, list(lhs)) else do.call(op, list(lhs, rhs)),
        rowSums = function(x, na.rm = FALSE, dims = 1L) base::rowSums(x, na.rm = na.rm),
        colSums = function(x, na.rm = FALSE, dims = 1L) base::colSums(x, na.rm = na.rm),
        resident_store = function(key, mat) assign(key, mat, envir = store),
        resident_has = function(key) exists(key, envir = store, inherits = FALSE),
        resident_drop = function(key) if (exists(key, envir = store, inherits = FALSE)) rm(list = key, envir = store),
        resident_materialize = function(key) get(key, envir = store, inherits = FALSE)
      )
      amatrix_register_backend("resident_handle_dispatch", backend, overwrite = TRUE)
      expected <- matrix(c(1, 2, 3, 4), nrow = 2, dimnames = list(c("r1", "r2"), c("c1", "c2")))
      capture <- function(expr) {
        tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "error_string"))
      }
      h <- resident_handle(expected, backend = "resident_handle_dispatch")
      list(
        direct = amatrix:::as.matrix.resident_handle(h),
        generic = capture(as.matrix(h)),
        generic_base = capture(base::as.matrix(h))
      )
    },
    args = list(getwd()),
    spinner = FALSE,
    show = FALSE
  )

  expect_equal(result$direct, matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("r1", "r2"), c("c1", "c2"))
  ))
  expect_equal(result$generic, result$direct)
  expect_equal(result$generic_base, result$direct)
})
