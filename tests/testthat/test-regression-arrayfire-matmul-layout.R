# Regression repro metadata
# Seed: 20260423
# Dimensions: 12x5 %*% 5x7 non-square dense matrices
# Backend / precision / dispatch: arrayfire, fast precision, cold matmul dispatch
# R version / platform: captured by CI sessionInfo() on failure
# Issues: amatrix-x3c.2

test_that("ArrayFire matmul preserves R column-major matrix layout", {
  skip_on_cran()
  spec <- optional_backend_specs()[[match(
    "arrayfire",
    vapply(optional_backend_specs(), `[[`, character(1), "backend")
  )]]
  skip_if_backend_package_missing(spec)

  ns <- optional_backend_namespace(spec$package)
  enable_probe <- get("amatrix_arrayfire_enable_probe", envir = ns, inherits = FALSE)
  register_backend <- get(spec$register_fun, envir = ns, inherits = FALSE)
  native_available <- try(enable_probe(register = FALSE), silent = TRUE)
  skip_if_not(isTRUE(native_available), "arrayfire native backend not available")

  old_available <- getOption("amatrix.arrayfire.available")
  old_threshold <- getOption("amatrix.arrayfire.matmul_min_dim")
  options(amatrix.arrayfire.available = TRUE, amatrix.arrayfire.matmul_min_dim = 1L)
  on.exit({
    options(
      amatrix.arrayfire.available = old_available,
      amatrix.arrayfire.matmul_min_dim = old_threshold
    )
  }, add = TRUE)

  register_backend(overwrite = TRUE)

  set.seed(20260423L)
  x_host <- matrix(rnorm(12L * 5L), nrow = 12L, ncol = 5L)
  y_host <- matrix(rnorm(5L * 7L), nrow = 5L, ncol = 7L)
  x <- adgeMatrix(x_host, preferred_backend = "arrayfire", precision = "fast")
  y <- adgeMatrix(y_host, preferred_backend = "arrayfire", precision = "fast")

  plan <- amatrix_backend_plan(x, "matmul", y = y)
  expect_identical(plan$chosen, "arrayfire")
  expect_equal(as.matrix(x %*% y), x_host %*% y_host, tolerance = 1e-4)
})
