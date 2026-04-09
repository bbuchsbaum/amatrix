.opencl_eigen_spec <- function() {
  specs <- optional_backend_specs()
  specs[[match("opencl", vapply(specs, `[[`, character(1), "backend"))]]
}

.opencl_register_backend <- function(spec) {
  ns <- optional_backend_namespace(spec$package)
  get(spec$register_fun, envir = ns, inherits = FALSE)
}

test_that("OpenCL symmetric dense eigen path is chosen and matches base eigen", {
  spec <- .opencl_eigen_spec()
  skip_if_backend_package_missing(spec)

  register_backend <- .opencl_register_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    set.seed(20260409L)
    Z <- matrix(rnorm(16 * 16), nrow = 16, ncol = 16)
    S <- crossprod(Z) + diag(16)
    X <- adgeMatrix(S, preferred_backend = "opencl", precision = "fast")

    expect_identical(amatrix_backend_plan(X, "eigen")$chosen, "opencl")

    fit <- eigen(X, symmetric = TRUE)
    ref <- base::eigen(S, symmetric = TRUE)

    expect_equal(fit$values, ref$values, tolerance = 5e-6)
    resid <- norm(S %*% fit$vectors - fit$vectors %*% diag(fit$values, nrow = length(fit$values)), "F") / norm(S, "F")
    expect_lt(resid, 5e-5)
  })
})

test_that("OpenCL nonsymmetric dense eigen keeps host-general semantics", {
  spec <- .opencl_eigen_spec()
  skip_if_backend_package_missing(spec)

  register_backend <- .opencl_register_backend(spec)

  with_optional_backend_available(spec, {
    register_backend(overwrite = TRUE)

    A <- matrix(c(2, 1, 0, 3), nrow = 2)
    X <- adgeMatrix(A, preferred_backend = "opencl", precision = "fast")

    fit <- eigen(X, symmetric = FALSE)
    ref <- base::eigen(A, symmetric = FALSE)

    expect_equal(sort(Re(fit$values)), sort(Re(ref$values)), tolerance = 1e-10)
  })
})
