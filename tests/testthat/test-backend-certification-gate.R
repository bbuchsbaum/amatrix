test_that("backend certification gate exposes the expected backend matrix", {
  skip_if_not(file.exists(file.path("..", "..", "tools", "check-backend-certification.R")),
              "backend certification tool not reachable")

  old_autorun <- Sys.getenv("AMATRIX_BACKEND_CERTIFICATION_AUTORUN", unset = NA_character_)
  Sys.setenv(AMATRIX_BACKEND_CERTIFICATION_AUTORUN = "0")
  on.exit({
    if (is.na(old_autorun)) {
      Sys.unsetenv("AMATRIX_BACKEND_CERTIFICATION_AUTORUN")
    } else {
      Sys.setenv(AMATRIX_BACKEND_CERTIFICATION_AUTORUN = old_autorun)
    }
  }, add = TRUE)

  env <- new.env(parent = globalenv())
  sys.source(file.path("..", "..", "tools", "check-backend-certification.R"), envir = env)

  gates <- env$backend_certification_gates()
  expect_named(gates, c("mlx", "arrayfire", "opencl", "metal"))
  expect_match(gates$mlx$filter, "backend-certification-mlx", fixed = TRUE)
  expect_match(gates$arrayfire$filter, "arrayfire-matmul-layout", fixed = TRUE)
  expect_identical(gates$opencl$env, c(AMATRIX_OPENCL_PROBE_GPU = "1"))
  expect_identical(gates$metal$env, c(AMATRIX_METAL_PROBE_GPU = "1"))
})

test_that("backend certification gate argument parser validates backend names", {
  skip_if_not(file.exists(file.path("..", "..", "tools", "check-backend-certification.R")),
              "backend certification tool not reachable")

  old_autorun <- Sys.getenv("AMATRIX_BACKEND_CERTIFICATION_AUTORUN", unset = NA_character_)
  Sys.setenv(AMATRIX_BACKEND_CERTIFICATION_AUTORUN = "0")
  on.exit({
    if (is.na(old_autorun)) {
      Sys.unsetenv("AMATRIX_BACKEND_CERTIFICATION_AUTORUN")
    } else {
      Sys.setenv(AMATRIX_BACKEND_CERTIFICATION_AUTORUN = old_autorun)
    }
  }, add = TRUE)

  env <- new.env(parent = globalenv())
  sys.source(file.path("..", "..", "tools", "check-backend-certification.R"), envir = env)

  parsed <- env$parse_backend_certification_args(c("--backends=mlx,metal", "--allow-skips"))
  expect_identical(parsed$backends, c("mlx", "metal"))
  expect_true(parsed$allow_skips)
  expect_error(
    env$parse_backend_certification_args("--backends=bogus"),
    "unknown backend gate"
  )
})
