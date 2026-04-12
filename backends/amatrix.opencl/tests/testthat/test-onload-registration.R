test_that("opencl onLoad registration no longer depends on enable option", {
  backend_state <- amatrix:::.amatrix_state$backends
  onload <- get(".onLoad", envir = asNamespace("amatrix.opencl"), inherits = FALSE)
  old_options <- options(amatrix.enable_opencl = FALSE)
  had_backend <- exists("opencl", envir = backend_state, inherits = FALSE)
  saved_backend <- if (had_backend) get("opencl", envir = backend_state, inherits = FALSE) else NULL

  on.exit({
    options(old_options)
    if (exists("opencl", envir = backend_state, inherits = FALSE)) {
      rm("opencl", envir = backend_state)
    }
    if (had_backend) {
      assign("opencl", saved_backend, envir = backend_state)
    }
  }, add = TRUE)

  if (exists("opencl", envir = backend_state, inherits = FALSE)) {
    rm("opencl", envir = backend_state)
  }

  expect_false("opencl" %in% ls(envir = backend_state, all.names = FALSE))
  expect_no_error(onload(NULL, "amatrix.opencl"))
  expect_true("opencl" %in% ls(envir = backend_state, all.names = FALSE))
})
