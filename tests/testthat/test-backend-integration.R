for (spec in optional_backend_specs()) {
  test_that(sprintf("core harness integrates with optional backend package %s", spec$package), {
    skip_if_backend_package_missing(spec)
    # This block executes real backend compute (matmul/qr/solve/... below), so a
    # loadable package is not sufficient: the native runtime must be functional.
    # arrayfire is the one backend whose package builds and loads from source
    # even when its native runtime (libaf) is absent (mlx/opencl fail to build in
    # that case and are already skipped above), so the ops would error at run
    # time. Probe the native backend and skip cleanly if it is unavailable,
    # exactly as test-regression-arrayfire-matmul-layout.R does.
    if (identical(spec$backend, "arrayfire")) {
      af_ns <- optional_backend_namespace(spec$package)
      enable_probe <- get("amatrix_arrayfire_enable_probe", envir = af_ns, inherits = FALSE)
      skip_if_not(
        isTRUE(tryCatch(enable_probe(register = FALSE), error = function(e) FALSE)),
        "arrayfire native backend not available"
      )
    }

    ns <- optional_backend_namespace(spec$package)
    register_backend <- get(spec$register_fun, envir = ns, inherits = FALSE)

    with_optional_backend_available(spec, {
      capabilities <- backend_package_capabilities(spec)
      backend_available <- backend_package_available(spec)
      register_backend(overwrite = TRUE)
      status <- amatrix_backend_status(spec$backend)

      expect_true(spec$backend %in% amatrix_backend_names())
      expect_identical(amatrix_backend_capabilities(spec$backend), capabilities)
      expect_identical(status$name, spec$backend)
      expect_identical(status$available, backend_available)
      expect_identical(status$precision_modes, "fast")
      expect_true(nzchar(status$features))
      expect_identical(status$capabilities, paste(capabilities, collapse = ","))

      if (identical(spec$backend, "arrayfire")) {
        old_af_backend <- amatrix.arrayfire:::amatrix_arrayfire_active_backend()
        amatrix.arrayfire:::amatrix_arrayfire_set_backend("cpu")
        on.exit(
          amatrix.arrayfire:::amatrix_arrayfire_set_backend(if (identical(old_af_backend, 4L)) "opencl" else "cpu"),
          add = TRUE
        )
      }

      x_host <- matrix(c(1, 2, 3, 4), nrow = 2)
      dense <- adgeMatrix(x_host, preferred_backend = spec$backend)
      dense_fast <- adgeMatrix(x_host, preferred_backend = spec$backend, precision = "fast")
      sparse <- adgCMatrix(matrix(c(1, 0, 0, 1), nrow = 2), preferred_backend = spec$backend)

      dense_plan <- amatrix_backend_plan(dense, "matmul", y = diag(2))
      dense_fast_plan <- amatrix_backend_plan(dense_fast, "matmul", y = diag(2))
      qr_fast_plan <- amatrix_backend_plan(dense_fast, "qr")
      unsupported_plan <- amatrix_backend_plan(dense, "solve")
      sparse_plan <- amatrix_backend_plan(sparse, "matmul", y = diag(2))
      summary <- amatrix_backend_matrix(
        dense_fast,
        ops = c("matmul", "qr", "solve"),
        y_map = list(matmul = diag(2))
      )
      fac <- qr(dense_fast)
      fac_base <- base::qr(x_host)
      backend_impl <- amatrix:::.amatrix_get_backend(spec$backend)
      matmul_expected <- if (isTRUE(backend_impl$supports("matmul", dense_fast, y = diag(2)))) spec$backend else "cpu"
      qr_expected <- if (isTRUE(backend_impl$supports("qr", dense_fast))) spec$backend else "cpu"

      expect_identical(dense_plan$chosen, "cpu")
      expect_identical(dense_fast_plan$chosen, matmul_expected)
      expect_identical(qr_fast_plan$chosen, qr_expected)
      expect_identical(unsupported_plan$chosen, "cpu")
      expect_identical(sparse_plan$chosen, "cpu")
      solve_expected <- if (isTRUE(backend_impl$supports("solve", dense_fast))) spec$backend else "cpu"
      expect_identical(summary$chosen, c(matmul_expected, qr_expected, solve_expected))
      expect_identical(summary$cpu_fallback, c(identical(matmul_expected, "cpu"), identical(qr_expected, "cpu"), identical(solve_expected, "cpu")))

      expect_s4_class(dense_fast %*% diag(2), "adgeMatrix")
      expect_s3_class(fac, "amQR")
      expect_equal(as.matrix(dense_fast %*% diag(2)), x_host %*% diag(2), tolerance = 1e-10)
      expect_equal(as.matrix(crossprod(dense_fast)), crossprod(x_host), tolerance = 1e-10)
      expect_equal(as.matrix(qr.Q(fac)), base::qr.Q(fac_base), tolerance = 1e-4)
      expect_equal(as.matrix(qr.R(fac)), base::qr.R(fac_base), tolerance = 1e-4)
      expect_equal(as.matrix(qr.solve(fac, x_host)), base::qr.solve(fac_base, x_host), tolerance = 1e-4)
      expect_equal(as.matrix(qr.coef(fac, x_host)), base::qr.coef(fac_base, x_host), tolerance = 1e-4)
      expect_equal(rowSums(dense_fast), rowSums(x_host), tolerance = 1e-10)
      expect_equal(colSums(dense_fast), colSums(x_host), tolerance = 1e-10)
      expect_equal(as.matrix(solve(dense)), solve(x_host), tolerance = 1e-10)
      expect_equal(as.matrix(sparse %*% diag(2)), as.matrix(as(sparse, "dgCMatrix") %*% diag(2)), tolerance = 1e-10)
    })
  })

  test_that(sprintf("optional backend %s auto-registers on demand", spec$backend), {
    skip_if_backend_package_missing(spec)

    ns <- optional_backend_namespace(spec$package)
    register_backend <- get(spec$register_fun, envir = ns, inherits = FALSE)

    had_backend <- exists(spec$backend, envir = amatrix:::.amatrix_state$backends, inherits = FALSE)
    saved_backend <- if (had_backend) {
      get(spec$backend, envir = amatrix:::.amatrix_state$backends, inherits = FALSE)
    } else {
      NULL
    }

    if (had_backend) {
      rm(list = spec$backend, envir = amatrix:::.amatrix_state$backends)
    }

    on.exit({
      if (exists(spec$backend, envir = amatrix:::.amatrix_state$backends, inherits = FALSE)) {
        rm(list = spec$backend, envir = amatrix:::.amatrix_state$backends)
      }
      if (!is.null(saved_backend)) {
        amatrix_register_backend(spec$backend, saved_backend, overwrite = TRUE)
      } else {
        register_backend(overwrite = TRUE)
      }
    }, add = TRUE)

    expect_false(exists(spec$backend, envir = amatrix:::.amatrix_state$backends, inherits = FALSE))

    backend <- amatrix:::.amatrix_get_backend(spec$backend)
    status <- amatrix_backend_status(spec$backend)

    expect_true(is.list(backend))
    expect_true(exists(spec$backend, envir = amatrix:::.amatrix_state$backends, inherits = FALSE))
    expect_identical(status$name, spec$backend)
  })

  test_that(sprintf("optional backend %s respects explicit per-backend disable", spec$backend), {
    skip_if_backend_package_missing(spec)

    skip_if(is.null(spec$disable_option), "backend has no per-backend disable option")

    old_disable <- getOption(spec$disable_option)
    options(structure(list(TRUE), names = spec$disable_option))
    on.exit(options(structure(list(old_disable), names = spec$disable_option)), add = TRUE)

    had_backend <- exists(spec$backend, envir = amatrix:::.amatrix_state$backends, inherits = FALSE)
    saved_backend <- if (had_backend) {
      get(spec$backend, envir = amatrix:::.amatrix_state$backends, inherits = FALSE)
    } else {
      NULL
    }

    if (had_backend) {
      rm(list = spec$backend, envir = amatrix:::.amatrix_state$backends)
    }

    on.exit({
      if (exists(spec$backend, envir = amatrix:::.amatrix_state$backends, inherits = FALSE)) {
        rm(list = spec$backend, envir = amatrix:::.amatrix_state$backends)
      }
      if (!is.null(saved_backend)) {
        amatrix_register_backend(spec$backend, saved_backend, overwrite = TRUE)
      }
    }, add = TRUE)

    expect_false(amatrix:::.amatrix_try_register_optional_backend(spec$backend))
    expect_false(spec$backend %in% amatrix_backend_names())
    expect_error(amatrix:::.amatrix_get_backend(spec$backend), "not registered")
  })

  test_that(sprintf("optional backend %s respects disabled global auto-registration", spec$backend), {
    skip_if_backend_package_missing(spec)

    old_opt <- getOption("amatrix.optional_backends")
    options(amatrix.optional_backends = FALSE)
    on.exit(options(amatrix.optional_backends = old_opt), add = TRUE)

    had_backend <- exists(spec$backend, envir = amatrix:::.amatrix_state$backends, inherits = FALSE)
    saved_backend <- if (had_backend) {
      get(spec$backend, envir = amatrix:::.amatrix_state$backends, inherits = FALSE)
    } else {
      NULL
    }

    if (had_backend) {
      rm(list = spec$backend, envir = amatrix:::.amatrix_state$backends)
    }

    on.exit({
      if (exists(spec$backend, envir = amatrix:::.amatrix_state$backends, inherits = FALSE)) {
        rm(list = spec$backend, envir = amatrix:::.amatrix_state$backends)
      }
      if (!is.null(saved_backend)) {
        amatrix_register_backend(spec$backend, saved_backend, overwrite = TRUE)
      }
    }, add = TRUE)

    expect_false(amatrix:::.amatrix_try_register_optional_backend(spec$backend))
    expect_false(spec$backend %in% amatrix_backend_names())
    expect_error(amatrix:::.amatrix_get_backend(spec$backend), "not registered")
  })
}
