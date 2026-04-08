test_that("random dense arithmetic and products match host behavior", {
  set.seed(1)

  for (iter in seq_len(10)) {
    case <- random_dense_case()
    x_host <- case$x
    y_host <- case$y
    v_host <- case$rhs

    x <- adgeMatrix(x_host)
    y <- adgeMatrix(y_host)

    expect_equal(as.matrix(x + y), x_host + y_host, tolerance = 1e-10)
    expect_equal(as.matrix(x - 2), x_host - 2, tolerance = 1e-10)
    expect_equal(as.matrix(x * y), x_host * y_host, tolerance = 1e-10)
    expect_equal(as.matrix(x %*% v_host), x_host %*% v_host, tolerance = 1e-10)
    expect_equal(as.matrix(crossprod(x)), crossprod(x_host), tolerance = 1e-10)
    expect_equal(as.matrix(tcrossprod(x)), tcrossprod(x_host), tolerance = 1e-10)
    expect_equal(rowSums(x), rowSums(x_host), tolerance = 1e-10)
    expect_equal(colSums(x), colSums(x_host), tolerance = 1e-10)
  }
})

test_that("random sparse operations stay Matrix-backed and numerically aligned", {
  set.seed(2)

  for (iter in seq_len(10)) {
    case <- random_sparse_case()
    x <- adgCMatrix(case$raw)
    host <- case$host
    rhs <- case$rhs

    expect_equal(as.matrix(x %*% rhs), as.matrix(host %*% rhs), tolerance = 1e-10)
    expect_equal(rowSums(x), Matrix::rowSums(host), tolerance = 1e-10)
    expect_equal(colSums(x), Matrix::colSums(host), tolerance = 1e-10)

    y <- x
    y[1, 1] <- 7
    host[1, 1] <- 7
    expect_true(inherits(y, "aMatrix"))
    expect_equal(as.matrix(y), as.matrix(host), tolerance = 1e-10)
  }
})

test_that("random dense solve and chol match host behavior on SPD inputs", {
  set.seed(3)

  for (iter in seq_len(10)) {
    n <- sample(2:6, 1)
    z <- matrix(rnorm(n * n), nrow = n)
    spd <- crossprod(z) + diag(n) * 0.5
    rhs <- matrix(rnorm(n * sample(1:4, 1)), nrow = n)

    x <- adgeMatrix(spd)

    expect_equal(as.matrix(chol(x)), base::chol(spd), tolerance = 1e-10)
    expect_equal(as.matrix(solve(x)), base::solve(spd), tolerance = 1e-10)
    expect_equal(as.matrix(solve(x, rhs)), base::solve(spd, rhs), tolerance = 1e-10)
  }
})

test_that("wrapper path honors a registered non-cpu backend when supported", {
  counter <- new.env(parent = emptyenv())
  counter$matmul <- 0L
  counter$ewise <- 0L
  counter$rowSums <- 0L
  counter$colSums <- 0L
  counter$crossprod <- 0L
  counter$tcrossprod <- 0L
  counter$solve <- 0L
  counter$chol <- 0L
  counter$qr <- 0L
  counter$svd <- 0L
  counter$eigen <- 0L
  counter$diag <- 0L

  with_registered_backend("recording", make_recording_backend(counter), {
    x_host <- matrix(1:4, nrow = 2)
    x <- adgeMatrix(x_host, preferred_backend = "recording")

    expect_s4_class(x %*% diag(2), "adgeMatrix")
    expect_s4_class(x + 1, "adgeMatrix")
    expect_equal(rowSums(x), rowSums(x_host))
    expect_equal(colSums(x), colSums(x_host))

    expect_equal(counter$matmul, 1L)
    expect_equal(counter$ewise, 1L)
    expect_equal(counter$rowSums, 1L)
    expect_equal(counter$colSums, 1L)
    expect_equal(counter$crossprod, 0L)
  })
})

test_that("unsupported operations fall back to cpu semantics without hitting the backend", {
  counter <- new.env(parent = emptyenv())
  counter$matmul <- 0L
  counter$ewise <- 0L
  counter$rowSums <- 0L
  counter$colSums <- 0L
  counter$crossprod <- 0L
  counter$tcrossprod <- 0L
  counter$solve <- 0L
  counter$chol <- 0L
  counter$qr <- 0L
  counter$svd <- 0L
  counter$eigen <- 0L
  counter$diag <- 0L

  with_registered_backend("recording_fallback", make_recording_backend(counter, supported_ops = c("matmul")), {
    x_host <- matrix(c(4, 1, 1, 3), nrow = 2)
    rhs <- matrix(c(1, 2), nrow = 2)
    x <- adgeMatrix(x_host, preferred_backend = "recording_fallback")

    expect_equal(as.matrix(x %*% rhs), x_host %*% rhs)
    expect_equal(as.matrix(crossprod(x)), crossprod(x_host))
    expect_equal(as.matrix(chol(x)), chol(x_host))
    expect_equal(as.matrix(solve(x)), solve(x_host))

    expect_equal(counter$matmul, 1L)
    expect_equal(counter$crossprod, 0L)
    expect_equal(counter$chol, 0L)
    expect_equal(counter$solve, 0L)
  })
})

test_that("backend plan explains chosen and skipped candidates", {
  counter <- new.env(parent = emptyenv())
  counter$matmul <- 0L
  counter$ewise <- 0L
  counter$rowSums <- 0L
  counter$colSums <- 0L
  counter$crossprod <- 0L
  counter$tcrossprod <- 0L
  counter$solve <- 0L
  counter$chol <- 0L
  counter$qr <- 0L
  counter$svd <- 0L
  counter$eigen <- 0L
  counter$diag <- 0L

  with_registered_backend("recording_plan", make_recording_backend(counter, supported_ops = c("matmul")), {
    x <- adgeMatrix(matrix(1:4, nrow = 2), preferred_backend = "recording_plan")

    matmul_plan <- amatrix_backend_plan(x, "matmul", y = diag(2))
    solve_plan <- amatrix_backend_plan(x, "solve")

    expect_identical(matmul_plan$chosen, "recording_plan")
    expect_identical(solve_plan$chosen, "cpu")

    matmul_candidates <- setNames(matmul_plan$candidates, vapply(matmul_plan$candidates, `[[`, character(1), "name"))
    solve_candidates <- setNames(solve_plan$candidates, vapply(solve_plan$candidates, `[[`, character(1), "name"))

    expect_true(matmul_candidates$recording_plan$registered)
    expect_true(matmul_candidates$recording_plan$available)
    expect_true(matmul_candidates$recording_plan$supported)
    expect_true(matmul_candidates$recording_plan$chosen)

    expect_true(solve_candidates$recording_plan$registered)
    expect_true(solve_candidates$recording_plan$available)
    expect_false(solve_candidates$recording_plan$supported)
    expect_false(solve_candidates$recording_plan$chosen)
    expect_true(solve_candidates$cpu$chosen)
  })
})

test_that("backend matrix summarizes chosen backends across operations", {
  counter <- new.env(parent = emptyenv())
  counter$matmul <- 0L
  counter$ewise <- 0L
  counter$rowSums <- 0L
  counter$colSums <- 0L
  counter$crossprod <- 0L
  counter$tcrossprod <- 0L
  counter$solve <- 0L
  counter$chol <- 0L
  counter$qr <- 0L
  counter$svd <- 0L
  counter$eigen <- 0L
  counter$diag <- 0L

  with_registered_backend("recording_matrix", make_recording_backend(counter, supported_ops = c("matmul", "ewise")), {
    x <- adgeMatrix(matrix(1:4, nrow = 2), preferred_backend = "recording_matrix")
    summary <- amatrix_backend_matrix(
      x,
      ops = c("matmul", "ewise", "solve"),
      y_map = list(matmul = diag(2), ewise = 1)
    )

    expect_identical(summary$op, c("matmul", "ewise", "solve"))
    expect_identical(summary$precision, c("strict", "strict", "strict"))
    expect_identical(summary$chosen, c("recording_matrix", "recording_matrix", "cpu"))
    expect_identical(summary$chosen_path, c("cold", "cold", "cold"))
    expect_identical(summary$resident_reuse, c(FALSE, FALSE, FALSE))
    expect_identical(summary$cpu_fallback, c(FALSE, FALSE, TRUE))
    expect_true(grepl("recording_matrix\\[RAP-C-KSX\\]", summary$candidate_summary[[1]]))
    expect_true(grepl("cpu\\[RAP-C-KSX\\]", summary$candidate_summary[[3]]))
  })
})

test_that("backend capabilities and status are explicit", {
  counter <- new.env(parent = emptyenv())
  counter$matmul <- 0L
  counter$ewise <- 0L
  counter$rowSums <- 0L
  counter$colSums <- 0L
  counter$crossprod <- 0L
  counter$tcrossprod <- 0L
  counter$solve <- 0L
  counter$chol <- 0L
  counter$qr <- 0L
  counter$svd <- 0L
  counter$eigen <- 0L
  counter$diag <- 0L

  with_registered_backend("recording_status", make_recording_backend(counter, supported_ops = c("matmul", "ewise")), {
    expect_true("cpu" %in% amatrix_backend_names())
    expect_identical(
      amatrix_backend_capabilities("recording_status"),
      c("ewise", "matmul")
    )
    expect_identical(amatrix_backend_features("recording_status"), "dense_f64")

    status <- amatrix_backend_status(c("cpu", "recording_status"))
    expect_identical(status$name, c("cpu", "recording_status"))
    expect_identical(status$available, c(TRUE, TRUE))
    expect_identical(status$precision_modes, c("strict,fast", "strict,fast"))
    expect_true(grepl("dense_f64", status$features[[1]]))
    expect_identical(status$features[[2]], "dense_f64")
    expect_identical(status$residency_capable, c(FALSE, TRUE))
    expect_true(grepl("matmul", status$capabilities[[1]]))
    expect_identical(status$capabilities[[2]], "ewise,matmul")
  })
})

test_that("strict precision falls back to cpu when backend is fast-only", {
  counter <- new.env(parent = emptyenv())
  counter$matmul <- 0L

  with_registered_backend(
    "fast_only",
    make_recording_backend(counter, supported_ops = c("matmul"), precision_modes = "fast"),
    {
      x <- adgeMatrix(matrix(1:4, nrow = 2), preferred_backend = "fast_only", precision = "strict")
      plan <- amatrix_backend_plan(x, "matmul", y = diag(2))

      expect_identical(plan$requested_precision, "strict")
      expect_identical(plan$chosen, "cpu")

      candidates <- setNames(plan$candidates, vapply(plan$candidates, `[[`, character(1), "name"))
      expect_false(candidates$fast_only$precision_compatible)
      expect_false(candidates$fast_only$supported)

      expect_equal(as.matrix(x %*% diag(2)), matrix(1:4, nrow = 2))
      expect_equal(counter$matmul, 0L)
    }
  )
})

test_that("fast precision enables fast-only backend routing", {
  counter <- new.env(parent = emptyenv())
  counter$matmul <- 0L

  with_registered_backend(
    "fast_only_enabled",
    make_recording_backend(counter, supported_ops = c("matmul"), precision_modes = "fast"),
    {
      x <- adgeMatrix(matrix(1:4, nrow = 2), preferred_backend = "fast_only_enabled", precision = "fast")
      plan <- amatrix_backend_plan(x, "matmul", y = diag(2))

      expect_identical(plan$requested_precision, "fast")
      expect_identical(plan$chosen, "fast_only_enabled")
      expect_equal(counter$matmul, 0L)

      expect_equal(as.matrix(x %*% diag(2)), matrix(1:4, nrow = 2))
      expect_equal(counter$matmul, 1L)
    }
  )
})

test_that("backend selection respects preferred backend before policy fallback", {
  preferred_counter <- new.env(parent = emptyenv())
  policy_counter <- new.env(parent = emptyenv())

  for (counter in list(preferred_counter, policy_counter)) {
    counter$matmul <- 0L
    counter$ewise <- 0L
    counter$rowSums <- 0L
    counter$colSums <- 0L
    counter$crossprod <- 0L
    counter$tcrossprod <- 0L
    counter$solve <- 0L
    counter$chol <- 0L
    counter$qr <- 0L
    counter$svd <- 0L
    counter$eigen <- 0L
    counter$diag <- 0L
  }

  with_registered_backend("mlx", make_recording_backend(policy_counter, supported_ops = c("matmul")), {
    with_registered_backend("arrayfire", make_recording_backend(preferred_counter, supported_ops = c("matmul")), {
      x <- adgeMatrix(
        matrix(1:4, nrow = 2),
        preferred_backend = "arrayfire",
        policy = "mlx"
      )

      plan <- amatrix_backend_plan(x, "matmul", y = diag(2))
      result <- x %*% diag(2)

      expect_s4_class(result, "adgeMatrix")
      expect_identical(plan$chosen, "arrayfire")
      expect_equal(preferred_counter$matmul, 1L)
      expect_equal(policy_counter$matmul, 0L)
    })
  })
})

test_that("dense chaining can reuse resident backend state", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend("resident_recording", make_recording_backend(counter, supported_ops = c("matmul", "ewise")), {
    x <- adgeMatrix(matrix(1:4, nrow = 2), preferred_backend = "resident_recording")
    y <- diag(2)

    first <- x %*% y
    second <- first * 2
    third <- second %*% y

    expect_s4_class(first, "adgeMatrix")
    expect_s4_class(second, "adgeMatrix")
    expect_s4_class(third, "adgeMatrix")
    expect_equal(as.matrix(third), ((matrix(1:4, nrow = 2) %*% y) * 2) %*% y)

    expect_true(counter$resident_store >= 2L)
    expect_true(counter$matmul_resident >= 2L)
    expect_true(counter$ewise_resident >= 1L)

    info <- amatrix_residency_info(second)
    expect_identical(info$backend[[1]], "resident_recording")
    expect_true(info$live[[1]])
  })
})

test_that("resident dense crossprod chaining can reuse backend state", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend("resident_products", make_recording_backend(counter, supported_ops = c("matmul", "crossprod", "tcrossprod", "ewise")), {
    x_host <- matrix(1:9, nrow = 3)
    x <- adgeMatrix(x_host, preferred_backend = "resident_products")

    cp <- crossprod(x)
    tcp <- tcrossprod(x)
    mixed <- cp * 2

    expect_s4_class(cp, "adgeMatrix")
    expect_s4_class(tcp, "adgeMatrix")
    expect_s4_class(mixed, "adgeMatrix")
    expect_equal(as.matrix(cp), crossprod(x_host))
    expect_equal(as.matrix(tcp), tcrossprod(x_host))
    expect_equal(as.matrix(mixed), crossprod(x_host) * 2)

    expect_true(counter$resident_store >= 1L)
    expect_true(counter$crossprod_resident >= 1L)
    expect_true(counter$tcrossprod_resident >= 1L)
    expect_true(counter$ewise_resident >= 1L)
  })
})

test_that("backend plan distinguishes cold support from resident reuse", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "resident_plan",
    make_recording_backend(
      counter,
      supported_ops = c("ewise"),
      cold_supported_ops = c("ewise"),
      resident_supported_ops = c("matmul", "ewise")
    ),
    {
    x <- adgeMatrix(matrix(1:4, nrow = 2), preferred_backend = "resident_plan")

    cold_plan <- amatrix_backend_plan(x, "matmul", y = diag(2))
    resident_x <- x * 2
    resident_plan <- amatrix_backend_plan(resident_x, "matmul", y = diag(2))
    summary <- amatrix_backend_matrix(resident_x, ops = c("matmul", "crossprod"), y_map = list(matmul = diag(2)))

    expect_identical(cold_plan$chosen, "cpu")
    expect_identical(cold_plan$chosen_path, "cold")
    expect_false(cold_plan$candidates[[1]]$supported_cold)
    expect_false(cold_plan$candidates[[1]]$supported_resident)
    expect_false(cold_plan$candidates[[1]]$resident_active)

    expect_identical(resident_plan$chosen, "resident_plan")
    expect_identical(resident_plan$chosen_path, "resident")
    expect_false(resident_plan$candidates[[1]]$supported_cold)
    expect_true(resident_plan$candidates[[1]]$supported_resident)
    expect_true(resident_plan$candidates[[1]]$resident_active)

    expect_identical(summary$chosen, c("resident_plan", "cpu"))
    expect_identical(summary$chosen_path, c("resident", "cold"))
    expect_identical(summary$resident_reuse, c(TRUE, FALSE))

    result <- resident_x %*% diag(2)
    expect_s4_class(result, "adgeMatrix")
    expect_equal(as.matrix(result), as.matrix(resident_x) %*% diag(2))
    expect_true(counter$matmul_resident >= 1L)
  })
})

test_that("execution info summarizes residency and operation planning", {
  counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "exec_info",
    make_recording_backend(
      counter,
      supported_ops = c("ewise"),
      cold_supported_ops = c("ewise"),
      resident_supported_ops = c("matmul", "ewise")
    ),
    {
    x <- adgeMatrix(matrix(1:4, nrow = 2), preferred_backend = "exec_info")
    resident_x <- x * 2
    info <- amatrix_execution_info(
      resident_x,
      ops = c("matmul", "crossprod"),
      y_map = list(matmul = diag(2))
    )

    expect_identical(info$object_id, resident_x@object_id)
    expect_identical(info$preferred_backend, "exec_info")
    expect_identical(info$pinned_backend, "exec_info")
    expect_identical(info$policy, resident_x@policy)
    expect_identical(info$precision, resident_x@precision)
    expect_true(is.data.frame(info$residency))
    expect_true(is.data.frame(info$plans))
    expect_identical(info$residency$backend[[1]], "exec_info")
    expect_identical(info$residency$pinned_backend[[1]], "exec_info")
    expect_true(info$residency$live[[1]])
    expect_identical(info$plans$chosen, c("exec_info", "cpu"))
    expect_identical(info$plans$pinned_backend, c("exec_info", "exec_info"))
    expect_identical(info$plans$resident_reuse, c(TRUE, FALSE))
  })
})

test_that("resident objects are pinned and do not hop to another accelerator backend", {
  resident_counter <- new.env(parent = emptyenv())
  other_counter <- new.env(parent = emptyenv())

  with_registered_backend(
    "pinned_backend",
    make_recording_backend(
      resident_counter,
      supported_ops = c("ewise"),
      cold_supported_ops = c("ewise"),
      resident_supported_ops = c("ewise")
    ),
    {
      with_registered_backend(
        "mlx",
        make_recording_backend(
          other_counter,
          supported_ops = c("matmul")
        ),
        {
          x <- adgeMatrix(
            matrix(1:4, nrow = 2),
            preferred_backend = "pinned_backend",
            policy = "mlx"
          )

          resident_x <- x * 2
          plan <- amatrix_backend_plan(resident_x, "matmul", y = diag(2))
          result <- resident_x %*% diag(2)

          expect_identical(plan$pinned_backend, "pinned_backend")
          expect_identical(plan$preferred, c("pinned_backend", "cpu"))
          expect_identical(plan$chosen, "cpu")
          expect_identical(plan$chosen_path, "cold")
          expect_equal(if (is.null(other_counter$matmul)) 0L else other_counter$matmul, 0L)
          expect_equal(as.matrix(result), as.matrix(resident_x) %*% diag(2))
        }
      )
    }
  )
})

test_that("lm_fit matches least-squares coefficients for shared-X many-Y", {
  set.seed(11)
  X_host <- cbind(1, matrix(rnorm(40), nrow = 10, ncol = 4))
  beta_host <- matrix(rnorm(15), nrow = 5, ncol = 3)
  Y_host <- X_host %*% beta_host + matrix(rnorm(30, sd = 1e-6), nrow = 10, ncol = 3)

  fit <- lm_fit(adgeMatrix(X_host), Y_host)
  coef_host <- solve(crossprod(X_host), crossprod(X_host, Y_host))

  expect_s3_class(fit, "lm_fit")
  expect_s4_class(fit$coefficients, "adgeMatrix")
  expect_equal(as.matrix(coef(fit)), coef_host, tolerance = 1e-8)
  expect_equal(as.matrix(fitted(fit)), X_host %*% coef_host, tolerance = 1e-8)
  expect_equal(as.matrix(residuals(fit)), Y_host - X_host %*% coef_host, tolerance = 1e-8)
  expect_identical(fit$precision, "strict")
})

test_that("lm_fit supports a qr-backed path", {
  set.seed(16)
  X_host <- cbind(1, matrix(rnorm(40), nrow = 10, ncol = 4))
  beta_host <- matrix(rnorm(15), nrow = 5, ncol = 3)
  Y_host <- X_host %*% beta_host + matrix(rnorm(30, sd = 1e-6), nrow = 10, ncol = 3)

  fit <- lm_fit(adgeMatrix(X_host), Y_host, method = "qr")
  coef_host <- base::qr.coef(base::qr(X_host), Y_host)

  expect_s3_class(fit, "lm_fit")
  expect_identical(fit$method, "qr")
  expect_null(fit$xtx)
  expect_null(fit$xty)
  expect_identical(fit$qr_representation, "base_qr")
  expect_identical(fit$qr_helper_path, "compact_factor")
  expect_true(fit$qr_compact_factor_available)
  expect_identical(fit$qr_compact_factor_source, "native")
  expect_equal(as.matrix(coef(fit)), coef_host, tolerance = 1e-8)
  expect_equal(as.matrix(fitted(fit)), base::qr.fitted(base::qr(X_host), Y_host), tolerance = 1e-8)
  expect_equal(as.matrix(residuals(fit)), base::qr.resid(base::qr(X_host), Y_host), tolerance = 1e-8)
})

test_that("lm_fit routes hot kernels through a preferred backend in fast mode", {
  counter <- new.env(parent = emptyenv())
  counter$crossprod <- 0L
  counter$matmul <- 0L
  counter$solve <- 0L

  with_registered_backend(
    "fit_backend",
    make_recording_backend(
      counter,
      supported_ops = c("crossprod", "matmul"),
      features = c("dense_f32")
    ),
    {
      set.seed(12)
      X <- adgeMatrix(matrix(rnorm(30), nrow = 10, ncol = 3), preferred_backend = "fit_backend", precision = "fast")
      Y <- matrix(rnorm(20), nrow = 10, ncol = 2)

      fit <- lm_fit(X, Y)

      expect_s3_class(fit, "lm_fit")
      expect_true(counter$crossprod >= 2L)
      expect_true(counter$matmul >= 1L)
      expect_identical(counter$solve, 0L)
    }
  )
})

test_that("lm_fit reuses cached shared-X work across repeated fits", {
  counter <- new.env(parent = emptyenv())
  counter$crossprod <- 0L
  counter$matmul <- 0L
  counter$solve <- 0L

  with_registered_backend(
    "fit_cache_backend",
    make_recording_backend(
      counter,
      supported_ops = c("crossprod"),
      features = c("dense_f32")
    ),
    {
      set.seed(13)
      X <- adgeMatrix(
        matrix(rnorm(36), nrow = 12, ncol = 3),
        preferred_backend = "fit_cache_backend",
        precision = "fast"
      )
      Y1 <- matrix(rnorm(24), nrow = 12, ncol = 2)
      Y2 <- matrix(rnorm(24), nrow = 12, ncol = 2)

      fit1 <- lm_fit(X, Y1, include_fitted = FALSE, include_residuals = FALSE)
      fit2 <- lm_fit(X, Y2, include_fitted = FALSE, include_residuals = FALSE)

      expect_false(fit1$cache_reused)
      expect_true(fit2$cache_reused)
      expect_identical(fit1$cache_key, fit2$cache_key)
      expect_equal(counter$crossprod, 3L)
      expect_equal(as.matrix(coef(fit1)), solve(crossprod(as.matrix(X)), crossprod(as.matrix(X), Y1)), tolerance = 1e-8)
      expect_equal(as.matrix(coef(fit2)), solve(crossprod(as.matrix(X)), crossprod(as.matrix(X), Y2)), tolerance = 1e-8)
    }
  )
})

test_that("lm_fit reuses cached qr work across repeated fits", {
  set.seed(17)
  X <- adgeMatrix(matrix(rnorm(48), nrow = 12, ncol = 4))
  Y1 <- matrix(rnorm(24), nrow = 12, ncol = 2)
  Y2 <- matrix(rnorm(24), nrow = 12, ncol = 2)

  fit1 <- lm_fit(X, Y1, include_fitted = FALSE, include_residuals = FALSE, method = "qr")
  fit2 <- lm_fit(X, Y2, include_fitted = FALSE, include_residuals = FALSE, method = "qr")

  expect_false(fit1$cache_reused)
  expect_true(fit2$cache_reused)
  expect_identical(fit1$cache_key, fit2$cache_key)
  expect_identical(fit1$method, "qr")
  expect_identical(fit2$method, "qr")
  expect_identical(fit1$qr_helper_path, "compact_factor")
  expect_identical(fit2$qr_helper_path, "compact_factor")
  expect_identical(fit1$qr_compact_factor_source, "native")
  expect_identical(fit2$qr_compact_factor_source, "native")
})

test_that("many_lm provides a repeated-response workflow surface", {
  set.seed(18)
  X_host <- cbind(1, matrix(rnorm(40), nrow = 10, ncol = 4))
  beta_host <- matrix(rnorm(15), nrow = 5, ncol = 3)
  Y_host <- X_host %*% beta_host + matrix(rnorm(30, sd = 1e-6), nrow = 10, ncol = 3)

  fit <- many_lm(adgeMatrix(X_host), Y_host, include_residuals = TRUE, method = "qr")
  coef_host <- base::qr.coef(base::qr(X_host), Y_host)
  resid_host <- base::qr.resid(base::qr(X_host), Y_host)
  rss_host <- colSums(resid_host^2)

  expect_s3_class(fit, "am_many_lm_fit")
  expect_identical(fit$responses, 3L)
  expect_identical(fit$observations, 10L)
  expect_identical(fit$method, "qr")
  expect_identical(fit$qr_representation, "base_qr")
  expect_identical(fit$qr_helper_path, "compact_factor")
  expect_true(fit$qr_compact_factor_available)
  expect_identical(fit$qr_compact_factor_source, "native")
  expect_equal(as.matrix(coef(fit)), coef_host, tolerance = 1e-8)
  expect_equal(as.matrix(residuals(fit)), resid_host, tolerance = 1e-8)
  expect_equal(fit$rss, rss_host, tolerance = 1e-8)
  expect_equal(fit$sigma2, rss_host / fit$df.residual, tolerance = 1e-8)
})

test_that("many_lm batched Y path matches per-column fits exactly", {
  set.seed(1807)
  X_host <- cbind(1, matrix(rnorm(200), nrow = 40, ncol = 5))
  Y_host <- matrix(rnorm(40 * 8), nrow = 40, ncol = 8)

  X_arg <- adgeMatrix(X_host)
  batched <- many_lm(X_arg, Y_host, method = "qr")
  coef_batched <- as.matrix(coef(batched))

  coef_per_col <- vapply(
    seq_len(ncol(Y_host)),
    function(j) {
      fit_j <- many_lm(X_arg, Y_host[, j, drop = FALSE], method = "qr")
      as.matrix(coef(fit_j))[, 1L]
    },
    numeric(ncol(X_host))
  )

  expect_identical(dim(coef_batched), c(ncol(X_host), ncol(Y_host)))
  expect_equal(coef_batched, coef_per_col, tolerance = .Machine$double.eps * 1e3)
})

test_that("many_lm supports weighted repeated-response fits", {
  set.seed(181)
  X_host <- cbind(1, matrix(rnorm(40), nrow = 10, ncol = 4))
  beta_host <- matrix(rnorm(15), nrow = 5, ncol = 3)
  Y_host <- X_host %*% beta_host + matrix(rnorm(30, sd = 1e-6), nrow = 10, ncol = 3)
  weights <- seq_len(nrow(X_host)) / nrow(X_host)
  sqrt_w <- sqrt(weights)
  Xw <- X_host * sqrt_w
  Yw <- Y_host * sqrt_w
  resid_host <- Y_host - X_host %*% base::qr.coef(base::qr(Xw), Yw)
  rss_host <- colSums(resid_host^2 * weights)

  fit <- many_lm(
    adgeMatrix(X_host),
    Y_host,
    weights = weights,
    include_residuals = TRUE,
    method = "qr"
  )

  expect_s3_class(fit, "am_many_lm_fit")
  expect_equal(as.matrix(coef(fit)), base::qr.coef(base::qr(Xw), Yw), tolerance = 1e-8)
  expect_equal(as.matrix(residuals(fit)), resid_host, tolerance = 1e-8)
  expect_equal(fit$weights, as.double(weights))
  expect_equal(fit$rss, rss_host, tolerance = 1e-8)
  expect_equal(fit$sigma2, rss_host / fit$df.residual, tolerance = 1e-8)
})

test_that("array_lm restores array-shaped response outputs", {
  set.seed(19)
  X_host <- cbind(1, matrix(rnorm(40), nrow = 10, ncol = 4))
  beta_host <- matrix(rnorm(20), nrow = 5, ncol = 4)
  Y_mat <- X_host %*% beta_host + matrix(rnorm(40, sd = 1e-6), nrow = 10, ncol = 4)
  Y_array <- array(Y_mat, dim = c(10, 2, 2))

  fit <- array_lm(
    adgeMatrix(X_host),
    Y_array,
    include_fitted = TRUE,
    include_residuals = TRUE,
    method = "qr"
  )

  expect_s3_class(fit, "am_array_lm_fit")
  expect_identical(fit$responses, 4L)
  expect_identical(fit$observations, 10L)
  expect_identical(fit$response_dims, c(2L, 2L))
  expect_identical(fit$qr_representation, "base_qr")
  expect_identical(fit$qr_helper_path, "compact_factor")
  expect_true(fit$qr_compact_factor_available)
  expect_identical(fit$qr_compact_factor_source, "native")
  expect_equal(dim(fitted(fit)), c(10L, 2L, 2L))
  expect_equal(dim(residuals(fit)), c(10L, 2L, 2L))
  expect_equal(dim(fit$rss), c(2L, 2L))
  expect_equal(dim(fit$sigma2), c(2L, 2L))
  expect_s4_class(coef(fit), "adgeMatrix")
  expect_equal(as.matrix(coef(fit)), base::qr.coef(base::qr(X_host), Y_mat), tolerance = 1e-8)
})

test_that("array_lm supports weighted array-response fits", {
  set.seed(191)
  X_host <- cbind(1, matrix(rnorm(40), nrow = 10, ncol = 4))
  beta_host <- matrix(rnorm(20), nrow = 5, ncol = 4)
  Y_mat <- X_host %*% beta_host + matrix(rnorm(40, sd = 1e-6), nrow = 10, ncol = 4)
  Y_array <- array(Y_mat, dim = c(10, 2, 2))
  weights <- seq_len(nrow(X_host)) / nrow(X_host)
  sqrt_w <- sqrt(weights)
  coef_host <- base::qr.coef(base::qr(X_host * sqrt_w), Y_mat * sqrt_w)

  fit <- array_lm(
    adgeMatrix(X_host),
    Y_array,
    weights = weights,
    include_fitted = TRUE,
    include_residuals = TRUE,
    method = "qr"
  )

  expect_s3_class(fit, "am_array_lm_fit")
  expect_equal(fit$weights, as.double(weights))
  expect_equal(dim(fitted(fit)), c(10L, 2L, 2L))
  expect_equal(dim(residuals(fit)), c(10L, 2L, 2L))
  expect_equal(dim(fit$rss), c(2L, 2L))
  expect_equal(dim(fit$sigma2), c(2L, 2L))
  expect_equal(as.matrix(coef(fit)), coef_host, tolerance = 1e-8)
})

test_that("ridge_fit matches penalized least-squares coefficients", {
  set.seed(14)
  X_host <- cbind(1, matrix(rnorm(48), nrow = 12, ncol = 4))
  beta_host <- matrix(rnorm(15), nrow = 5, ncol = 3)
  Y_host <- X_host %*% beta_host + matrix(rnorm(36, sd = 1e-6), nrow = 12, ncol = 3)
  lambda <- 0.75

  fit <- ridge_fit(adgeMatrix(X_host), Y_host, lambda = lambda, penalize_intercept = FALSE)
  penalty <- diag(c(0, rep(lambda, ncol(X_host) - 1L)))
  coef_host <- solve(crossprod(X_host) + penalty, crossprod(X_host, Y_host))

  expect_s3_class(fit, "ridge_fit")
  expect_s4_class(fit$coefficients, "adgeMatrix")
  expect_equal(as.matrix(coef(fit)), coef_host, tolerance = 1e-8)
  expect_equal(as.matrix(fitted(fit)), X_host %*% coef_host, tolerance = 1e-8)
  expect_equal(as.matrix(residuals(fit)), Y_host - X_host %*% coef_host, tolerance = 1e-8)
  expect_identical(fit$precision, "strict")
})

test_that("ridge_fit penalizes the first predictor when no intercept column is present", {
  set.seed(1401)
  X_host <- matrix(rnorm(48), nrow = 12, ncol = 4)
  beta_host <- matrix(rnorm(12), nrow = 4, ncol = 3)
  Y_host <- X_host %*% beta_host + matrix(rnorm(36, sd = 1e-6), nrow = 12, ncol = 3)
  lambda <- 0.75

  fit <- ridge_fit(adgeMatrix(X_host), Y_host, lambda = lambda, intercept = FALSE, penalize_intercept = FALSE)
  coef_host <- solve(crossprod(X_host) + diag(lambda, ncol(X_host)), crossprod(X_host, Y_host))

  expect_equal(as.matrix(coef(fit)), coef_host, tolerance = 1e-8)
})

test_that("ridge_fit leaves only the explicit intercept column unpenalized", {
  set.seed(1402)
  X_host <- matrix(rnorm(48), nrow = 12, ncol = 4)
  beta_host <- matrix(rnorm(15), nrow = 5, ncol = 3)
  X_design <- cbind(1, X_host)
  Y_host <- X_design %*% beta_host + matrix(rnorm(36, sd = 1e-6), nrow = 12, ncol = 3)
  lambda <- 0.5

  fit <- ridge_fit(adgeMatrix(X_host), Y_host, lambda = lambda, intercept = TRUE, penalize_intercept = FALSE)
  penalty <- diag(c(0, rep(lambda, ncol(X_host))))
  coef_host <- solve(crossprod(X_design) + penalty, crossprod(X_design, Y_host))

  expect_equal(as.matrix(coef(fit)), coef_host, tolerance = 1e-8)
})

test_that("ridge_fit reuses cached shared-X work across repeated fits", {
  counter <- new.env(parent = emptyenv())
  counter$crossprod <- 0L
  counter$ewise <- 0L
  counter$solve <- 0L

  with_registered_backend(
    "ridge_cache_backend",
    make_recording_backend(
      counter,
      supported_ops = c("crossprod", "ewise"),
      features = c("dense_f32")
    ),
    {
      set.seed(15)
      X <- adgeMatrix(
        matrix(rnorm(45), nrow = 15, ncol = 3),
        preferred_backend = "ridge_cache_backend",
        precision = "fast"
      )
      Y1 <- matrix(rnorm(30), nrow = 15, ncol = 2)
      Y2 <- matrix(rnorm(30), nrow = 15, ncol = 2)

      fit1 <- ridge_fit(X, Y1, lambda = 0.5, include_fitted = FALSE, include_residuals = FALSE)
      fit2 <- ridge_fit(X, Y2, lambda = 1.0, include_fitted = FALSE, include_residuals = FALSE)

      expect_false(fit1$cache_reused)
      expect_true(fit2$cache_reused)
      expect_identical(fit1$cache_key, fit2$cache_key)
      expect_equal(counter$crossprod, 3L)
      expect_true(counter$ewise >= 2L)
      expect_identical(counter$solve, 0L)
    }
  )
})

test_that("wls_fit matches weighted least-squares coefficients", {
  set.seed(20)
  X_host <- cbind(1, matrix(rnorm(48), nrow = 12, ncol = 4))
  beta_host <- matrix(rnorm(15), nrow = 5, ncol = 3)
  Y_host <- X_host %*% beta_host + matrix(rnorm(36, sd = 1e-6), nrow = 12, ncol = 3)
  weights <- seq_len(nrow(X_host)) / nrow(X_host)
  sqrt_w <- sqrt(weights)
  Xw <- X_host * sqrt_w
  Yw <- Y_host * sqrt_w

  fit <- wls_fit(adgeMatrix(X_host), Y_host, weights = weights, method = "qr")
  coef_host <- base::qr.coef(base::qr(Xw), Yw)

  expect_s3_class(fit, "wls_fit")
  expect_s4_class(fit$coefficients, "adgeMatrix")
  expect_identical(fit$qr_representation, "base_qr")
  expect_identical(fit$qr_helper_path, "compact_factor")
  expect_true(fit$qr_compact_factor_available)
  expect_identical(fit$qr_compact_factor_source, "native")
  expect_equal(as.matrix(coef(fit)), coef_host, tolerance = 1e-8)
  expect_identical(fit$precision, "strict")
})

test_that("wls_fit reuses cached weighted shared-X work across repeated fits", {
  set.seed(22)
  X <- adgeMatrix(matrix(rnorm(45), nrow = 15, ncol = 3))
  Y1 <- matrix(rnorm(30), nrow = 15, ncol = 2)
  Y2 <- matrix(rnorm(30), nrow = 15, ncol = 2)
  weights <- rep(1, nrow(X))

  fit1 <- wls_fit(X, Y1, weights = weights, include_fitted = FALSE, include_residuals = FALSE, method = "qr")
  fit2 <- wls_fit(X, Y2, weights = weights, include_fitted = FALSE, include_residuals = FALSE, method = "qr")

  expect_false(fit1$cache_reused)
  expect_true(fit2$cache_reused)
  expect_identical(fit1$cache_key, fit2$cache_key)
  expect_identical(fit1$method, "qr")
  expect_identical(fit1$qr_compact_factor_source, "native")
  expect_identical(fit2$qr_compact_factor_source, "native")
})

test_that("covariance matches stats::cov on dense inputs", {
  set.seed(23)
  X_host <- matrix(rnorm(72), nrow = 12, ncol = 6)
  fit <- covariance(adgeMatrix(X_host))

  expect_s4_class(fit, "adgeMatrix")
  expect_equal(as.matrix(fit), stats::cov(X_host), tolerance = 1e-8)
})

test_that("covariance supports weighted covariance", {
  set.seed(24)
  X_host <- matrix(rnorm(72), nrow = 12, ncol = 6)
  weights <- seq_len(nrow(X_host)) / nrow(X_host)
  means <- colSums(X_host * weights) / sum(weights)
  centered <- sweep(X_host, 2L, means, FUN = "-")
  cov_host <- crossprod(centered * sqrt(weights)) / (sum(weights) - 1)

  fit <- covariance(adgeMatrix(X_host), weights = weights)

  expect_s4_class(fit, "adgeMatrix")
  expect_equal(as.matrix(fit), cov_host, tolerance = 1e-8)
})

test_that("covariance supports blockwise evaluation", {
  set.seed(241)
  X_host <- matrix(rnorm(96), nrow = 12, ncol = 8)

  fit_full <- covariance(adgeMatrix(X_host))
  fit_block <- covariance(adgeMatrix(X_host), block_size = 3L)

  expect_s4_class(fit_block, "adgeMatrix")
  expect_equal(as.matrix(fit_block), as.matrix(fit_full), tolerance = 1e-8)
})

test_that("covariance supports weighted blockwise evaluation", {
  set.seed(242)
  X_host <- matrix(rnorm(96), nrow = 12, ncol = 8)
  weights <- seq_len(nrow(X_host)) / nrow(X_host)

  fit_full <- covariance(adgeMatrix(X_host), weights = weights)
  fit_block <- covariance(adgeMatrix(X_host), weights = weights, block_size = 3L)

  expect_s4_class(fit_block, "adgeMatrix")
  expect_equal(as.matrix(fit_block), as.matrix(fit_full), tolerance = 1e-8)
})

test_that("correlation matches stats::cor on dense inputs", {
  set.seed(25)
  X_host <- matrix(rnorm(72), nrow = 12, ncol = 6)
  fit <- correlation(adgeMatrix(X_host))

  expect_s4_class(fit, "adgeMatrix")
  expect_equal(as.matrix(fit), stats::cor(X_host), tolerance = 1e-8)
})

test_that("correlation supports blockwise evaluation", {
  set.seed(251)
  X_host <- matrix(rnorm(96), nrow = 12, ncol = 8)

  fit_full <- correlation(adgeMatrix(X_host))
  fit_block <- correlation(adgeMatrix(X_host), block_size = 3L)

  expect_s4_class(fit_block, "adgeMatrix")
  expect_equal(as.matrix(fit_block), as.matrix(fit_full), tolerance = 1e-8)
})

test_that("qr returns an amQR object with QR helper methods", {
  x_host <- matrix(c(1, 1, 1, 1, 2, 3), nrow = 3, ncol = 2)
  y_host <- matrix(c(1, 2, 4), nrow = 3, ncol = 1)
  x <- adgeMatrix(x_host)

  fac <- qr(x)
  fac_base <- base::qr(x_host)

  expect_s3_class(fac, "amQR")
  expect_s3_class(fac, "amDenseQR")
  expect_identical(fac$backend, "cpu")
  expect_identical(fac$precision, "strict")
  expect_identical(fac$representation, "base_qr")
  expect_identical(dim(fac), c(3L, 2L))
  expect_true(fac$thin)
  expect_false(fac$pivoted)
  expect_false(fac$q_materialized)
  expect_false(fac$r_materialized)

  info <- qr_info(fac)
  expect_identical(info$rank, 2L)
  expect_identical(info$dim, c(3L, 2L))
  expect_true(info$thin)
  expect_false(info$pivoted)
  expect_identical(info$representation, "base_qr")
  expect_identical(info$backend, "cpu")
  expect_identical(info$precision, "strict")
  expect_true(info$compact_factor_available)
  expect_identical(info$compact_factor_source, "native")
  expect_false(info$q_materialized)
  expect_false(info$r_materialized)

  q_fit <- qr.Q(fac)
  r_fit <- qr.R(fac)
  solve_fit <- qr.solve(fac, y_host)
  inv_fit <- qr.solve(qr(adgeMatrix(diag(2))))
  coef_fit <- qr.coef(fac, y_host)
  qty_fit <- qr.qty(fac, y_host)
  qy_fit <- qr.qy(fac, y_host)
  fitted_fit <- qr.fitted(fac, y_host)
  resid_fit <- qr.resid(fac, y_host)

  expect_s4_class(q_fit, "adgeMatrix")
  expect_s4_class(r_fit, "adgeMatrix")
  expect_s4_class(solve_fit, "adgeMatrix")
  expect_s4_class(inv_fit, "adgeMatrix")
  expect_s4_class(coef_fit, "adgeMatrix")
  expect_s4_class(qty_fit, "adgeMatrix")
  expect_s4_class(qy_fit, "adgeMatrix")
  expect_s4_class(fitted_fit, "adgeMatrix")
  expect_s4_class(resid_fit, "adgeMatrix")

  expect_equal(as.matrix(q_fit), base::qr.Q(fac_base), tolerance = 1e-10)
  expect_equal(as.matrix(r_fit), base::qr.R(fac_base), tolerance = 1e-10)
  expect_equal(as.matrix(solve_fit), base::qr.solve(fac_base, y_host), tolerance = 1e-10)
  expect_equal(as.matrix(inv_fit), base::qr.solve(base::qr(diag(2))), tolerance = 1e-10)
  expect_equal(as.matrix(coef_fit), base::qr.coef(fac_base, y_host), tolerance = 1e-10)
  expect_equal(as.matrix(qty_fit), base::qr.qty(fac_base, y_host), tolerance = 1e-10)
  expect_equal(as.matrix(qy_fit), base::qr.qy(fac_base, y_host), tolerance = 1e-10)
  expect_equal(as.matrix(fitted_fit), base::qr.fitted(fac_base, y_host), tolerance = 1e-10)
  expect_equal(as.matrix(resid_fit), base::qr.resid(fac_base, y_host), tolerance = 1e-10)
})

test_that("qr.solve on rectangular amQR matches base QR solve", {
  x_host <- matrix(c(1, 1, 1, 1, 2, 3), nrow = 3, ncol = 2)
  y_host <- matrix(c(1, 2, 4), nrow = 3, ncol = 1)
  fac <- qr(adgeMatrix(x_host))
  fac_base <- base::qr(x_host)

  expect_equal(as.matrix(qr.solve(fac, y_host)), base::qr.solve(fac_base, y_host), tolerance = 1e-10)
  expect_error(qr.solve(fac), "only square matrices can be inverted")
})

test_that("qr_info reports explicit backend QR metadata", {
  if (is.null(optional_backend_namespace("amatrix.mlx"))) {
    skip("Optional backend package 'amatrix.mlx' is not installed")
  }
  if (!amatrix_backend_status("mlx")$available) {
    skip("mlx backend not available")
  }

  old <- options(amatrix.mlx.available = TRUE)
  on.exit(options(old), add = TRUE)

  x_host <- matrix(rnorm(30), nrow = 6, ncol = 5)
  fac <- qr(adgeMatrix(x_host, preferred_backend = "mlx", precision = "fast"))
  info_before <- qr_info(fac)
  expect_identical(info_before$representation, "explicit_qr")
  expect_true(info_before$compact_factor_available)
  expect_identical(info_before$compact_factor_source, "reconstructable")
  expect_false(info_before$compact_factor_materialized)
  expect_identical(info_before$helper_path, "native_resident_backend")
  expect_false(info_before$q_materialized)
  expect_true(info_before$r_materialized)

  invisible(qr.qty(fac, matrix(1, nrow = 6, ncol = 1)))
  info_after <- qr_info(fac)

  expect_s3_class(fac, "amQR")
  expect_identical(info_after$representation, "explicit_qr")
  expect_identical(info_after$dim, c(6L, 5L))
  expect_true(info_after$thin)
  expect_false(info_after$pivoted)
  expect_identical(info_after$helper_path, "native_resident_backend")
  expect_false(info_after$q_materialized)
  expect_true(info_after$r_materialized)
  expect_true(info_after$compact_factor_available)
  expect_identical(info_after$compact_factor_source, "reconstructable")
  expect_false(info_after$compact_factor_materialized)
})

test_that("qr_info reports compact MLX QR representation when requested", {
  if (is.null(optional_backend_namespace("amatrix.mlx"))) {
    skip("Optional backend package 'amatrix.mlx' is not installed")
  }
  if (!amatrix_backend_status("mlx")$available) {
    skip("mlx backend not available")
  }

  old <- options(amatrix.mlx.available = TRUE, amatrix.mlx.qr_helper_mode = "compact")
  on.exit(options(old), add = TRUE)

  x_host <- matrix(rnorm(30), nrow = 6, ncol = 5)
  fac <- qr(adgeMatrix(x_host, preferred_backend = "mlx", precision = "fast"))
  info <- qr_info(fac)

  expect_identical(info$representation, "mlx_compact_qr")
  expect_identical(info$helper_path, "compact_mlx_factor")
  expect_true(info$compact_factor_available)
  expect_identical(info$compact_factor_source, "host_compact")
  expect_false(info$compact_factor_materialized)
  expect_false(info$q_materialized)
  expect_true(info$r_materialized)
  expect_equal(as.matrix(qr.coef(fac, x_host)), base::qr.coef(base::qr(x_host), x_host), tolerance = 1e-8)
  expect_equal(as.matrix(qr.qty(fac, x_host)), base::qr.qty(base::qr(x_host), x_host), tolerance = 1e-8)

  info_after <- qr_info(fac)
  expect_true(info_after$compact_factor_materialized)
})

test_that("qr-backed fits report bridge-compact provenance for explicit backend QR", {
  if (is.null(optional_backend_namespace("amatrix.mlx"))) {
    skip("Optional backend package 'amatrix.mlx' is not installed")
  }
  if (!amatrix_backend_status("mlx")$available) {
    skip("mlx backend not available")
  }

  old <- options(amatrix.mlx.available = TRUE)
  on.exit(options(old), add = TRUE)

  set.seed(902)
  X_host <- matrix(rnorm(30), nrow = 6, ncol = 5)
  Y_host <- matrix(rnorm(18), nrow = 6, ncol = 3)
  Y_array <- array(Y_host, dim = c(6, 3, 1))
  X <- adgeMatrix(X_host, preferred_backend = "mlx", precision = "fast")

  lm_fit <- lm_fit(X, Y_host, method = "qr")
  many_fit <- many_lm(X, Y_host, method = "qr")
  array_fit <- array_lm(X, Y_array, method = "qr")

  expect_identical(lm_fit$qr_representation, "explicit_qr")
  expect_identical(lm_fit$qr_helper_path, "native_resident_backend")
  expect_true(lm_fit$qr_compact_factor_available)
  expect_identical(lm_fit$qr_compact_factor_source, "reconstructable")

  expect_identical(many_fit$qr_representation, "explicit_qr")
  expect_identical(many_fit$qr_helper_path, "native_resident_backend")
  expect_true(many_fit$qr_compact_factor_available)
  expect_identical(many_fit$qr_compact_factor_source, "reconstructable")

  expect_identical(array_fit$qr_representation, "explicit_qr")
  expect_identical(array_fit$qr_helper_path, "native_resident_backend")
  expect_true(array_fit$qr_compact_factor_available)
  expect_identical(array_fit$qr_compact_factor_source, "reconstructable")
})

test_that("tall-skinny compact MLX QR uses tsqr-blocked provenance", {
  if (is.null(optional_backend_namespace("amatrix.mlx"))) {
    skip("Optional backend package 'amatrix.mlx' is not installed")
  }
  if (!amatrix_backend_status("mlx")$available) {
    skip("mlx backend not available")
  }

  old <- options(
    amatrix.mlx.available = TRUE,
    amatrix.mlx.qr_helper_mode = "compact",
    amatrix.mlx.qr_compact_method = "tsqr",
    amatrix.mlx.qr_tsqr_block_rows = 8L
  )
  on.exit(options(old), add = TRUE)

  set.seed(903)
  X_host <- matrix(rnorm(24 * 4), nrow = 24, ncol = 4)
  Y_host <- matrix(rnorm(24 * 3), nrow = 24, ncol = 3)
  X <- adgeMatrix(X_host, preferred_backend = "mlx", precision = "fast")

  fac <- qr(X)
  info <- qr_info(fac)
  many_fit <- many_lm(X, Y_host, method = "qr", cache = TRUE)
  fac_base <- base::qr(X_host)

  expect_identical(info$representation, "mlx_compact_qr")
  expect_identical(info$helper_path, "compact_mlx_factor")
  expect_identical(info$compact_factor_source, "tsqr_blocked")
  expect_false(info$q_materialized)
  expect_false(info$r_materialized)
  qty <- qr.qty(fac, Y_host)
  expect_equal(abs(diag(as.matrix(qr.R(fac)))), abs(diag(base::qr.R(fac_base))), tolerance = 1e-6)
  expect_equal(as.matrix(qr.coef(fac, Y_host)), base::qr.coef(fac_base, Y_host), tolerance = 1e-6)
  expect_equal(as.matrix(qr.fitted(fac, Y_host)), base::qr.fitted(fac_base, Y_host), tolerance = 1e-6)
  expect_equal(as.matrix(qr.resid(fac, Y_host)), base::qr.resid(fac_base, Y_host), tolerance = 1e-6)
  expect_identical(dim(qty), dim(Y_host))
  expect_equal(as.matrix(qr.qy(fac, qty)), Y_host, tolerance = 1e-8)

  expect_identical(many_fit$qr_representation, "mlx_compact_qr")
  expect_identical(many_fit$qr_helper_path, "compact_mlx_factor")
  expect_identical(many_fit$qr_compact_factor_source, "tsqr_blocked")
  expect_equal(as.matrix(coef(many_fit)), base::qr.coef(fac_base, Y_host), tolerance = 1e-6)
})

test_that("MLX QR cache keys distinguish native and compact strategies", {
  if (is.null(optional_backend_namespace("amatrix.mlx"))) {
    skip("Optional backend package 'amatrix.mlx' is not installed")
  }

  set.seed(904)
  X_host <- matrix(rnorm(1024 * 128), nrow = 1024, ncol = 128)
  Y_host <- matrix(rnorm(1024 * 32), nrow = 1024, ncol = 32)
  X <- adgeMatrix(X_host, preferred_backend = "mlx", precision = "fast")

  old <- options(
    amatrix.mlx.available = TRUE,
    amatrix.mlx.qr_tsqr_block_rows = 256L
  )
  on.exit(options(old), add = TRUE)

  options(amatrix.mlx.qr_helper_mode = "native")
  fit_native <- many_lm(X, Y_host, method = "qr", cache = TRUE, include_residuals = FALSE)

  options(amatrix.mlx.qr_helper_mode = "compact", amatrix.mlx.qr_compact_method = "tsqr")
  fit_compact <- many_lm(X, Y_host, method = "qr", cache = TRUE, include_residuals = FALSE)

  expect_identical(fit_native$qr_helper_path, "native_resident_backend")
  expect_identical(fit_native$qr_compact_factor_source, "reconstructable")
  expect_identical(fit_compact$qr_helper_path, "compact_mlx_factor")
  expect_identical(fit_compact$qr_compact_factor_source, "tsqr_blocked")
  expect_false(identical(fit_native$cache_key, fit_compact$cache_key))
  expect_equal(as.matrix(coef(fit_compact)), base::qr.coef(base::qr(X_host), Y_host), tolerance = 1e-6)
})

test_that("am_qr is the deterministic internal QR wrapper", {
  x_host <- matrix(c(1, 0, 0, 1), nrow = 2)
  x <- adgeMatrix(x_host, preferred_backend = "cpu")
  fac <- am_qr(x)

  expect_s3_class(fac, "amQR")
  expect_equal(as.matrix(qr.coef(fac, x_host)), diag(2), tolerance = 1e-10)
})

test_that("random dimnames and coercions round-trip through host representations", {
  set.seed(4)

  for (iter in seq_len(10)) {
    nr <- sample(2:5, 1)
    nc <- sample(2:5, 1)
    x <- adgeMatrix(matrix(rnorm(nr * nc), nrow = nr, ncol = nc))
    rn <- paste0("r", seq_len(nr), "_", iter)
    cn <- paste0("c", seq_len(nc), "_", iter)

    dimnames(x) <- list(rn, cn)

    expect_identical(rownames(as.matrix(x)), rn)
    expect_identical(colnames(as.matrix(x)), cn)
    expect_identical(dimnames(as.array(x)), list(rn, cn))
  }
})

test_that("zero-dimension matrices preserve semantics", {
  dense0 <- adgeMatrix(matrix(numeric(), nrow = 0, ncol = 3))
  sparse0 <- adgCMatrix(matrix(numeric(), nrow = 3, ncol = 0))

  expect_s4_class(dense0 + 1, "adgeMatrix")
  expect_equal(dim(dense0 + 1), c(0, 3))
  expect_equal(rowSums(dense0), numeric())
  expect_equal(colSums(dense0), numeric(3))
  expect_equal(dim(t(dense0)), c(3, 0))

  expect_true(inherits(sparse0, "aMatrix"))
  expect_equal(dim(sparse0), c(3, 0))
  expect_equal(rowSums(sparse0), numeric(3))
  expect_equal(colSums(sparse0), numeric())
  expect_equal(dim(t(sparse0)), c(0, 3))
})

test_that("Ops preserve NA and NaN semantics", {
  x_host <- matrix(c(1, NA, NaN, 4), nrow = 2)
  y_host <- matrix(c(2, 3, 4, NA), nrow = 2)
  x <- adgeMatrix(x_host)
  y <- adgeMatrix(y_host)

  expect_equal(as.matrix(x + y), x_host + y_host)
  expect_equal(as.matrix(x * y), x_host * y_host)
  expect_identical(as.matrix(x > y), x_host > y_host)
  expect_identical(as.matrix(x == y), x_host == y_host)
})
