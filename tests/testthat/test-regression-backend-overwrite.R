# Regression repros for backend-overwrite state invalidation.
# Seed: deterministic literal matrices (no RNG required)
# Shape: 2 x 2 / 4 x 4 dense matrices
# Backend: overwrite_state_backend
# Precision mode: strict
# Dispatch path: backend overwrite -> planner / health / model-cache reuse
# R/platform: R version 4.5.1 (2025-06-13) | aarch64-apple-darwin20
# Issues: amatrix-ubq, amatrix-2nh, amatrix-fqh

test_that("backend overwrite clears stale calibration thresholds for that backend", {
  state <- amatrix:::.amatrix_state
  counter <- new.env(parent = emptyenv())
  backend_name <- "overwrite_state_backend"

  old_calibration <- state$calibration
  on.exit(state$calibration <- old_calibration, add = TRUE)

  with_registered_backend(
    backend_name,
    make_recording_backend(counter, supported_ops = c("matmul")),
    {
      state$calibration <- list(
        version = "1",
        calibrated_at = Sys.time(),
        thresholds = setNames(list(list(gemm = Inf)), backend_name),
        results = data.frame(
          backend = backend_name,
          op = "gemm",
          op_base = "matmul",
          nrow = 4L,
          ncol = 4L,
          elements = 16L,
          nnz = NA_integer_,
          density = NA_real_,
          density_bucket = NA_character_,
          cpu_ms = 1,
          gpu_ms = 2,
          margin = NA_real_,
          gpu_wins = FALSE,
          stringsAsFactors = FALSE
        )
      )

      x <- adgeMatrix(matrix(1:16 * 1.0, 4, 4), preferred_backend = backend_name)
      expect_identical(amatrix_backend_plan(x, "matmul", y = diag(4))$chosen, "cpu")

      amatrix_register_backend(
        backend_name,
        make_recording_backend(new.env(parent = emptyenv()), supported_ops = c("matmul")),
        overwrite = TRUE
      )

      plan <- amatrix_backend_plan(x, "matmul", y = diag(4))
      expect_identical(plan$chosen, backend_name)
      expect_null(state$calibration$thresholds[[backend_name]])
      expect_false(any(state$calibration$results$backend == backend_name))
    }
  )
})

test_that("backend overwrite clears stale backend health state", {
  backend_name <- "overwrite_health_backend"

  with_registered_backend(
    backend_name,
    make_recording_backend(new.env(parent = emptyenv()), supported_ops = c("matmul")),
    {
      amatrix:::.amatrix_backend_health_mark(backend_name, "unhealthy", "stale")
      expect_identical(amatrix:::.amatrix_backend_health_get(backend_name)$status, "unhealthy")

      amatrix_register_backend(
        backend_name,
        make_recording_backend(new.env(parent = emptyenv()), supported_ops = c("matmul")),
        overwrite = TRUE
      )

      rec <- amatrix:::.amatrix_backend_health_get(backend_name)
      expect_identical(rec$status, "unprobed")
      expect_true(is.na(rec$reason))
    }
  )
})

test_that("backend overwrite invalidates cached chol and svd factors", {
  backend_name <- "overwrite_cache_backend"
  spd <- matrix(c(5, 1, 1, 4) * 1.0, 2, 2)
  counter_a <- new.env(parent = emptyenv())
  counter_b <- new.env(parent = emptyenv())
  cache_env <- amatrix:::.amatrix_state$model_cache

  with_registered_backend(
    backend_name,
    make_recording_backend(counter_a, supported_ops = c("chol", "svd")),
    {
      x <- adgeMatrix(spd, preferred_backend = backend_name)

      fac_1 <- chol_factor(x)
      svd_1 <- svd_factor(x, k = 2L)
      expect_true(length(ls(envir = cache_env, all.names = FALSE)) > 0L)

      fac_cached <- chol_factor(x)
      svd_cached <- svd_factor(x, k = 2L)
      expect_identical(fac_cached@backend, backend_name)
      expect_identical(svd_cached@backend, backend_name)

      amatrix_register_backend(
        backend_name,
        make_recording_backend(counter_b, supported_ops = c("chol", "svd")),
        overwrite = TRUE
      )
      expect_length(ls(envir = cache_env, all.names = FALSE), 0L)

      fac_2 <- chol_factor(x)
      svd_2 <- svd_factor(x, k = 2L)
      expect_identical(counter_b$chol, 1L)
      expect_identical(counter_b$svd, 1L)
      expect_equal(as.matrix(fac_2), base::chol(spd), tolerance = 1e-10)
      expect_equal(svd_2@d, base::svd(spd)$d[1:2], tolerance = 1e-10)
      expect_false(identical(fac_1@factor_obj, fac_2@factor_obj))
      expect_true(length(ls(envir = cache_env, all.names = FALSE)) > 0L)
    }
  )
})
