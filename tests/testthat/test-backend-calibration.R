test_that("matmul calibration distinguishes gemm from gemv", {
  counter <- new.env(parent = emptyenv())
  state <- get(".amatrix_state", envir = asNamespace("amatrix"))

  with_registered_backend("recording_calibration", make_recording_backend(counter, supported_ops = c("matmul")), {
    old_calibration <- state$calibration
    on.exit(state$calibration <- old_calibration, add = TRUE)

    state$calibration <- list(
      version = "1",
      calibrated_at = Sys.time(),
      thresholds = list(
        recording_calibration = list(
          gemm = 0L,
          gemv = Inf
        )
      ),
      results = data.frame()
    )

    x <- adgeMatrix(matrix(1:16, nrow = 4), preferred_backend = "recording_calibration")

    gemm_plan <- amatrix_backend_plan(x, "matmul", y = diag(4))
    gemv_plan <- amatrix_backend_plan(x, "matmul", y = matrix(1, nrow = 4, ncol = 1))

    expect_identical(gemm_plan$chosen, "recording_calibration")
    expect_identical(gemv_plan$chosen, "cpu")
  })
})

test_that("amatrix_calibrate benchmarks gemv and reductions separately", {
  counter <- new.env(parent = emptyenv())
  state <- get(".amatrix_state", envir = asNamespace("amatrix"))

  with_registered_backend(
    "recording_calibration_ops",
    make_recording_backend(counter, supported_ops = c("matmul", "rowSums", "colSums")),
    {
      old_calibration <- state$calibration
      on.exit(state$calibration <- old_calibration, add = TRUE)

      cal <- amatrix_calibrate(
        backend = "recording_calibration_ops",
        ops = c("matmul", "gemv", "rowSums", "colSums"),
        sizes = list(c(4L, 4L)),
        n_reps = 1L,
        persist = FALSE,
        quiet = TRUE
      )

      expect_true(all(c("gemm", "gemv", "rowSums", "colSums") %in% cal$results$op))
      expect_true(all(c("gemm", "gemv", "rowSums", "colSums") %in% names(cal$thresholds$recording_calibration_ops)))
    }
  )
})

test_that("amatrix_calibrate benchmarks solve and svd separately", {
  counter <- new.env(parent = emptyenv())
  state <- get(".amatrix_state", envir = asNamespace("amatrix"))

  with_registered_backend(
    "recording_factor_calibration",
    make_recording_backend(counter, supported_ops = c("solve", "svd")),
    {
      old_calibration <- state$calibration
      on.exit(state$calibration <- old_calibration, add = TRUE)

      cal <- amatrix_calibrate(
        backend = "recording_factor_calibration",
        ops = c("solve", "svd"),
        sizes = list(c(8L, 4L)),
        n_reps = 1L,
        persist = FALSE,
        quiet = TRUE
      )

      expect_true(all(c("solve", "svd") %in% cal$results$op))
      expect_true(all(c("solve", "svd") %in% names(cal$thresholds$recording_factor_calibration)))
    }
  )
})

test_that("sparse matmul calibration distinguishes density buckets", {
  counter <- new.env(parent = emptyenv())
  state <- get(".amatrix_state", envir = asNamespace("amatrix"))

  with_registered_backend(
    "recording_sparse_calibration",
    make_recording_backend(counter, supported_ops = c("matmul"), supports_sparse_matmul = TRUE),
    {
      old_calibration <- state$calibration
      on.exit(state$calibration <- old_calibration, add = TRUE)

      sparse_keys <- setNames(
        list(0L, Inf, 0L, Inf),
        c("spmv:ultra_sparse", "spmv:semi_dense", "spmm:ultra_sparse", "spmm:semi_dense")
      )

      state$calibration <- list(
        version = "1",
        calibrated_at = Sys.time(),
        thresholds = list(recording_sparse_calibration = sparse_keys),
        results = data.frame()
      )

      x_sparse <- as_adgCMatrix(
        amatrix:::.amatrix_sparse_benchmark_matrix(128L, 128L, density = 0.01),
        preferred_backend = "recording_sparse_calibration"
      )
      x_semidense <- as_adgCMatrix(
        amatrix:::.amatrix_sparse_benchmark_matrix(128L, 128L, density = 0.10),
        preferred_backend = "recording_sparse_calibration"
      )

      spmv_sparse <- amatrix_backend_plan(x_sparse, "matmul", y = matrix(1, nrow = 128, ncol = 1))
      spmv_semidense <- amatrix_backend_plan(x_semidense, "matmul", y = matrix(1, nrow = 128, ncol = 1))
      spmm_sparse <- amatrix_backend_plan(x_sparse, "matmul", y = matrix(1, nrow = 128, ncol = 4))
      spmm_semidense <- amatrix_backend_plan(x_semidense, "matmul", y = matrix(1, nrow = 128, ncol = 4))

      expect_identical(spmv_sparse$chosen, "recording_sparse_calibration")
      expect_identical(spmv_semidense$chosen, "cpu")
      expect_identical(spmm_sparse$chosen, "recording_sparse_calibration")
      expect_identical(spmm_semidense$chosen, "cpu")
    }
  )
})

test_that("amatrix_calibrate benchmarks spmv and spmm separately", {
  counter <- new.env(parent = emptyenv())
  state <- get(".amatrix_state", envir = asNamespace("amatrix"))

  with_registered_backend(
    "recording_sparse_ops",
    make_recording_backend(counter, supported_ops = c("matmul"), supports_sparse_matmul = TRUE),
    {
      old_calibration <- state$calibration
      on.exit(state$calibration <- old_calibration, add = TRUE)

      cal <- amatrix_calibrate(
        backend = "recording_sparse_ops",
        ops = c("spmv", "spmm"),
        sizes = list(c(64L, 64L)),
        sparse_densities = c(0.01, 0.10),
        n_reps = 1L,
        persist = FALSE,
        quiet = TRUE
      )

      expect_true(any(startsWith(cal$results$op, "spmv:")))
      expect_true(any(startsWith(cal$results$op, "spmm:")))
      expect_true(any(startsWith(names(cal$thresholds$recording_sparse_ops), "spmv:")))
      expect_true(any(startsWith(names(cal$thresholds$recording_sparse_ops), "spmm:")))
    }
  )
})
