# Invariant-driven planning and calibration checks.
#
# Focus:
# - threshold derivation picks the earliest monotone tail win
# - dispatch signatures / workloads are shape-stable and respect rhs width
# - calibration gating affects cold paths but never resident reuse

suppressPackageStartupMessages(library(amatrix))

test_that("threshold derivation returns earliest monotone winning tail", {
  cases <- list(
    list(
      name = "all_win",
      results = data.frame(op = "gemm", elements = c(32L, 64L, 128L), gpu_wins = c(TRUE, TRUE, TRUE)),
      expected = 32L
    ),
    list(
      name = "non_monotone_then_tail",
      results = data.frame(op = "gemm", elements = c(32L, 64L, 128L, 256L), gpu_wins = c(TRUE, FALSE, TRUE, TRUE)),
      expected = 128L
    ),
    list(
      name = "tail_only_last",
      results = data.frame(op = "gemm", elements = c(32L, 64L, 128L), gpu_wins = c(TRUE, FALSE, TRUE)),
      expected = 128L
    ),
    list(
      name = "never",
      results = data.frame(op = "gemm", elements = c(32L, 64L, 128L), gpu_wins = c(FALSE, FALSE, FALSE)),
      expected = Inf
    )
  )

  for (case in cases) {
    thresholds <- amatrix:::.amatrix_derive_thresholds(case$results, "gemm")
    expect_identical(thresholds$gemm, case$expected, info = case$name)
  }
})

test_that("dispatch signature and workload distinguish gemv/gemm and sparse buckets", {
  dense <- as_adgeMatrix(matrix(rnorm(12L), nrow = 3L, ncol = 4L))
  sparse_ultra <- as_adgCMatrix(amatrix:::.amatrix_sparse_benchmark_matrix(32L, 32L, density = 0.01))
  sparse_semi <- as_adgCMatrix(amatrix:::.amatrix_sparse_benchmark_matrix(32L, 32L, density = 0.10))

  gemv_rhs <- matrix(1, nrow = 4L, ncol = 1L)
  gemm_rhs <- matrix(1, nrow = 4L, ncol = 3L)

  expect_identical(amatrix:::.amatrix_dispatch_signature(dense, "matmul", y = gemv_rhs), "gemv")
  expect_identical(amatrix:::.amatrix_dispatch_signature(dense, "matmul", y = gemm_rhs), "gemm")
  expect_true(amatrix:::.amatrix_dispatch_workload(dense, "matmul", y = gemm_rhs) >
                amatrix:::.amatrix_dispatch_workload(dense, "matmul", y = gemv_rhs))

  expect_identical(amatrix:::.amatrix_dispatch_signature(sparse_ultra, "matmul", y = gemv_rhs), "spmv:sparse")
  expect_identical(amatrix:::.amatrix_dispatch_signature(sparse_ultra, "matmul", y = gemm_rhs), "spmm:sparse")
  expect_identical(amatrix:::.amatrix_dispatch_signature(sparse_semi, "matmul", y = gemv_rhs), "spmv:semi_dense")
  expect_identical(amatrix:::.amatrix_dispatch_signature(sparse_semi, "matmul", y = gemm_rhs), "spmm:semi_dense")
})

test_that("calibration gates cold path but not resident reuse", {
  counter <- new.env(parent = emptyenv())
  state <- get(".amatrix_state", envir = asNamespace("amatrix"))

  with_registered_backend(
    "planning_invariant_backend",
    make_recording_backend(
      counter,
      supported_ops = c("matmul"),
      cold_supported_ops = c("matmul"),
      resident_supported_ops = c("matmul")
    ),
    {
      old_calibration <- state$calibration
      on.exit(state$calibration <- old_calibration, add = TRUE)

      state$calibration <- list(
        version = "1",
        calibrated_at = Sys.time(),
        thresholds = list(planning_invariant_backend = list(gemm = Inf)),
        results = data.frame()
      )

      x <- adgeMatrix(matrix(1:16, nrow = 4L), preferred_backend = "planning_invariant_backend")
      y <- diag(4L)

      cold_plan <- amatrix_backend_plan(x, "matmul", y = y)
      resident_x <- amatrix_bind_resident(x, backend = "planning_invariant_backend")
      resident_plan <- amatrix_backend_plan(resident_x, "matmul", y = y)

      expect_identical(cold_plan$chosen, "cpu")
      expect_identical(cold_plan$chosen_path, "cold")
      expect_identical(resident_plan$chosen, "planning_invariant_backend")
      expect_identical(resident_plan$chosen_path, "resident")
    }
  )
})
