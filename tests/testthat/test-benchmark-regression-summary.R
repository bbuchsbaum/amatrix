benchmark_regression_context <- function() {
  # tools/ is not shipped in the installed package; skip when R CMD check
  # runs tests from an installed-pkg temp dir.
  repo_root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    winslash = "/", mustWork = FALSE
  )
  helper_path <- file.path(repo_root, "tools", "benchmark-helpers.R")
  script_path <- file.path(repo_root, "tools", "benchmark-regression.R")
  testthat::skip_if_not(
    file.exists(helper_path),
    "tools/benchmark-helpers.R not reachable (installed-pkg context)"
  )

  env <- new.env(parent = globalenv())
  sys.source(helper_path, envir = env)

  script_lines <- readLines(script_path, warn = FALSE)
  cutoff <- grep("^benchmark_regression_main <- function", script_lines)[[1L]] - 1L
  eval(parse(text = paste(script_lines[seq_len(cutoff)], collapse = "\n")), envir = env)

  env
}

test_that("benchmark summary separates availability support routing and performance states", {
  ctx <- benchmark_regression_context()

  rows <- data.frame(
    suite = c(rep("dense", 5L), "sparse"),
    op = c(rep("matmul", 5L), "block_lanczos"),
    size_label = c(rep("small", 5L), "medium"),
    variant = c(rep("cold", 5L), "warm"),
    requested_backend = c("cpu", "opencl", "arrayfire", "metal", "opencl", "opencl"),
    dispatch_probe_op = c(rep("matmul", 5L), "block_lanczos"),
    requested_supported = c(NA, TRUE, TRUE, FALSE, FALSE, TRUE),
    requested_support_reason = c(NA, "cold supported", "cold supported", "op unsupported", "backend unavailable", "iterative supported"),
    dispatch_backend = c("cpu", "opencl", "cpu", "cpu", NA, "opencl"),
    dispatch_path = c("cold", "cold", "cold", "cold", NA, "iterative"),
    status = c("ok", "ok", "ok", "unsupported", "unavailable", "ok"),
    error_message = rep(NA_character_, 6L),
    nrow = c(rep(256L, 5L), 4000L),
    ncol = c(rep(32L, 5L), 1000L),
    rhs_width = c(rep(0L, 5L), 8L),
    nnz = c(rep(0L, 5L), 2000L),
    density = c(rep(0, 5L), 0.05),
    density_bucket = c(rep("dense", 5L), "sparse"),
    reps = rep(1L, 6L),
    median_ms = c(10, 5, 12, NA, NA, 4),
    cpu_reference_ms = c(rep(10, 5L), 8),
    baseline_ms = rep(NA_real_, 6L),
    ratio_vs_baseline = rep(NA_real_, 6L),
    stringsAsFactors = FALSE
  )

  summary <- ctx$summarize_results(rows)
  key <- paste(
    summary$requested_backend,
    summary$op,
    summary$status,
    ifelse(is.na(summary$dispatch_backend), "NA", summary$dispatch_backend),
    sep = "|"
  )
  summary <- summary[match(
    c(
      "cpu|matmul|ok|cpu",
      "opencl|matmul|ok|opencl",
      "arrayfire|matmul|ok|cpu",
      "metal|matmul|unsupported|cpu",
      "opencl|matmul|unavailable|NA",
      "opencl|block_lanczos|ok|opencl"
    ),
    key
  ), ]

  expect_identical(summary$selected_backend, summary$dispatch_backend)
  expect_identical(summary$availability_state, c("cpu_baseline", "available", "available", "available", "unavailable", "available"))
  expect_identical(summary$support_state, c("cpu_baseline", "supported", "supported", "unsupported", "backend_unavailable", "supported"))
  expect_identical(summary$execution_state, c("executed", "executed", "executed", "unsupported", "unavailable", "executed"))
  expect_identical(summary$routing_state, c("cpu_baseline", "accelerated", "cpu_fallback", "not_run", "not_run", "accelerated"))
  expect_identical(summary$dispatch_state, summary$routing_state)
  expect_identical(
    summary$performance_state,
    c("cpu_baseline", "accelerated_faster_than_cpu", "not_accelerated", "not_run", "not_run", "accelerated_faster_than_cpu")
  )
})

test_that("canonical worker backend gating enables MLX by default on Apple Silicon and respects overrides", {
  ctx <- benchmark_regression_context()

  withr::local_envvar(c(AMATRIX_BENCHMARK_MLX_WORKERS = NA_character_))
  expect_identical(
    ctx$benchmark_worker_backend_allowed("mlx", include_arrayfire = FALSE),
    ctx$benchmark_running_on_apple_silicon()
  )

  withr::local_envvar(c(AMATRIX_BENCHMARK_MLX_WORKERS = "1"))
  expect_true(ctx$benchmark_worker_backend_allowed("mlx", include_arrayfire = FALSE))
  expect_true(is.na(ctx$benchmark_worker_backend_block_reason("mlx", include_arrayfire = FALSE)))

  withr::local_envvar(c(AMATRIX_BENCHMARK_MLX_WORKERS = "0"))
  expect_false(ctx$benchmark_worker_backend_allowed("mlx", include_arrayfire = FALSE))
  if (ctx$benchmark_running_on_apple_silicon()) {
    expect_match(
      ctx$benchmark_worker_backend_block_reason("mlx", include_arrayfire = FALSE),
      "MLX disabled via AMATRIX_BENCHMARK_MLX_WORKERS=0"
    )
  } else {
    expect_match(
      ctx$benchmark_worker_backend_block_reason("mlx", include_arrayfire = FALSE),
      "MLX canonical workers require Apple Silicon or AMATRIX_BENCHMARK_MLX_WORKERS=1"
    )
  }
})

test_that("canonical worker backend gating hard-gates ArrayFire on Apple unless unsafe override is set", {
  ctx <- benchmark_regression_context()

  withr::local_envvar(c(
    AMATRIX_BENCHMARK_ARRAYFIRE = "1",
    AMATRIX_ARRAYFIRE_PROBE_GPU = NA_character_,
    AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE = NA_character_
  ))

  if (ctx$benchmark_running_on_apple_silicon()) {
    expect_false(ctx$benchmark_worker_backend_allowed("arrayfire", include_arrayfire = TRUE))
    expect_match(
      ctx$benchmark_worker_backend_block_reason("arrayfire", include_arrayfire = TRUE),
      "ArrayFire is diagnostic-only on Apple benchmark workers"
    )

    withr::local_envvar(c(AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE = "1"))
    expect_true(ctx$benchmark_worker_backend_allowed("arrayfire", include_arrayfire = TRUE))
  } else {
    expect_true(ctx$benchmark_worker_backend_allowed("arrayfire", include_arrayfire = TRUE))
  }
})

test_that("Apple ArrayFire report bundles switch to the direct CLI entry path and surface policy notes", {
  ctx <- benchmark_regression_context()

  rows <- data.frame(
    suite = c("dense", "dense"),
    op = c("matmul", "matmul"),
    size_label = c("small", "small"),
    variant = c("cold", "warm"),
    requested_backend = c("arrayfire", "arrayfire"),
    dispatch_probe_op = c("matmul", "matmul"),
    requested_supported = c(TRUE, TRUE),
    requested_support_reason = c("cold supported", "warm supported"),
    dispatch_backend = c("arrayfire", "arrayfire"),
    dispatch_path = c("cold", "warm"),
    status = c("ok", "ok"),
    error_message = c(NA_character_, NA_character_),
    nrow = c(128L, 128L),
    ncol = c(128L, 128L),
    rhs_width = c(128L, 128L),
    nnz = c(0L, 0L),
    density = c(0, 0),
    density_bucket = c("dense", "dense"),
    reps = c(3L, 3L),
    median_ms = c(4, 2),
    sd_ms = c(0.2, 0.1),
    cpu_reference_ms = c(8, 4),
    rel_err = c(1e-5, 1e-5),
    baseline_ms = c(5, 2.5),
    ratio_vs_baseline = c(0.8, 0.8),
    stringsAsFactors = FALSE
  )

  summary <- ctx$summarize_results(rows)
  metadata <- list(
    created_at = Sys.time(),
    hostname = "apple-host",
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    platform = "aarch64-apple-darwin20",
    sysname = "Darwin",
    arch = "arm64",
    timezone = "America/Toronto"
  )

  bundle <- ctx$benchmark_report_bundle(summary, "tools/baseline.csv", "tools/benchmark-results/test-arrayfire", metadata)

  expect_true("policy_notes" %in% names(bundle$tables))
  expect_match(bundle$snippets$run_bundle, "benchmark-regression-cli\\.R")
  expect_match(bundle$snippets$arrayfire_debug, "debug-arrayfire-gpu\\.R")
  expect_true(any(bundle$tables$policy_notes$area == "Apple ArrayFire"))
  expect_true(any(grepl("tools/benchmark-kmeans\\.R", bundle$tables$policy_notes$detail, fixed = FALSE)))
})
