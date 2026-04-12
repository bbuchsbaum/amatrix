benchmark_regression_context <- function() {
  repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/", mustWork = TRUE)
  helper_path <- file.path(repo_root, "tools", "benchmark-helpers.R")
  script_path <- file.path(repo_root, "tools", "benchmark-regression.R")

  env <- new.env(parent = globalenv())
  sys.source(helper_path, envir = env)

  script_lines <- readLines(script_path, warn = FALSE)
  cutoff <- grep("^args <- parse_args\\(", script_lines)[[1L]] - 1L
  eval(parse(text = paste(script_lines[seq_len(cutoff)], collapse = "\n")), envir = env)

  env
}

test_that("benchmark summary separates availability support routing and performance states", {
  ctx <- benchmark_regression_context()

  rows <- data.frame(
    suite = rep("dense", 5L),
    op = rep("matmul", 5L),
    size_label = rep("small", 5L),
    variant = rep("cold", 5L),
    requested_backend = c("cpu", "opencl", "arrayfire", "metal", "opencl"),
    dispatch_probe_op = rep("matmul", 5L),
    requested_supported = c(NA, TRUE, TRUE, FALSE, FALSE),
    requested_support_reason = c(NA, "cold supported", "cold supported", "op unsupported", "backend unavailable"),
    dispatch_backend = c("cpu", "opencl", "cpu", "cpu", NA),
    dispatch_path = c("cold", "cold", "cold", "cold", NA),
    status = c("ok", "ok", "ok", "unsupported", "unavailable"),
    error_message = rep(NA_character_, 5L),
    nrow = rep(256L, 5L),
    ncol = rep(32L, 5L),
    rhs_width = rep(0L, 5L),
    nnz = rep(0L, 5L),
    density = rep(0, 5L),
    density_bucket = rep("dense", 5L),
    reps = rep(1L, 5L),
    median_ms = c(10, 5, 12, NA, NA),
    cpu_reference_ms = rep(10, 5L),
    baseline_ms = rep(NA_real_, 5L),
    ratio_vs_baseline = rep(NA_real_, 5L),
    stringsAsFactors = FALSE
  )

  summary <- ctx$summarize_results(rows)
  key <- paste(summary$requested_backend, summary$status, ifelse(is.na(summary$dispatch_backend), "NA", summary$dispatch_backend), sep = "|")
  summary <- summary[match(
    c(
      "cpu|ok|cpu",
      "opencl|ok|opencl",
      "arrayfire|ok|cpu",
      "metal|unsupported|cpu",
      "opencl|unavailable|NA"
    ),
    key
  ), ]

  expect_identical(summary$selected_backend, summary$dispatch_backend)
  expect_identical(summary$availability_state, c("cpu_baseline", "available", "available", "available", "unavailable"))
  expect_identical(summary$support_state, c("cpu_baseline", "supported", "supported", "unsupported", "backend_unavailable"))
  expect_identical(summary$execution_state, c("executed", "executed", "executed", "unsupported", "unavailable"))
  expect_identical(summary$routing_state, c("cpu_baseline", "accelerated", "cpu_fallback", "not_run", "not_run"))
  expect_identical(summary$dispatch_state, summary$routing_state)
  expect_identical(
    summary$performance_state,
    c("cpu_baseline", "accelerated_faster_than_cpu", "not_accelerated", "not_run", "not_run")
  )
})

test_that("canonical worker backend gating disables MLX by default", {
  ctx <- benchmark_regression_context()

  withr::local_envvar(c(AMATRIX_BENCHMARK_MLX_WORKERS = NA_character_))
  expect_false(ctx$benchmark_worker_backend_allowed("mlx", include_arrayfire = FALSE))
  expect_match(
    ctx$benchmark_worker_backend_block_reason("mlx", include_arrayfire = FALSE),
    "MLX is disabled in canonical regression workers"
  )

  withr::local_envvar(c(AMATRIX_BENCHMARK_MLX_WORKERS = "1"))
  expect_true(ctx$benchmark_worker_backend_allowed("mlx", include_arrayfire = FALSE))
  expect_true(is.na(ctx$benchmark_worker_backend_block_reason("mlx", include_arrayfire = FALSE)))
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
      "ArrayFire is disabled on Apple benchmark workers"
    )

    withr::local_envvar(c(AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE = "1"))
    expect_true(ctx$benchmark_worker_backend_allowed("arrayfire", include_arrayfire = TRUE))
  } else {
    expect_true(ctx$benchmark_worker_backend_allowed("arrayfire", include_arrayfire = TRUE))
  }
})
