benchmark_helper_context <- function() {
  repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/", mustWork = TRUE)
  helper_path <- file.path(repo_root, "tools", "benchmark-helpers.R")
  helper_env <- new.env(parent = globalenv())
  sys.source(helper_path, envir = helper_env)
  list(
    repo_root = repo_root,
    script_path = file.path(repo_root, "tools", "benchmark-regression.R"),
    helper_env = helper_env
  )
}

benchmark_summary_classification <- function(path) {
  summary <- read.csv(path, stringsAsFactors = FALSE)
  keys <- c(
    "suite", "op", "size_label", "variant", "requested_backend",
    "status", "selected_backend", "availability_state",
    "support_state", "execution_state", "routing_state"
  )
  summary <- summary[, keys, drop = FALSE]
  summary[do.call(order, summary[c("suite", "op", "size_label", "variant", "requested_backend")]), , drop = FALSE]
}

run_benchmark_regression_entry <- function(ctx, mode, output_dir, args) {
  helper_env <- ctx$helper_env
  script_path <- ctx$script_path
  rscript <- file.path(R.home("bin"), "Rscript")

  cli_args <- c(args, paste0("--output-dir=", output_dir))
  if (identical(mode, "direct-file")) {
    return(helper_env$benchmark_system2_capture(rscript, c(script_path, cli_args)))
  }

  if (identical(mode, "source")) {
    return(helper_env$benchmark_system2_capture(
      rscript,
      helper_env$benchmark_rscript_source_args(
        script_path,
        working_dir = ctx$repo_root,
        main_call = "benchmark_regression_main(commandArgs(trailingOnly = TRUE))",
        args = cli_args
      )
    ))
  }

  stop(sprintf("Unknown benchmark entry mode: %s", mode), call. = FALSE)
}

test_that("benchmark regression worker activates OpenCL groups without blanket unavailability", {
  ctx <- benchmark_helper_context()
  repo_root <- ctx$repo_root
  helper_env <- ctx$helper_env

  plan <- list(list(
    group_id = "dense-opencl-rsvd",
    suite = "dense",
    requested_backend = "opencl",
    op = "rsvd"
  ))
  plan_path <- tempfile("benchmark-regression-opencl-plan-", fileext = ".rds")
  out_path <- tempfile("benchmark-regression-opencl-out-", fileext = ".rds")
  on.exit(unlink(c(plan_path, out_path)), add = TRUE)
  saveRDS(plan, plan_path)

  script_path <- file.path(repo_root, "tools", "benchmark-regression.R")
  launch <- helper_env$benchmark_system2_capture(
    file.path(R.home("bin"), "Rscript"),
    helper_env$benchmark_rscript_source_args(
      script_path,
      working_dir = repo_root,
      main_call = "benchmark_regression_main(commandArgs(trailingOnly = TRUE))",
      args = c(
        "--worker",
        paste0("--plan=", plan_path),
        "--group-id=dense-opencl-rsvd",
        paste0("--out=", out_path)
      )
    )
  )

  expect_equal(launch$status, 0L, info = paste(launch$output, collapse = "\n"))
  expect_true(file.exists(out_path), info = paste(launch$output, collapse = "\n"))

  rows <- readRDS(out_path)
  expect_false(
    all(rows$status == "unavailable"),
    info = paste(unique(stats::na.omit(rows$error_message)), collapse = "\n")
  )
  expect_true(
    any(rows$status == "ok"),
    info = paste(unique(stats::na.omit(rows$error_message)), collapse = "\n")
  )
})

test_that("direct-file and sourced benchmark entrypoints classify OpenCL dense rows identically", {
  ctx <- benchmark_helper_context()

  direct_dir <- tempfile("benchmark-regression-direct-opencl-")
  source_dir <- tempfile("benchmark-regression-source-opencl-")
  unlink(c(direct_dir, source_dir), recursive = TRUE, force = TRUE)
  on.exit(unlink(c(direct_dir, source_dir), recursive = TRUE, force = TRUE), add = TRUE)

  args <- c("--backends=opencl", "--suites=dense")

  source_launch <- run_benchmark_regression_entry(ctx, "source", source_dir, args)
  expect_equal(source_launch$status, 0L, info = paste(source_launch$output, collapse = "\n"))

  source_summary_path <- file.path(source_dir, "summary.csv")
  expect_true(file.exists(source_summary_path), info = paste(source_launch$output, collapse = "\n"))

  source_summary <- benchmark_summary_classification(source_summary_path)
  if (!any(source_summary$status == "ok")) {
    skip("OpenCL benchmark entry is unavailable in this environment")
  }

  direct_launch <- run_benchmark_regression_entry(ctx, "direct-file", direct_dir, args)
  expect_equal(direct_launch$status, 0L, info = paste(direct_launch$output, collapse = "\n"))

  direct_summary_path <- file.path(direct_dir, "summary.csv")
  expect_true(file.exists(direct_summary_path), info = paste(direct_launch$output, collapse = "\n"))
  direct_summary <- benchmark_summary_classification(direct_summary_path)

  expect_identical(direct_summary, source_summary)
})

test_that("benchmark regression worker marks dense Metal groups unsupported instead of timing CPU fallback", {
  ctx <- benchmark_helper_context()
  repo_root <- ctx$repo_root
  helper_env <- ctx$helper_env

  metal_spec <- helper_env$.benchmark_optional_backend_specs(include_arrayfire = FALSE)[["metal"]]
  if (!isTRUE(helper_env$.benchmark_enable_backend(metal_spec))) {
    skip("Metal benchmark helper is unavailable in this environment")
  }

  plan <- list(list(
    group_id = "dense-metal-matmul",
    suite = "dense",
    requested_backend = "metal",
    op = "matmul"
  ))
  plan_path <- tempfile("benchmark-regression-metal-plan-", fileext = ".rds")
  out_path <- tempfile("benchmark-regression-metal-out-", fileext = ".rds")
  on.exit(unlink(c(plan_path, out_path)), add = TRUE)
  saveRDS(plan, plan_path)

  script_path <- file.path(repo_root, "tools", "benchmark-regression.R")
  launch <- helper_env$benchmark_system2_capture(
    file.path(R.home("bin"), "Rscript"),
    helper_env$benchmark_rscript_source_args(
      script_path,
      working_dir = repo_root,
      main_call = "benchmark_regression_main(commandArgs(trailingOnly = TRUE))",
      args = c(
        "--worker",
        paste0("--plan=", plan_path),
        "--group-id=dense-metal-matmul",
        paste0("--out=", out_path)
      )
    )
  )

  expect_equal(launch$status, 0L, info = paste(launch$output, collapse = "\n"))
  expect_true(file.exists(out_path), info = paste(launch$output, collapse = "\n"))

  rows <- readRDS(out_path)
  expect_true(all(rows$status == "unsupported"))
  expect_true(all(rows$requested_backend == "metal"))
  expect_true(all(rows$error_message %in% c("op unsupported", "calibration rejected", "resident op unsupported")))
})
