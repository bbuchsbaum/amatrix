benchmark_helper_context <- function() {
  repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/", mustWork = TRUE)
  helper_path <- file.path(repo_root, "tools", "benchmark-helpers.R")
  helper_env <- new.env(parent = globalenv())
  sys.source(helper_path, envir = helper_env)
  list(repo_root = repo_root, helper_env = helper_env)
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
