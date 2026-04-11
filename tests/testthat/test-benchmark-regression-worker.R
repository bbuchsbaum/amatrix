opencl_worker_spec <- function() {
  specs <- optional_backend_specs()
  specs[[match("amatrix.opencl", vapply(specs, `[[`, character(1), "package"))]]
}

test_that("benchmark regression worker activates OpenCL groups without blanket unavailability", {
  spec <- opencl_worker_spec()
  skip_if_backend_package_missing(spec)

  repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/", mustWork = TRUE)
  helper_path <- file.path(repo_root, "tools", "benchmark-helpers.R")
  helper_env <- new.env(parent = globalenv())
  sys.source(helper_path, envir = helper_env)

  opencl_spec <- helper_env$.benchmark_optional_backend_specs(include_arrayfire = FALSE)[["opencl"]]
  if (!isTRUE(helper_env$.benchmark_enable_backend(opencl_spec))) {
    skip("OpenCL benchmark helper is unavailable in this environment")
  }

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
