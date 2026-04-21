# Repro metadata for amatrix-vmc:
# - Seeds: dense matmul worker uses X seed = 1 and RHS seed = 3 in tools/benchmark-regression.R
# - Dimensions: 256x32 %*% 32x32, 1024x128 %*% 128x128, 4096x128 %*% 128x128, 4096x1024 %*% 1024x1024
# - Backend / precision: arrayfire, precision = fast
# - Dispatch path: benchmark regression worker, dense matmul group, cold + warm variants
# - Failure before fix: worker aborted with `cl::Error: clGetDeviceIDs` before rows were classified
# - R / platform: R 4.5.1, aarch64-apple-darwin20, macOS Sonoma 14.3
# - Issue: amatrix-vmc

arrayfire_regression_context <- function() {
  repo_root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    winslash = "/", mustWork = FALSE
  )
  helper_path <- file.path(repo_root, "tools", "benchmark-helpers.R")
  testthat::skip_if_not(
    file.exists(helper_path),
    "tools/benchmark-helpers.R not reachable (installed-pkg context)"
  )
  helper_env <- new.env(parent = globalenv())
  sys.source(helper_path, envir = helper_env)
  list(
    repo_root = repo_root,
    script_path = file.path(repo_root, "tools", "benchmark-regression.R"),
    helper_env = helper_env
  )
}

test_that("amatrix-vmc: ArrayFire defaults to the cpu runtime on Apple Silicon", {
  testthat::skip_if_not(
    identical(Sys.info()[["sysname"]], "Darwin") &&
      grepl("arm64|aarch64", R.version$arch, ignore.case = TRUE),
    "Apple Silicon specific regression"
  )

  ctx <- arrayfire_regression_context()
  helper_env <- ctx$helper_env
  spec <- helper_env$.benchmark_optional_backend_specs(include_arrayfire = TRUE)[["arrayfire"]]
  testthat::skip_if(is.null(spec), "ArrayFire backend spec not available")
  launch <- helper_env$benchmark_system2_capture(
    file.path(R.home("bin"), "Rscript"),
    c(
      "-e",
      paste0(
        "setwd(", shQuote(ctx$repo_root), ");",
        "Sys.setenv(AMATRIX_ARRAYFIRE_PROBE_GPU='1');",
        "Sys.unsetenv('AMATRIX_ARRAYFIRE_BACKEND');",
        "options(amatrix.arrayfire.backend=NULL);",
        "source('tools/benchmark-helpers.R', local=globalenv());",
        "spec <- .benchmark_optional_backend_specs(include_arrayfire=TRUE)[['arrayfire']];",
        "stopifnot(.benchmark_enable_backend(spec));",
        "ns <- ensure_optional_backend_namespace('amatrix.arrayfire', repo_dir='backends/amatrix.arrayfire');",
        "diag <- get('amatrix_arrayfire_diagnostics', envir=ns, inherits=FALSE)();",
        "cat(sprintf('ACTIVE_BACKEND=%s\\n', diag$active_backend))"
      )
    )
  )

  expect_equal(launch$status, 0L, info = paste(launch$output, collapse = "\n"))
  expect_true(any(grepl("^ACTIVE_BACKEND=1$", launch$output)), info = paste(launch$output, collapse = "\n"))
})

test_that("amatrix-vmc: ArrayFire benchmark worker no longer crashes on Apple Silicon", {
  testthat::skip_if_not(
    identical(Sys.info()[["sysname"]], "Darwin") &&
      grepl("arm64|aarch64", R.version$arch, ignore.case = TRUE),
    "Apple Silicon specific regression"
  )

  ctx <- arrayfire_regression_context()
  repo_root <- ctx$repo_root
  helper_env <- ctx$helper_env
  spec <- helper_env$.benchmark_optional_backend_specs(include_arrayfire = TRUE)[["arrayfire"]]
  testthat::skip_if(is.null(spec), "ArrayFire backend spec not available")
  testthat::skip_if_not(dir.exists(file.path(repo_root, spec$repo_dir)), "ArrayFire backend source directory not found")

  old_benchmark <- Sys.getenv("AMATRIX_BENCHMARK_ARRAYFIRE", unset = "")
  old_unsafe <- Sys.getenv("AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE", unset = "")
  old_probe <- Sys.getenv("AMATRIX_ARRAYFIRE_PROBE_GPU", unset = "")
  old_backend <- Sys.getenv("AMATRIX_ARRAYFIRE_BACKEND", unset = "")
  on.exit({
    if (nzchar(old_benchmark)) Sys.setenv(AMATRIX_BENCHMARK_ARRAYFIRE = old_benchmark) else Sys.unsetenv("AMATRIX_BENCHMARK_ARRAYFIRE")
    if (nzchar(old_unsafe)) Sys.setenv(AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE = old_unsafe) else Sys.unsetenv("AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE")
    if (nzchar(old_probe)) Sys.setenv(AMATRIX_ARRAYFIRE_PROBE_GPU = old_probe) else Sys.unsetenv("AMATRIX_ARRAYFIRE_PROBE_GPU")
    if (nzchar(old_backend)) Sys.setenv(AMATRIX_ARRAYFIRE_BACKEND = old_backend) else Sys.unsetenv("AMATRIX_ARRAYFIRE_BACKEND")
  }, add = TRUE)

  Sys.setenv(
    AMATRIX_BENCHMARK_ARRAYFIRE = "1",
    AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE = "1",
    AMATRIX_ARRAYFIRE_PROBE_GPU = "1"
  )
  Sys.unsetenv("AMATRIX_ARRAYFIRE_BACKEND")

  plan <- list(list(
    group_id = "dense-arrayfire-matmul",
    suite = "dense",
    requested_backend = "arrayfire",
    op = "matmul"
  ))
  plan_path <- tempfile("benchmark-regression-af-crash-plan-", fileext = ".rds")
  out_path <- tempfile("benchmark-regression-af-crash-out-", fileext = ".rds")
  on.exit(unlink(c(plan_path, out_path)), add = TRUE)
  saveRDS(plan, plan_path)

  launch <- helper_env$benchmark_system2_capture(
    file.path(R.home("bin"), "Rscript"),
    helper_env$benchmark_rscript_source_args(
      ctx$script_path,
      working_dir = repo_root,
      main_call = "benchmark_regression_main(commandArgs(trailingOnly = TRUE))",
      args = c(
        "--worker",
        paste0("--plan=", plan_path),
        "--group-id=dense-arrayfire-matmul",
        paste0("--out=", out_path)
      )
    )
  )

  expect_equal(launch$status, 0L, info = paste(launch$output, collapse = "\n"))
  expect_true(file.exists(out_path), info = paste(launch$output, collapse = "\n"))

  rows <- readRDS(out_path)
  expect_true(all(rows$requested_backend == "arrayfire"))
  expect_false(any(rows$status == "crash"), info = paste(unique(stats::na.omit(rows$error_message)), collapse = "\n"))
  expect_true(any(rows$status == "ok"), info = paste(unique(stats::na.omit(rows$error_message)), collapse = "\n"))
})

test_that("amatrix-vmc: explicit cpu,arrayfire discovery avoids unrelated backend startup", {
  testthat::skip_if_not(
    identical(Sys.info()[["sysname"]], "Darwin") &&
      grepl("arm64|aarch64", R.version$arch, ignore.case = TRUE),
    "Apple Silicon specific regression"
  )

  ctx <- arrayfire_regression_context()
  repo_root <- ctx$repo_root
  helper_env <- ctx$helper_env
  spec <- helper_env$.benchmark_optional_backend_specs(include_arrayfire = TRUE)[["arrayfire"]]
  testthat::skip_if(is.null(spec), "ArrayFire backend spec not available")
  testthat::skip_if_not(dir.exists(file.path(repo_root, spec$repo_dir)), "ArrayFire backend source directory not found")

  old_benchmark <- Sys.getenv("AMATRIX_BENCHMARK_ARRAYFIRE", unset = "")
  old_unsafe <- Sys.getenv("AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE", unset = "")
  old_probe <- Sys.getenv("AMATRIX_ARRAYFIRE_PROBE_GPU", unset = "")
  old_backend <- Sys.getenv("AMATRIX_ARRAYFIRE_BACKEND", unset = "")
  on.exit({
    if (nzchar(old_benchmark)) Sys.setenv(AMATRIX_BENCHMARK_ARRAYFIRE = old_benchmark) else Sys.unsetenv("AMATRIX_BENCHMARK_ARRAYFIRE")
    if (nzchar(old_unsafe)) Sys.setenv(AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE = old_unsafe) else Sys.unsetenv("AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE")
    if (nzchar(old_probe)) Sys.setenv(AMATRIX_ARRAYFIRE_PROBE_GPU = old_probe) else Sys.unsetenv("AMATRIX_ARRAYFIRE_PROBE_GPU")
    if (nzchar(old_backend)) Sys.setenv(AMATRIX_ARRAYFIRE_BACKEND = old_backend) else Sys.unsetenv("AMATRIX_ARRAYFIRE_BACKEND")
  }, add = TRUE)

  Sys.setenv(
    AMATRIX_BENCHMARK_ARRAYFIRE = "1",
    AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE = "1",
    AMATRIX_ARRAYFIRE_PROBE_GPU = "1"
  )
  Sys.unsetenv("AMATRIX_ARRAYFIRE_BACKEND")

  launch <- helper_env$benchmark_system2_capture(
    file.path(R.home("bin"), "Rscript"),
    c(
      "-e",
      paste0(
        "setwd(", shQuote(repo_root), ");",
        "Sys.setenv(AMATRIX_BENCHMARK_ARRAYFIRE='1', AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE='1', AMATRIX_ARRAYFIRE_PROBE_GPU='1');",
        "Sys.unsetenv('AMATRIX_ARRAYFIRE_BACKEND');",
        "source('tools/benchmark-regression.R', local=globalenv());",
        "initialize_regression_benchmark_context();",
        "specs <- canonical_backend_specs(include_arrayfire=TRUE, only=c('cpu','arrayfire'));",
        "cat(sprintf('DISCOVERED=%s\\n', paste(names(specs), collapse=',')))"
      )
    )
  )

  expect_equal(launch$status, 0L, info = paste(launch$output, collapse = "\n"))
  expect_true(any(grepl("^DISCOVERED=cpu,arrayfire$", launch$output)), info = paste(launch$output, collapse = "\n"))
})
