# Repro metadata:
# - Seed(s): none; fixed synthetic benchmark rows
# - Dimensions: sparse medium suite cells (4000 x 1000) with rhs_width 8 and 32
# - Backend / precision / dispatch: cpu / strict / cold + resident summary rows
# - R / platform: captured at runtime by sessionInfo()
# - Issue: benchmark report plotting crashed on duplicate display labels when
#   sparse drift or warm-state rows differed only by density / rhs_width

benchmark_plot_regression_context <- function() {
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

test_that("benchmark report plots tolerate duplicate sparse display labels [amatrix-drd]", {
  ctx <- benchmark_plot_regression_context()

  results <- data.frame(
    suite = rep("sparse", 4L),
    op = rep("spmm", 4L),
    size_label = rep("medium", 4L),
    variant = c("cold", "resident", "cold", "resident"),
    requested_backend = rep("cpu", 4L),
    dispatch_probe_op = rep("matmul", 4L),
    requested_supported = rep(TRUE, 4L),
    requested_support_reason = rep("cold supported", 4L),
    dispatch_backend = rep("cpu", 4L),
    dispatch_path = rep("cold", 4L),
    status = rep("ok", 4L),
    error_message = rep(NA_character_, 4L),
    nrow = rep(4000L, 4L),
    ncol = rep(1000L, 4L),
    rhs_width = c(8L, 8L, 32L, 32L),
    density = c(0.001, 0.001, 0.005, 0.005),
    density_bucket = rep("sparse", 4L),
    nnz = c(4000L, 4000L, 20000L, 20000L),
    reps = rep(5L, 4L),
    median_ms = c(12, 8, 16, 10),
    mean_ms = c(12, 8, 16, 10),
    sd_ms = c(0.5, 0.5, 0.8, 0.8),
    p05_ms = c(11.5, 7.5, 15.2, 9.3),
    p95_ms = c(12.6, 8.4, 16.9, 10.7),
    n_reps = rep(5L, 4L),
    rel_err = rep(NA_real_, 4L),
    cpu_reference_ms = c(12, 8, 16, 10),
    baseline_ms = c(10, 7, 12, 9),
    ratio_vs_baseline = c(1.2, 8 / 7, 16 / 12, 10 / 9),
    stringsAsFactors = FALSE
  )

  out_dir <- tempfile("benchmark-plot-duplicate-labels-")
  unlink(out_dir, recursive = TRUE, force = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  baseline_path <- tempfile(fileext = ".csv")
  utils::write.csv(
    results[, c(
      "suite", "op", "size_label", "variant", "requested_backend",
      "dispatch_backend", "dispatch_path", "nrow", "ncol", "rhs_width",
      "density", "density_bucket", "status", "median_ms"
    )],
    baseline_path,
    row.names = FALSE
  )
  on.exit(unlink(baseline_path), add = TRUE)

  expect_no_error({
    paths <- ctx$write_outputs(results, out_dir, baseline_path = baseline_path)
    expect_true(file.exists(file.path(paths$plots_dir, "baseline-drift.png")))
    expect_true(file.exists(file.path(paths$plots_dir, "warm-cold-gains.png")))
  })
})
