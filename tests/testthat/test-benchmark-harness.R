benchmark_harness_context <- function() {
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

make_synthetic_results <- function() {
  data.frame(
    suite             = c("dense", "dense", "dense", "dense", "dense"),
    op                = c("matmul", "crossprod", "rsvd", "matmul", "dist"),
    size_label        = c("medium", "medium", "medium", "small", "large"),
    variant           = c("cold",   "cold",      "cold",  "warm",   "cold"),
    requested_backend = "cpu",
    dispatch_backend  = "cpu",
    dispatch_path     = "cpu",
    nrow              = c(1024L, 1024L, 1024L, 256L, 4096L),
    ncol              = c(128L,  128L,  128L,  32L,  128L),
    rhs_width         = NA_integer_,
    density           = NA_real_,
    density_bucket    = "dense",
    status            = "ok",
    median_ms         = c(5.0, 9.0, 13.0, 0.9, 2.5),
    stringsAsFactors  = FALSE
  )
}

write_canonical_baseline_fixture <- function(path, ctx) {
  cols <- intersect(
    c(ctx$key_columns, "status", "median_ms"),
    names(make_synthetic_results())
  )
  rows <- make_synthetic_results()[, cols, drop = FALSE]
  utils::write.csv(rows, path, row.names = FALSE)
  path
}

write_legacy_baseline_fixture <- function(path) {
  rows <- data.frame(
    op          = c("matmul", "crossprod", "rsvd", "matmul", "dist"),
    size        = c("1024x128", "1024x128", "1024x128", "256x32", "4096x128"),
    backend     = "cpu",
    variant     = c("cold", "cold", "cold", "warm", "cold"),
    median_ms   = c(5.1, 9.1, 13.2, 0.95, 2.55),
    stringsAsFactors = FALSE
  )
  utils::write.csv(rows, path, row.names = FALSE)
  path
}

test_that("add_baseline_compare joins legacy baseline.csv without dropping rows", {
  ctx <- benchmark_harness_context()

  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  write_legacy_baseline_fixture(tmp)

  results <- make_synthetic_results()
  merged <- ctx$add_baseline_compare(results, tmp)

  expect_true("baseline_ms" %in% names(merged))
  expect_equal(nrow(merged), nrow(results))
  expect_true(all(!is.na(merged$baseline_ms)),
              info = "legacy baseline schema should fully join after normalization")
  expect_true(all(!is.na(merged$ratio_vs_baseline)))
})

test_that("add_baseline_compare joins canonical baseline.csv with every row matched", {
  ctx <- benchmark_harness_context()

  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  write_canonical_baseline_fixture(tmp, ctx)

  results <- make_synthetic_results()
  merged <- ctx$add_baseline_compare(results, tmp)

  expect_true(all(!is.na(merged$baseline_ms)))
  expect_equal(nrow(merged), nrow(results))
  expected_key <- with(results, paste(op, size_label, variant, requested_backend, sep = "|"))
  merged_key   <- with(merged,  paste(op, size_label, variant, requested_backend, sep = "|"))
  expect_setequal(merged_key, expected_key)
})

test_that("add_baseline_compare join is stable when nnz drifts between runs", {
  ctx <- benchmark_harness_context()

  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)

  baseline_rows <- data.frame(
    suite             = "sparse",
    op                = "spmm",
    size_label        = "medium",
    variant           = "cold",
    requested_backend = "cpu",
    dispatch_backend  = "cpu",
    dispatch_path     = "cpu",
    nrow              = 4000L,
    ncol              = 1000L,
    rhs_width         = 8L,
    density           = 0.01,
    density_bucket    = "sparse",
    status            = "ok",
    median_ms         = 12.3,
    nnz               = 40000L,
    stringsAsFactors  = FALSE
  )
  utils::write.csv(baseline_rows, tmp, row.names = FALSE)

  current <- baseline_rows
  current$median_ms <- 12.9
  current$nnz <- 40137L

  merged <- ctx$add_baseline_compare(current, tmp)
  expect_false("nnz" %in% ctx$key_columns,
               info = "nnz must not be part of the baseline join key")
  expect_true(all(!is.na(merged$baseline_ms)))
  expect_equal(merged$baseline_ms, 12.3)
})

test_that("tools/baseline.csv joins against every synthetic CPU result", {
  ctx <- benchmark_harness_context()
  repo_root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    winslash = "/", mustWork = FALSE
  )
  baseline_path <- file.path(repo_root, "tools", "baseline.csv")
  skip_if_not(file.exists(baseline_path), "tools/baseline.csv missing")

  results <- make_synthetic_results()
  merged <- ctx$add_baseline_compare(results, baseline_path)

  expect_true(all(!is.na(merged$baseline_ms)),
              info = "tools/baseline.csv must successfully join all synthetic rows")
})

test_that("write_outputs creates a readable benchmark report bundle", {
  ctx <- benchmark_harness_context()

  results <- data.frame(
    suite = c("dense", "dense", "dense", "dense"),
    op = c("matmul", "matmul", "crossprod", "crossprod"),
    size_label = c("small", "small", "medium", "medium"),
    variant = c("cold", "warm", "cold", "warm"),
    requested_backend = c("cpu", "cpu", "opencl", "opencl"),
    dispatch_probe_op = c("matmul", "matmul", "crossprod", "crossprod"),
    requested_supported = c(TRUE, TRUE, TRUE, TRUE),
    requested_support_reason = c("cold supported", "cold supported", "cold supported", "cold supported"),
    dispatch_backend = c("cpu", "cpu", "opencl", "opencl"),
    dispatch_path = c("cold", "cold", "cold", "cold"),
    status = c("ok", "ok", "ok", "ok"),
    error_message = NA_character_,
    nrow = c(256L, 256L, 1024L, 1024L),
    ncol = c(32L, 32L, 128L, 128L),
    rhs_width = c(32L, 32L, 0L, 0L),
    density = c(0, 0, 0, 0),
    density_bucket = c("dense", "dense", "dense", "dense"),
    nnz = c(0L, 0L, 0L, 0L),
    reps = c(7L, 7L, 7L, 7L),
    median_ms = c(2.0, 1.0, 5.0, 3.0),
    mean_ms = c(2.1, 1.1, 5.1, 3.2),
    sd_ms = c(0.2, 0.1, 0.4, 0.2),
    p05_ms = c(1.8, 0.9, 4.8, 2.9),
    p95_ms = c(2.3, 1.2, 5.5, 3.4),
    n_reps = c(7L, 7L, 7L, 7L),
    rel_err = c(NA_real_, NA_real_, 1e-5, 1e-5),
    cpu_reference_ms = c(2.0, 1.0, 6.0, 4.5),
    baseline_ms = c(1.7, 1.0, 4.0, 2.0),
    ratio_vs_baseline = c(2.0 / 1.7, 1.0 / 1.0, 5.0 / 4.0, 3.0 / 2.0),
    stringsAsFactors = FALSE
  )

  out_dir <- tempfile("benchmark-report-bundle-")
  unlink(out_dir, recursive = TRUE, force = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  baseline_path <- tempfile(fileext = ".csv")
  write.csv(results[, c("op", "size_label", "variant", "requested_backend", "median_ms")], baseline_path, row.names = FALSE)
  on.exit(unlink(baseline_path), add = TRUE)

  paths <- ctx$write_outputs(results, out_dir, baseline_path = baseline_path)

  expect_true(file.exists(paths$raw))
  expect_true(file.exists(paths$summary))
  expect_true(file.exists(paths$regressions))
  expect_true(file.exists(paths$warm_ratios))
  expect_true(file.exists(paths$routing_summary))
  expect_true(file.exists(paths$summary_md))
  expect_true(file.exists(paths$report_data))
  expect_true(file.exists(paths$report_qmd))
  expect_true(file.exists(paths$report_css))
  expect_true(file.exists(paths$report_tex))
  expect_true(dir.exists(paths$plots_dir))
  expect_true(file.exists(file.path(paths$plots_dir, "baseline-drift.png")))
  expect_true(file.exists(file.path(paths$plots_dir, "warm-cold-gains.png")))
  expect_true(file.exists(file.path(paths$plots_dir, "routing-overview.png")))

  report <- readRDS(paths$report_data)
  expect_true("policy_notes" %in% names(report$tables))
  expect_true("backend_overview" %in% names(report$tables))
  expect_true("suite_overview" %in% names(report$tables))
  expect_true("op_coverage" %in% names(report$tables))
  expect_gt(nrow(report$tables$backend_overview), 0L)
  expect_gt(nrow(report$tables$suite_overview), 0L)
  expect_gt(nrow(report$tables$op_coverage), 0L)
  expect_gt(nrow(report$tables$warm_pairs), 0L)
  expect_true("snippets" %in% names(report))
  expect_match(report$snippets$render_pdf, "benchmark-report\\.pdf")

  if (!is.na(paths$report_html)) {
    expect_true(file.exists(paths$report_html))
  }
})
