benchmark_harness_context <- function() {
  repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/", mustWork = TRUE)
  helper_path <- file.path(repo_root, "tools", "benchmark-helpers.R")
  script_path <- file.path(repo_root, "tools", "benchmark-regression.R")

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
  repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), winslash = "/", mustWork = TRUE)
  baseline_path <- file.path(repo_root, "tools", "baseline.csv")
  skip_if_not(file.exists(baseline_path), "tools/baseline.csv missing")

  results <- make_synthetic_results()
  merged <- ctx$add_baseline_compare(results, baseline_path)

  expect_true(all(!is.na(merged$baseline_ms)),
              info = "tools/baseline.csv must successfully join all synthetic rows")
})
