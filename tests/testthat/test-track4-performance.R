# Track 4 — performance under discipline
#
# Tests the calibration + benchmark-report surface that enforces the Track 4
# speed contract:
#
#   * .amatrix_sys_hash produces a deterministic, non-empty fingerprint
#   * Calibration cache is invalidated when the hash changes (simulated)
#   * amatrix_benchmark_report() reads tools/baseline.csv and produces a
#     cold/warm/ratio frame plus a calibration frame
#
# These tests stay in the source-tree-only regime (read tools/baseline.csv
# via a relative path) and are skipped when the file is absent.

test_that(".amatrix_sys_hash returns a stable non-empty string", {
  h1 <- amatrix:::.amatrix_sys_hash()
  h2 <- amatrix:::.amatrix_sys_hash()

  expect_type(h1, "character")
  expect_true(nzchar(h1))
  expect_identical(h1, h2)                    # deterministic
  expect_true(startsWith(h1, "v1|"))          # versioned
  expect_true(grepl(R.version$platform, h1, fixed = TRUE))
})

test_that("calibration cache is invalidated on hash mismatch", {
  skip_on_cran()

  st <- amatrix:::.amatrix_state
  ns <- asNamespace("amatrix")

  # Save any pre-existing session state and restore on exit.
  prior <- st$calibration
  on.exit(st$calibration <- prior, add = TRUE)

  tmp_cache <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_cache), add = TRUE)

  bogus_cal <- list(
    version       = "2",
    calibrated_at = Sys.time(),
    sys_hash      = "v1|FAKE|bogus|xx|fake-platform|99.99",
    thresholds    = list(mlx = list(gemm = 4096L)),
    results       = data.frame()
  )
  saveRDS(bogus_cal, tmp_cache)

  st$calibration <- NULL

  # Stub the path helper just for this test.
  orig_path_fn <- get(".amatrix_calibration_path", envir = ns)
  unlockBinding(".amatrix_calibration_path", ns)
  on.exit({
    assign(".amatrix_calibration_path", orig_path_fn, envir = ns)
    lockBinding(".amatrix_calibration_path", ns)
  }, add = TRUE)
  assign(".amatrix_calibration_path", function() tmp_cache, envir = ns)

  # Load should REJECT the cache because its sys_hash doesn't match.
  amatrix:::.amatrix_load_calibration()
  expect_null(st$calibration,
              info = "stale cache must be rejected on hash mismatch")

  # Rewrite the cache with the current hash — load should accept.
  good_cal <- bogus_cal
  good_cal$sys_hash <- amatrix:::.amatrix_sys_hash()
  saveRDS(good_cal, tmp_cache)
  st$calibration <- NULL
  amatrix:::.amatrix_load_calibration()
  expect_false(is.null(st$calibration),
               info = "cache with matching hash must be accepted")
  expect_identical(st$calibration$version, "2")
})

test_that("amatrix_benchmark_report reads baseline.csv and pivots cold/warm", {
  skip_on_cran()

  baseline_path <- testthat::test_path("..", "..", "tools", "baseline.csv")
  skip_if_not(file.exists(baseline_path), "tools/baseline.csv not shipped")

  rep <- amatrix_benchmark_report(baseline_path = baseline_path)

  expect_named(rep, c("baseline", "calibration"))
  expect_s3_class(rep$baseline, "data.frame")
  expect_s3_class(rep$calibration, "data.frame")
  expect_true(all(c("op", "size", "backend", "cold_ms", "warm_ms",
                    "warm_vs_cold_ratio", "speedup_vs_cpu") %in% names(rep$baseline)))
  expect_gt(nrow(rep$baseline), 0L)

  # cpu rows should all have speedup_vs_cpu == 1 (base reference).
  cpu_rows <- rep$baseline[rep$baseline$backend == "cpu", , drop = FALSE]
  expect_true(all(is.finite(cpu_rows$cold_ms) | is.finite(cpu_rows$warm_ms)))
  expect_true(all(cpu_rows$speedup_vs_cpu == 1, na.rm = TRUE))

  # warm_vs_cold_ratio is NA when either variant is missing, finite otherwise.
  have_both <- is.finite(rep$baseline$cold_ms) & is.finite(rep$baseline$warm_ms)
  if (any(have_both)) {
    expect_true(all(is.finite(rep$baseline$warm_vs_cold_ratio[have_both])))
  }
})

test_that("amatrix_benchmark_report gracefully handles a missing baseline", {
  rep <- amatrix_benchmark_report(baseline_path = tempfile(fileext = ".csv"))
  expect_named(rep, c("baseline", "calibration"))
  expect_identical(nrow(rep$baseline), 0L)
})

test_that("amatrix_benchmark_report gracefully handles baseline_path = NULL", {
  rep <- amatrix_benchmark_report(baseline_path = NULL)
  expect_named(rep, c("baseline", "calibration"))
  expect_identical(nrow(rep$baseline), 0L)
})
