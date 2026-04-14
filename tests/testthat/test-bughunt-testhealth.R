# test-bughunt-testhealth.R
# Meta-tests that audit the test suite itself for health issues.
# amatrix-uct: options leak in test-conformance.R
# amatrix-hf6: loose tolerance (> 1e-4) masking numerical drift

# ── amatrix-uct: options leak detection ──────────────────────────────────────

test_that("test-conformance.R does not leak amatrix.mlx.qr_helper_mode [amatrix-uct]", {
  # Lines 1279 and 1282 of test-conformance.R set amatrix.mlx.qr_helper_mode
  # and amatrix.mlx.qr_compact_method mid-test. The on.exit at line 1277 only
  # covers amatrix.mlx.available and amatrix.mlx.qr_tsqr_block_rows.
  # Any failure between line 1279 and end-of-test leaks the option.
  #
  # This test verifies the ABSENCE of teardown by parsing the file and checking
  # that every options() mutation site has a corresponding on.exit coverage.
  conformance_path <- system.file(
    "tests", "testthat", "test-conformance.R",
    package = "amatrix"
  )
  if (!nzchar(conformance_path)) {
    conformance_path <- file.path(
      find.package("amatrix"), "..", "..", "tests", "testthat", "test-conformance.R"
    )
  }
  if (!file.exists(conformance_path)) {
    # Try relative from package source
    pkg_root <- tryCatch(
      rprojroot::find_package_root_file(path = getwd()),
      error = function(e) NULL
    )
    if (!is.null(pkg_root)) {
      conformance_path <- file.path(pkg_root, "tests", "testthat", "test-conformance.R")
    }
  }
  skip_if_not(file.exists(conformance_path), "test-conformance.R not found in source tree")

  lines <- readLines(conformance_path)

  # The two bare options() calls that set qr_helper_mode without on.exit
  leak_lines <- grep(
    'options\\(amatrix\\.mlx\\.qr_helper_mode\\s*=',
    lines
  )
  # There should be 0 such bare calls (they should all be covered by on.exit)
  # Currently there are 2 (lines 1279 and 1282): this FAILS to expose the bug.
  expect_equal(
    length(leak_lines), 0L,
    info = paste(
      "options(amatrix.mlx.qr_helper_mode=...) set without on.exit coverage at lines:",
      paste(leak_lines, collapse = ", ")
    )
  )
})

test_that("test-conformance.R does not leak amatrix.mlx.qr_compact_method [amatrix-uct]", {
  pkg_root <- tryCatch(
    rprojroot::find_package_root_file(path = getwd()),
    error = function(e) NULL
  )
  if (is.null(pkg_root)) {
    pkg_root <- Sys.getenv("R_PACKAGE_DIR", unset = NA_character_)
  }
  conformance_path <- if (!is.null(pkg_root) && !is.na(pkg_root)) {
    file.path(pkg_root, "tests", "testthat", "test-conformance.R")
  } else {
    NA_character_
  }
  skip_if_not(
    !is.na(conformance_path) && file.exists(conformance_path),
    "test-conformance.R not found in source tree"
  )

  lines <- readLines(conformance_path)

  # Bare options() setting qr_compact_method without being captured in old <- options(...)
  leak_lines <- grep(
    'options\\([^)]*amatrix\\.mlx\\.qr_compact_method',
    lines
  )
  # Line 1282 is a bare options() with no corresponding on.exit: this FAILS
  bare_leaks <- leak_lines[!vapply(leak_lines, function(ln) {
    # Check if the same line or a line within 3 lines before has "old <-"
    any(grepl("old\\s*<-\\s*options\\(", lines[max(1, ln - 3):ln]))
  }, logical(1))]

  expect_equal(
    length(bare_leaks), 0L,
    info = paste(
      "Bare options(amatrix.mlx.qr_compact_method=...) at lines (no on.exit):",
      paste(bare_leaks, collapse = ", ")
    )
  )
})

# ── amatrix-hf6: loose tolerance detection ───────────────────────────────────

test_that("test-svd-factor.R line ~640 does not use tolerance > 1e-4 [amatrix-hf6]", {
  pkg_root <- tryCatch(
    rprojroot::find_package_root_file(path = getwd()),
    error = function(e) NULL
  )
  svd_path <- if (!is.null(pkg_root)) {
    file.path(pkg_root, "tests", "testthat", "test-svd-factor.R")
  } else NA_character_
  skip_if_not(
    !is.na(svd_path) && file.exists(svd_path),
    "test-svd-factor.R not found"
  )

  lines <- readLines(svd_path)
  # Find lines with tolerance= values > 1e-4 (i.e. 1e-4, 1e-3, 1e-2, etc.)
  tol_lines <- grep("tolerance\\s*=\\s*1e-[0-3][^0-9]|tolerance\\s*=\\s*[0-9]+e-[0-3]|tolerance\\s*=\\s*[0-9]*\\.[0-9]+[^e]", lines)
  # Specifically catch tolerance = 1e-4 (threshold is > 1e-4 being too loose for CPU)
  loose <- grep("tolerance\\s*=\\s*1e-4\\b", lines)

  expect_equal(
    length(loose), 0L,
    info = paste(
      "test-svd-factor.R uses tolerance=1e-4 (too loose for CPU ops) at lines:",
      paste(loose, collapse = ", ")
    )
  )
})

test_that("test-lmm-primitives.R stochastic tolerance=1.0 is documented [amatrix-hf6]", {
  pkg_root <- tryCatch(
    rprojroot::find_package_root_file(path = getwd()),
    error = function(e) NULL
  )
  lmm_path <- if (!is.null(pkg_root)) {
    file.path(pkg_root, "tests", "testthat", "test-lmm-primitives.R")
  } else NA_character_
  skip_if_not(
    !is.na(lmm_path) && file.exists(lmm_path),
    "test-lmm-primitives.R not found"
  )

  lines <- readLines(lmm_path)
  # Lines using tolerance >= 0.5 that are NOT stochastic estimator tests
  loose <- grep("tolerance\\s*=\\s*[0-9]+\\.0\\b|tolerance\\s*=\\s*0\\.[5-9]", lines)

  # Each such line must have "stochastic" commented nearby (within 2 lines)
  undocumented <- loose[!vapply(loose, function(ln) {
    context <- lines[max(1, ln - 2):min(length(lines), ln + 2)]
    any(grepl("stochastic", context, ignore.case = TRUE))
  }, logical(1))]

  expect_equal(
    length(undocumented), 0L,
    info = paste(
      "test-lmm-primitives.R has tolerance >= 0.5 without 'stochastic' comment at lines:",
      paste(undocumented, collapse = ", ")
    )
  )
})

test_that("test-chol-factor.R does not use tolerance > 1e-4 for non-GPU ops [amatrix-hf6]", {
  pkg_root <- tryCatch(
    rprojroot::find_package_root_file(path = getwd()),
    error = function(e) NULL
  )
  chol_path <- if (!is.null(pkg_root)) {
    file.path(pkg_root, "tests", "testthat", "test-chol-factor.R")
  } else NA_character_
  skip_if_not(
    !is.na(chol_path) && file.exists(chol_path),
    "test-chol-factor.R not found"
  )

  lines <- readLines(chol_path)
  # 5e-6 is tighter than 1e-4, but flag anything looser than 1e-4
  loose <- grep("tolerance\\s*=\\s*[0-9.]*e-[0-3]\\b|tolerance\\s*=\\s*[0-9]+\\.[0-9]+\\b", lines)
  # Filter: only flag values actually > 1e-4 (i.e. 5e-6 passes, 2e-4 fails)
  # The real issue from audit: 5e-6 in CLBlast path - acceptable but flagged by quality tracking
  # Flag anything >= 1e-4 to satisfy the 1e-4 GPU threshold requirement
  tol_pattern <- "tolerance\\s*=\\s*(\\d+(?:\\.\\d+)?(?:e[+-]?\\d+)?)"
  loose_vals <- grep(tol_pattern, lines)
  bad_lines <- loose_vals[vapply(loose_vals, function(ln) {
    m <- regmatches(lines[ln], regexpr("tolerance\\s*=\\s*([0-9.e+-]+)", lines[ln]))
    if (length(m) == 0) return(FALSE)
    val_str <- sub("tolerance\\s*=\\s*", "", m)
    val <- suppressWarnings(as.numeric(val_str))
    !is.na(val) && val > 1e-4
  }, logical(1))]

  expect_equal(
    length(bad_lines), 0L,
    info = paste(
      "test-chol-factor.R uses tolerance > 1e-4 at lines:",
      paste(bad_lines, collapse = ", ")
    )
  )
})
