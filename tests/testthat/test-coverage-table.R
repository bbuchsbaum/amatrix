# Track 1 PR gate: enforces that every exported operation appears in the
# coverage matrix (§4) of planning_docs/quality-tracking.md.
#
# Always enforced: row existence for every non-excluded export.
#
# Strict mode (AMATRIX_COVERAGE_STRICT=1): additionally enforces that no row
# has '—' in Oracle / Metamorphic / Adversarial / Regression columns. Track 3
# is the work to close those gaps; Track 2's nightly gate sets the env var.
#
# Runs from the source tree only; skipped when planning_docs/ is not shipped
# (installed package, CRAN).

local({
  excluded_patterns <- c(
    "^amatrix_",
    "^with_amatrix$",
    "^adgeMatrix$", "^adgCMatrix$",
    "^as_adgeMatrix$", "^as_adgCMatrix$",
    "^as\\.",
    "^ncol\\.", "^nrow\\.",
    "^resident_handle$",
    "^rh_rowSums$", "^rh_colSums$",
    "^kron_matrix$"
  )

  is_excluded <- function(name) {
    any(vapply(
      excluded_patterns,
      function(p) grepl(p, name, perl = TRUE),
      logical(1)
    ))
  }

  find_pkg_root <- function() {
    root <- tryCatch(
      testthat::test_path("..", ".."),
      error = function(e) NA_character_
    )
    if (is.na(root) || !dir.exists(root)) NA_character_ else root
  }

  parse_namespace_exports <- function(namespace_path) {
    lines <- readLines(namespace_path, warn = FALSE)
    export_lines <- grep("^export\\([^)]+\\)$", lines, value = TRUE)
    sub("^export\\(([^)]+)\\)$", "\\1", export_lines)
  }

  read_matrix_region <- function(doc_path) {
    lines <- readLines(doc_path, warn = FALSE, encoding = "UTF-8")
    start_idx <- grep("^## 4\\. Coverage Matrix", lines)
    end_idx <- grep("^## 5\\. ", lines)
    if (length(start_idx) != 1 || length(end_idx) != 1) return(NULL)
    lines[seq(start_idx[1], end_idx[1] - 1)]
  }

  extract_backticked <- function(line) {
    m <- gregexpr("`([^`]+)`", line, perl = TRUE)[[1]]
    if (m[1] == -1L) return(character(0))
    starts <- as.integer(m)
    lens <- attr(m, "match.length")
    vapply(
      seq_along(starts),
      function(i) substr(line, starts[i] + 1L, starts[i] + lens[i] - 2L),
      character(1)
    )
  }

  test_that("every exported op has a row in the coverage matrix", {
    skip_on_cran()

    pkg_root <- find_pkg_root()
    skip_if(is.na(pkg_root), "coverage-table test requires source tree")

    namespace_path <- file.path(pkg_root, "NAMESPACE")
    doc_path <- file.path(pkg_root, "planning_docs", "quality-tracking.md")

    skip_if_not(file.exists(namespace_path), "NAMESPACE not found")
    skip_if_not(
      file.exists(doc_path),
      "planning_docs/quality-tracking.md not found (source tree only)"
    )

    exports <- parse_namespace_exports(namespace_path)
    ops <- exports[!vapply(exports, is_excluded, logical(1))]

    region <- read_matrix_region(doc_path)
    expect_false(
      is.null(region),
      info = "Could not find '## 4. Coverage Matrix' section in quality-tracking.md"
    )

    row_lines <- grep("^\\|\\s*`", region, value = TRUE)
    covered_names <- unique(unlist(lapply(row_lines, extract_backticked)))

    missing_rows <- setdiff(ops, covered_names)

    expect_true(
      length(missing_rows) == 0L,
      info = paste0(
        "Exports missing from planning_docs/quality-tracking.md §4:\n  - ",
        paste(missing_rows, collapse = "\n  - "),
        "\n\nAdd a row with test-type cells marked '\u2014' (em dash). New ",
        "exports must be classified as either an operation (add to \u00a74.2.x) ",
        "or an excluded API (update \u00a74.1 and the excluded_patterns list ",
        "in tests/testthat/test-coverage-table.R)."
      )
    )
  })

  test_that("coverage matrix has no gaps (strict mode)", {
    skip_on_cran()
    skip_if_not(
      identical(Sys.getenv("AMATRIX_COVERAGE_STRICT", "0"), "1"),
      "strict coverage check disabled (set AMATRIX_COVERAGE_STRICT=1 to enforce)"
    )

    pkg_root <- find_pkg_root()
    skip_if(is.na(pkg_root), "coverage-table test requires source tree")

    doc_path <- file.path(pkg_root, "planning_docs", "quality-tracking.md")
    skip_if_not(file.exists(doc_path), "quality doc not found")

    region <- read_matrix_region(doc_path)
    skip_if(is.null(region), "matrix region not found")

    row_lines <- grep("^\\|\\s*`", region, value = TRUE)
    gap_marker <- "\u2014"  # em dash U+2014
    gaps <- character(0)

    for (line in row_lines) {
      cells <- trimws(strsplit(line, "|", fixed = TRUE)[[1]])
      cells <- cells[nzchar(cells)]
      if (length(cells) < 6L) next
      op_label <- cells[1]
      test_cells <- cells[2:5]  # Oracle / Metamorphic / Adversarial / Regression
      if (any(grepl(gap_marker, test_cells, fixed = TRUE))) {
        gaps <- c(gaps, op_label)
      }
    }

    expect_true(
      length(gaps) == 0L,
      info = paste0(
        "Coverage matrix rows with gaps (strict mode). Fill Oracle / ",
        "Metamorphic / Adversarial / Regression columns to close these:\n  - ",
        paste(gaps, collapse = "\n  - ")
      )
    )
  })
})
