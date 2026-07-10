# Sync the coverage matrix (planning_docs/quality-tracking.md §4) against
# NAMESPACE exports and the benchmark harness op list (mote amatrix-zxn).
#
# Usage:
#   Rscript tools/sync-coverage-table.R          # check: exit 1 with a diff on drift
#   Rscript tools/sync-coverage-table.R --fix    # append stub rows / upgrade
#                                                # Benchmark cells, then re-check
#
# Checks:
#   1. Every non-excluded export has a row in the matrix (same contract the
#      PR-gate test tests/testthat/test-coverage-table.R enforces; this tool
#      exists to FIX drift, the test only detects it).
#   2. Every benchmark-harness op maps to a row whose Benchmark column is
#      populated (✓ or ○). With --fix, '—' cells for harness-measured ops are
#      upgraded to '✓'; '○' cells are reported but never auto-upgraded.
#   3. Simple single-token rows whose export no longer exists are reported as
#      orphans (never auto-deleted).

sync_coverage_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  fix <- "--fix" %in% args
  root <- normalizePath(".", mustWork = TRUE)
  ns_path <- file.path(root, "NAMESPACE")
  doc_path <- file.path(root, "planning_docs", "quality-tracking.md")
  stopifnot(file.exists(ns_path), file.exists(doc_path))

  # Keep in lockstep with tests/testthat/test-coverage-table.R.
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
    any(vapply(excluded_patterns, function(p) grepl(p, name, perl = TRUE), logical(1)))
  }

  # Harness op -> Export token in the matrix. Validated below so renames on
  # either side are caught instead of silently ignored. Keep in lockstep with
  # group_plan() / sparse_iterative_ops() in tools/benchmark-regression.R.
  harness_map <- c(
    matmul = "matmul", crossprod = "crossprod", covariance = "covariance",
    dist = "dist_matrix", chol = "chol", solve_rhs = "solve",
    eigen_sym = "eigh", many_lm = "many_lm", rsvd = "rsvd", svd = "svd",
    sinkhorn = "sinkhorn", spmv = "sparse %*%", spmm = "sparse %*%",
    block_lanczos = "block_lanczos", svd_factor_subspace = "svd_factor"
  )

  lines <- readLines(doc_path, warn = FALSE, encoding = "UTF-8")
  start_idx <- grep("^## 4\\. Coverage Matrix", lines)
  end_idx <- grep("^## 5\\. ", lines)
  stopifnot(length(start_idx) == 1L, length(end_idx) == 1L)
  region <- seq(start_idx, end_idx - 1L)

  row_idx <- region[grepl("^\\| `", lines[region])]
  first_token <- function(line) {
    m <- regmatches(line, regexpr("`[^`]+`", line))
    if (length(m) == 0L) NA_character_ else gsub("`", "", m)
  }
  row_tokens <- vapply(lines[row_idx], first_token, character(1), USE.NAMES = FALSE)

  ns_lines <- readLines(ns_path, warn = FALSE)
  exports <- sub("^export\\(([^)]+)\\)$", "\\1",
                 grep("^export\\([^)]+\\)$", ns_lines, value = TRUE))
  ops <- exports[!vapply(exports, is_excluded, logical(1))]

  problems <- character(0)

  # 1. Missing rows.
  missing <- setdiff(ops, row_tokens)
  if (length(missing) > 0L && fix) {
    stub <- sprintf("| `%s` | — | — | — | — | — | C? M? A? O? X— |", missing)
    header_needed <- !any(grepl("^#### 4\\.2\\.98 Auto-added", lines))
    insert_at <- end_idx - 1L
    block <- c(
      if (header_needed) c(
        "",
        "#### 4.2.98 Auto-added by tools/sync-coverage-table.R (needs triage)",
        "",
        "| Export | Oracle | Metamorphic | Adversarial | Regression | Benchmark | Tiers |",
        "|---|---|---|---|---|---|---|"
      ),
      stub
    )
    if (header_needed) {
      lines <- append(lines, block, after = insert_at)
    } else {
      auto_hdr <- grep("^#### 4\\.2\\.98 Auto-added", lines)
      auto_rows <- which(grepl("^\\| `", lines) & seq_along(lines) > auto_hdr)
      after <- if (length(auto_rows) > 0L) max(auto_rows) else auto_hdr + 2L
      lines <- append(lines, stub, after = after)
    }
    message(sprintf("FIXED: appended %d stub row(s): %s",
                    length(missing), paste(missing, collapse = ", ")))
    missing <- character(0)
    end_idx <- grep("^## 5\\. ", lines) # region shifted
  }
  if (length(missing) > 0L) {
    problems <- c(problems, sprintf(
      "Missing coverage row for exported op(s): %s (run with --fix to append stubs)",
      paste(missing, collapse = ", ")
    ))
  }

  # 2. Benchmark column vs harness ops.
  set_cell <- function(line, col, value) {
    parts <- strsplit(line, "|", fixed = TRUE)[[1]]
    parts[col + 1L] <- sprintf(" %s ", value)
    paste(parts, collapse = "|")
  }
  benchmark_col <- 5L # Export=1, Oracle=2, Metamorphic=3, Adversarial=4, Regression=5? no:
  # Columns: | Export | Oracle | Metamorphic | Adversarial | Regression | Benchmark | Tiers |
  benchmark_col <- 6L
  for (op in unique(names(harness_map))) {
    token <- harness_map[[op]]
    hit <- row_idx[which(row_tokens == token)]
    # Re-locate after any --fix mutation.
    if (length(hit) == 0L || any(hit > length(lines)) ||
        !identical(first_token(lines[hit[1L]]), token)) {
      all_rows <- which(grepl("^\\| `", lines))
      hit <- all_rows[vapply(lines[all_rows], function(l) identical(first_token(l), token), logical(1))]
    }
    if (length(hit) == 0L) {
      problems <- c(problems, sprintf(
        "Harness op '%s' maps to table token '%s' which has no row — fix harness_map or the table",
        op, token
      ))
      next
    }
    line <- lines[hit[1L]]
    cells <- strsplit(line, "|", fixed = TRUE)[[1]]
    bench <- trimws(cells[benchmark_col + 1L])
    if (identical(bench, "—")) {
      if (fix) {
        lines[hit[1L]] <- set_cell(line, benchmark_col, "✓")
        message(sprintf("FIXED: Benchmark cell for '%s' (harness op '%s') — -> ✓", token, op))
      } else {
        problems <- c(problems, sprintf(
          "Harness measures op '%s' but Benchmark cell for '%s' is '—' (run with --fix)",
          op, token
        ))
      }
    } else if (identical(bench, "○")) {
      message(sprintf("NOTE: Benchmark cell for '%s' (harness op '%s') is '○' — review manually", token, op))
    }
  }

  # 3. Orphan rows (single plain-token rows only; alias/method rows skipped).
  plain <- row_tokens[grepl("^[A-Za-z0-9_.]+$", row_tokens) &
                        !grepl(" method|/", lines[row_idx])]
  orphans <- setdiff(plain, exports)
  if (length(orphans) > 0L) {
    message(sprintf("NOTE: table rows with no matching export (review, not auto-removed): %s",
                    paste(unique(orphans), collapse = ", ")))
  }

  if (fix) {
    writeLines(lines, doc_path, useBytes = TRUE)
  }

  if (length(problems) > 0L) {
    cat("Coverage-table sync FAILED:\n", paste0("  - ", problems, collapse = "\n"), "\n", sep = "")
    quit(save = "no", status = 1L)
  }
  cat("Coverage-table sync OK.\n")
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  sync_coverage_main()
}
