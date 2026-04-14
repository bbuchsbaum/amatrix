#!/usr/bin/env Rscript
# tools/lint-anti-patterns.R
#
# Maintainable lint script derived from round-2 bug-hunt anti-pattern greps.
# Run in CI or interactively. Exits non-zero if HIGH-severity new hits found.
#
# Usage:
#   Rscript tools/lint-anti-patterns.R                  # report mode
#   Rscript tools/lint-anti-patterns.R --baseline       # seed / refresh baseline
#   Rscript tools/lint-anti-patterns.R --json /tmp/out.json  # write JSON to path
#   Rscript tools/lint-anti-patterns.R --no-exit        # never exit non-zero (useful locally)
#
# CI integration:
#   Add to .github/workflows as a step:
#     - name: Anti-pattern lint
#       run: Rscript tools/lint-anti-patterns.R
#
# Dependencies: base R only (jsonlite used only if installed, for prettier JSON)

suppressPackageStartupMessages({
  has_jsonlite <- requireNamespace("jsonlite", quietly = TRUE)
})

# ---------------------------------------------------------------------------
# CLI parsing
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
opt_baseline   <- "--baseline"   %in% args
opt_no_exit    <- "--no-exit"    %in% args
opt_json_idx   <- which(args == "--json")
opt_json_path  <- if (length(opt_json_idx) > 0 && length(args) >= opt_json_idx + 1)
                    args[[opt_json_idx + 1]] else NULL

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
script_dir     <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) "tools")
root           <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
r_glob         <- file.path(root, "R")
baseline_path  <- file.path(root, "tools", "lint-anti-patterns-baseline.json")

# Collect all R source files
all_r_files <- list.files(r_glob, pattern = "\\.R$", full.names = TRUE)

# ---------------------------------------------------------------------------
# Pattern definitions
# Each entry:
#   name        - short id used in baseline keys
#   description - one-line human description
#   regex       - PCRE pattern passed to grep(perl=TRUE)
#   glob        - file filter (character vector of full paths)
#   excludes    - file basenames that are intentionally exempt
#   severity    - "HIGH" | "MED" | "LOW"
#   explanation - why this is a smell
#   fix         - suggested fix
# ---------------------------------------------------------------------------
patterns <- list(

  list(
    name        = "P01_unclassed_stop",
    description = "stop() with bare string message — no condition class",
    regex       = 'stop\\("',
    glob        = all_r_files,
    excludes    = character(0),
    severity    = "HIGH",
    explanation = paste(
      "Callers cannot catch GPU-specific failures by class.",
      "tryCatch(expr, amatrix_error = function(e) ...) will never fire.",
      "All 90+ sites use fragile message-string matching instead."
    ),
    fix         = "Use stop(structure(list(message=msg, call=sys.call()), class=c('amatrix_error','error','condition')))"
  ),

  list(
    name        = "P02_trycatch_swallow_null",
    description = "tryCatch error handler silently returns NULL",
    regex       = 'error\\s*=\\s*function\\s*\\([^)]*\\)\\s*NULL',
    glob        = all_r_files,
    # backend-warmup.R and backend-calibration.R and memory-stats.R are
    # intentional best-effort probes. backend-calibration.R warmup loop
    # also intentional. All other swallows should be reviewed.
    excludes    = c("backend-warmup.R", "backend-calibration.R", "memory-stats.R"),
    severity    = "HIGH",
    explanation = paste(
      "Silent NULL return hides errors from callers.",
      "Intentional sites: backend probes in .amatrix_get_backend().",
      "Bug sites: wrappers.R double-drop, chol-factor.R out_key leak, qr.R NULL Q/R."
    ),
    fix         = "Log with message() or signal an amatrix_fallback condition before returning NULL"
  ),

  list(
    name        = "P03_host_cache_valid_set_true",
    description = "host_cache_valid <- TRUE (never a matching FALSE site)",
    regex       = 'host_cache_valid\\s*<-\\s*TRUE',
    glob        = all_r_files,
    excludes    = character(0),
    severity    = "HIGH",
    explanation = paste(
      "host_cache_valid is SET TRUE at bind time but NEVER SET FALSE.",
      "In-place resident mutations leave cache flag stale;",
      "next amatrix_materialize_dense() returns stale host data."
    ),
    fix         = "Set host_cache_valid <- FALSE in every path that mutates the resident buffer in place"
  ),

  list(
    name        = "P04_host_cache_valid_never_false",
    description = "host_cache_valid <- FALSE absent (sentinel check — expect 0 hits is BAD)",
    regex       = 'host_cache_valid\\s*<-\\s*FALSE',
    glob        = all_r_files,
    excludes    = character(0),
    severity    = "HIGH",
    explanation = paste(
      "Complement to P03: if this count is 0 the cache invalidation path is missing.",
      "A count of 0 here IS the bug."
    ),
    fix         = "Add host_cache_valid <- FALSE in broadcast_ewise_resident_inplace_key and any other in-place mutator"
  ),

  list(
    name        = "P05_double_drop_pattern",
    description = "resident_drop called inside tryCatch AND again unconditionally after",
    # Matches lines where resident_drop is called with out_key or scaled_key
    regex       = 'resident_drop\\([^)]*(?:out_key|scaled_key|temp_key)[^)]*\\)',
    glob        = all_r_files,
    excludes    = character(0),
    severity    = "HIGH",
    explanation = paste(
      "Key dropped inside tryCatch error handler, then unconditionally dropped again after.",
      "Double-drop corrupts the backend key registry.",
      "Confirmed in chol-factor.R:175-191, wrappers.R:1283-1299, wrappers.R:1345-1359."
    ),
    fix         = "Use on.exit(backend$resident_drop(key), add=TRUE) and remove manual drops"
  ),

  list(
    name        = "P06_resident_key_alloc_no_on_exit",
    description = "resident key allocated without on.exit guard",
    regex       = '\\.amatrix_next_resident_key\\(|resident_next_key\\(',
    glob        = all_r_files,
    # irlba.R explicitly uses on.exit for resident handles; wrappers.R is the problem
    excludes    = c("irlba.R"),
    severity    = "MED",
    explanation = paste(
      "~40 sites in wrappers.R pre-allocate out_key/scaled_key inside tryCatch blocks",
      "but rely on manual try-drop after the block.",
      "If a second error fires between tryCatch return and manual drop, key leaks permanently."
    ),
    fix         = "Register on.exit(try(backend$resident_drop(key)), add=TRUE) immediately after allocation"
  ),

  list(
    name        = "P07_as_matrix_dimnames_drop",
    description = "as.matrix() call that may silently drop dimnames",
    regex       = 'as\\.matrix\\(',
    glob        = all_r_files,
    # backend-cpu.R discards dimnames intentionally (re-wrapped by caller)
    excludes    = c("backend-cpu.R"),
    severity    = "LOW",
    explanation = paste(
      "as.matrix() on adgeMatrix/aTransposeView returns plain matrix without dimnames",
      "unless the S4 method is dispatched.",
      "High-severity sites: residency.R:396,403 (deferred path), qr.R:726 (RHS backsolve)."
    ),
    fix         = "Use amatrix_materialize_host() or explicitly copy x@Dimnames after as.matrix()"
  ),

  list(
    name        = "P08_missing_call_false_stop",
    description = "stop() without call.=FALSE leaking internal call stack",
    regex       = 'stop\\([^)]+\\)(?!.*call\\.\\s*=\\s*FALSE)',
    glob        = all_r_files,
    excludes    = character(0),
    severity    = "LOW",
    explanation = paste(
      "stop() without call.=FALSE attaches the internal R call to the condition.",
      "Users see amatrix internal function names in their error messages."
    ),
    fix         = "Add call.=FALSE to all user-facing stop() calls"
  ),

  list(
    name        = "P09_s4_nextmethod_fallthrough",
    description = "NextMethod() in S4 method — likely S3/S4 dispatch confusion",
    regex       = '(?<!call)NextMethod\\(\\)',
    glob        = all_r_files,
    excludes    = character(0),
    severity    = "MED",
    explanation = paste(
      "NextMethod() is an S3 idiom. In S4 methods, callNextMethod() is correct.",
      "Using NextMethod() in an S4 context silently falls through to the wrong method."
    ),
    fix         = "Replace NextMethod() with callNextMethod() inside setMethod() blocks"
  ),

  list(
    name        = "P10_trycatch_c_bridge_swallow",
    description = ".Call() inside tryCatch(error=function(...) NULL) — C errors silently swallowed",
    regex       = '\\.Call\\(',
    glob        = all_r_files,
    excludes    = character(0),
    severity    = "MED",
    explanation = paste(
      ".Call() to C bridges can partially mutate internal state before throwing.",
      "Swallowing the error leaves R-side state inconsistent.",
      "Confirmed: irlba.R Lanczos bridge, amatrix_compile_product in svd-factor.R."
    ),
    fix         = "Let .Call() errors propagate or wrap in withCallingHandlers() that logs before re-signalling"
  )

)

# ---------------------------------------------------------------------------
# Grep engine
# ---------------------------------------------------------------------------
grep_pattern <- function(pat, files) {
  hits <- list()
  for (f in files) {
    bname <- basename(f)
    if (bname %in% pat$excludes) next
    lines <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
    matched <- grep(pat$regex, lines, perl = TRUE)
    for (ln in matched) {
      hits[[length(hits) + 1L]] <- list(
        file    = f,
        relpath = sub(paste0(root, "/"), "", f, fixed = TRUE),
        line    = ln,
        context = trimws(lines[[ln]])
      )
    }
  }
  hits
}

# ---------------------------------------------------------------------------
# Run all patterns
# ---------------------------------------------------------------------------
cat("# amatrix anti-pattern lint\n\n")
cat(sprintf("Scanning %d R source files in %s\n\n", length(all_r_files), r_glob))

results <- list()
t0 <- proc.time()

for (pat in patterns) {
  hits <- grep_pattern(pat, pat$glob)
  results[[pat$name]] <- list(
    pattern     = pat,
    hits        = hits,
    hit_count   = length(hits)
  )
}

elapsed <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("Scan completed in %.2f seconds\n\n", elapsed))

# ---------------------------------------------------------------------------
# Load baseline if it exists
# ---------------------------------------------------------------------------
baseline <- list()
if (file.exists(baseline_path)) {
  raw <- tryCatch(
    if (has_jsonlite) jsonlite::fromJSON(baseline_path, simplifyVector = FALSE)
    else              {
      # minimal JSON reader for simple numeric values
      txt <- paste(readLines(baseline_path, warn = FALSE), collapse = "\n")
      # Parse as named list of integers via crude extraction
      keys <- regmatches(txt, gregexpr('"[A-Za-z0-9_]+"\\s*:\\s*[0-9]+', txt))[[1]]
      parsed <- list()
      for (kv in keys) {
        parts <- strsplit(kv, ":\\s*", perl = TRUE)[[1]]
        k <- gsub('"', '', parts[[1]])
        v <- as.integer(trimws(parts[[2]]))
        parsed[[k]] <- v
      }
      parsed
    },
    error = function(e) list()
  )
  baseline <- raw
}

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
cat("---\n\n")
cat("## Results by Pattern\n\n")

total_hits    <- 0L
high_new_hits <- 0L
summary_rows  <- list()

for (nm in names(results)) {
  r   <- results[[nm]]
  pat <- r$pattern
  cnt <- r$hit_count
  baseline_cnt <- if (!is.null(baseline[[nm]])) as.integer(baseline[[nm]]) else NA_integer_

  net_new <- if (!is.na(baseline_cnt)) cnt - baseline_cnt else NA_integer_
  is_new  <- !is.na(net_new) && net_new > 0

  total_hits <- total_hits + cnt
  if (is_new && pat$severity == "HIGH") high_new_hits <- high_new_hits + net_new

  # Special logic for P04: 0 hits IS the problem
  p04_flag <- nm == "P04_host_cache_valid_never_false" && cnt == 0

  cat(sprintf("### %s [%s]\n", pat$name, pat$severity))
  cat(sprintf("**%s**\n\n", pat$description))
  cat(sprintf("- Hits: **%d**", cnt))
  if (!is.na(baseline_cnt)) {
    delta_str <- if (net_new >= 0) paste0("+", net_new) else as.character(net_new)
    cat(sprintf(" | Baseline: %d | Net-new: %s", baseline_cnt, delta_str))
  }
  cat("\n")
  cat(sprintf("- Severity: %s | Excludes: %s\n",
              pat$severity,
              if (length(pat$excludes) == 0) "none" else paste(pat$excludes, collapse = ", ")))
  cat(sprintf("- Explanation: %s\n", pat$explanation))
  cat(sprintf("- Fix: `%s`\n\n", pat$fix))

  if (p04_flag) {
    cat("**ALARM: 0 hits for P04 means host_cache_valid is NEVER set FALSE — confirmed bug.**\n\n")
  }

  if (cnt > 0 && cnt <= 20) {
    cat("| File | Line | Context |\n|------|------|---------|\n")
    for (h in r$hits) {
      ctx <- substr(h$context, 1, 80)
      ctx <- gsub("|", "\\|", ctx, fixed = TRUE)
      cat(sprintf("| `%s` | %d | `%s` |\n", h$relpath, h$line, ctx))
    }
    cat("\n")
  } else if (cnt > 20) {
    cat(sprintf("_(showing first 10 of %d hits)_\n\n", cnt))
    cat("| File | Line | Context |\n|------|------|---------|\n")
    for (h in r$hits[seq_len(min(10L, cnt))]) {
      ctx <- substr(h$context, 1, 80)
      ctx <- gsub("|", "\\|", ctx, fixed = TRUE)
      cat(sprintf("| `%s` | %d | `%s` |\n", h$relpath, h$line, ctx))
    }
    cat("\n")
  }

  summary_rows[[length(summary_rows) + 1L]] <- list(
    pattern     = nm,
    severity    = pat$severity,
    hits        = cnt,
    baseline    = baseline_cnt,
    net_new     = net_new,
    p04_alarm   = p04_flag
  )
}

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
cat("---\n\n## Summary\n\n")
cat("| Pattern | Severity | Hits | Baseline | Net-new |\n")
cat("|---------|----------|------|----------|---------|\n")
for (row in summary_rows) {
  base_str <- if (is.na(row$baseline)) "—" else as.character(row$baseline)
  new_str  <- if (is.na(row$net_new))  "—" else {
    if (row$net_new > 0) paste0("**+", row$net_new, "**")
    else as.character(row$net_new)
  }
  cat(sprintf("| %s | %s | %d | %s | %s |\n",
              row$pattern, row$severity, row$hits, base_str, new_str))
}
cat(sprintf("\n**Total hits: %d | High-severity net-new: %d**\n\n",
            total_hits, high_new_hits))

# ---------------------------------------------------------------------------
# Baseline mode: write current counts as new baseline
# ---------------------------------------------------------------------------
if (opt_baseline) {
  baseline_data <- setNames(
    lapply(results, function(r) r$hit_count),
    names(results)
  )
  json_str <- if (has_jsonlite) {
    jsonlite::toJSON(baseline_data, auto_unbox = TRUE, pretty = TRUE)
  } else {
    # hand-roll minimal JSON
    pairs <- mapply(function(k, v) sprintf('  "%s": %d', k, v),
                    names(baseline_data), unlist(baseline_data))
    paste0("{\n", paste(pairs, collapse = ",\n"), "\n}")
  }
  writeLines(json_str, baseline_path)
  cat(sprintf("Baseline written to %s\n", baseline_path))
}

# ---------------------------------------------------------------------------
# Optional JSON output for machine consumers
# ---------------------------------------------------------------------------
json_out_path <- opt_json_path %||% file.path(tempdir(), "lint-anti-patterns.json")
`%||%` <- function(a, b) if (!is.null(a)) a else b
json_out_path <- if (!is.null(opt_json_path)) opt_json_path else file.path(tempdir(), "lint-anti-patterns.json")

machine_out <- list(
  generated   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  total_hits  = total_hits,
  high_new    = high_new_hits,
  patterns    = lapply(results, function(r) list(
    name      = r$pattern$name,
    severity  = r$pattern$severity,
    hit_count = r$hit_count,
    hits      = lapply(r$hits, function(h) list(
      file    = h$relpath,
      line    = h$line,
      context = h$context
    ))
  ))
)

json_str_out <- if (has_jsonlite) {
  jsonlite::toJSON(machine_out, auto_unbox = TRUE, pretty = TRUE)
} else {
  paste0('{"generated":"', machine_out$generated,
         '","total_hits":', machine_out$total_hits,
         ',"high_new":', machine_out$high_new, '}')
}
writeLines(json_str_out, json_out_path)
cat(sprintf("\nMachine-readable JSON written to: %s\n", json_out_path))

# ---------------------------------------------------------------------------
# Exit code
# ---------------------------------------------------------------------------
if (!opt_no_exit && !opt_baseline && high_new_hits > 0) {
  cat(sprintf("\nFAIL: %d HIGH-severity net-new hit(s) found. Fix or update baseline.\n",
              high_new_hits))
  quit(status = 1L, save = "no")
} else if (!opt_baseline) {
  cat("\nPASS: No HIGH-severity net-new hits (or --no-exit set).\n")
}
