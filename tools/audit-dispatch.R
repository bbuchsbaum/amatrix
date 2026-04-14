#!/usr/bin/env Rscript
#
# tools/audit-dispatch.R — Track 3 dispatch hermeticity audit
#
# Enforces that every dispatch-sensitive generic has an explicit setMethod
# registration for every mixed (plain LHS, amatrix RHS) and (amatrix LHS,
# plain RHS) signature. Without these, S4 dispatch falls through to base R
# or Matrix generics which coerce the amatrix operand to a plain matrix and
# silently destroy GPU residency.
#
# Target generics: `%*%`, `crossprod`, `tcrossprod`.
#
# Target classes:
#   plain:   matrix, numeric, dgeMatrix, dgCMatrix
#   amatrix: adgeMatrix, adgCMatrix
#
# Required grid: every mixed pair (exactly one plain side, one amatrix side).
#
# The audit walks R/*.R with the R parser (NOT regex) to extract setMethod
# registrations, then cross-checks against the required grid. "ANY" and
# ancestor-class signatures count: if a method is registered for signature
# (aMatrix, ANY) and the target is (adgeMatrix, matrix), the method covers it.
#
# Exit codes:
#   0 — all required signatures present
#   1 — one or more gaps detected (PR-blocking)
#   2 — parser error or missing R/ directory
#
# Invocation:
#   Rscript tools/audit-dispatch.R
#   Rscript tools/audit-dispatch.R --json    # machine-readable output
#   Rscript tools/audit-dispatch.R --verbose # list every required pair

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
output_json <- "--json" %in% args
verbose <- "--verbose" %in% args

# Resolve repo root from the --file= argument so the audit can be invoked
# from anywhere (including CI without a cd step).
raw_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", raw_args, value = TRUE)
script_path <- if (length(file_arg) >= 1L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path("tools", "audit-dispatch.R"), mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
r_dir <- file.path(repo_root, "R")

if (!dir.exists(r_dir)) {
  message("FATAL: R/ directory not found at ", r_dir)
  quit(status = 2L)
}

# ---------------------------------------------------------------------------
# Parse R/*.R and collect every setMethod(generic, signature, body) call.
# ---------------------------------------------------------------------------

parse_signature <- function(sig_expr) {
  # Accepts:
  #   "classname"                      — single positional
  #   c("c1", "c2")                    — positional multi
  #   signature(x = "c1", y = "c2")    — named
  if (is.character(sig_expr) && length(sig_expr) == 1L) {
    return(sig_expr)
  }
  if (!is.call(sig_expr)) return(NULL)

  head_name <- as.character(sig_expr[[1L]])
  if (!(head_name %in% c("signature", "c"))) return(NULL)

  parts <- as.list(sig_expr)[-1L]
  out <- vapply(parts, function(p) {
    if (is.character(p) && length(p) == 1L) p
    else if (is.name(p)) as.character(p)
    else NA_character_
  }, character(1))

  if (any(is.na(out))) return(NULL)

  # Reorder named signatures to canonical (x, y, e1, e2) order so we can
  # compare positionally with target pairs.
  nm <- names(out)
  if (!is.null(nm) && any(nzchar(nm))) {
    canonical <- intersect(c("x", "y", "e1", "e2"), nm)
    rest <- setdiff(nm, canonical)
    out <- out[c(canonical, rest)]
  }
  unname(out)
}

extract_setMethods <- function(file) {
  exprs <- tryCatch(
    parse(file = file, keep.source = FALSE),
    error = function(e) {
      message(sprintf("WARN: could not parse %s: %s", file, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(exprs)) return(list())

  results <- list()
  for (e in as.list(exprs)) {
    if (!is.call(e)) next
    head <- e[[1L]]
    if (!identical(head, as.name("setMethod"))) next
    if (length(e) < 3L) next

    generic <- tryCatch(
      eval(e[[2L]], envir = baseenv()),
      error = function(err) NULL
    )
    if (is.null(generic) || !is.character(generic) || length(generic) != 1L) next

    sig <- parse_signature(e[[3L]])
    if (is.null(sig)) next

    results[[length(results) + 1L]] <- list(
      file = basename(file),
      generic = generic,
      sig = sig
    )
  }
  results
}

r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE, recursive = FALSE)
all_methods <- unlist(lapply(r_files, extract_setMethods), recursive = FALSE)
method_index <- split(all_methods, vapply(all_methods, `[[`, character(1), "generic"))

# ---------------------------------------------------------------------------
# Class hierarchy for ancestor-match coverage.
# ---------------------------------------------------------------------------

class_ancestors <- list(
  adgeMatrix     = c("adgeMatrix", "aMatrix", "dgeMatrix", "Matrix", "ANY"),
  adgCMatrix     = c("adgCMatrix", "aMatrix", "dgCMatrix", "Matrix", "ANY"),
  aTransposeView = c("aTransposeView", "aMatrix", "ANY"),
  aMatrix        = c("aMatrix", "ANY"),
  matrix         = c("matrix", "ANY"),
  numeric        = c("numeric", "vector", "ANY"),
  dgeMatrix      = c("dgeMatrix", "Matrix", "ANY"),
  dgCMatrix      = c("dgCMatrix", "Matrix", "ANY"),
  Matrix         = c("Matrix", "ANY"),
  ANY            = c("ANY")
)

class_covers_target <- function(registered, target) {
  anc <- class_ancestors[[target]] %||% c(target, "ANY")
  registered %in% anc
}

method_covers_signature <- function(method_sig, target_sig) {
  if (length(method_sig) != length(target_sig)) return(FALSE)
  all(vapply(
    seq_along(target_sig),
    function(i) class_covers_target(method_sig[i], target_sig[i]),
    logical(1)
  ))
}

is_signature_covered <- function(generic, target_sig) {
  methods <- method_index[[generic]]
  if (is.null(methods)) return(FALSE)
  for (m in methods) {
    if (method_covers_signature(m$sig, target_sig)) return(TRUE)
  }
  FALSE
}

# ---------------------------------------------------------------------------
# Build the required grid and audit.
# ---------------------------------------------------------------------------

target_generics <- c("%*%", "crossprod", "tcrossprod")
plain_classes   <- c("matrix", "numeric", "dgeMatrix", "dgCMatrix")
amatrix_classes <- c("adgeMatrix", "adgCMatrix")

target_pairs <- list()
for (lhs in c(plain_classes, amatrix_classes)) {
  for (rhs in c(plain_classes, amatrix_classes)) {
    is_lhs_plain <- lhs %in% plain_classes
    is_rhs_plain <- rhs %in% plain_classes
    if (xor(is_lhs_plain, is_rhs_plain)) {  # exactly one plain side
      target_pairs[[length(target_pairs) + 1L]] <- c(lhs, rhs)
    }
  }
}

findings <- list()
for (generic in target_generics) {
  for (pair in target_pairs) {
    covered <- is_signature_covered(generic, pair)
    findings[[length(findings) + 1L]] <- list(
      generic = generic,
      lhs = pair[1L],
      rhs = pair[2L],
      covered = covered
    )
  }
}

missing_entries <- Filter(function(f) !f$covered, findings)

# ---------------------------------------------------------------------------
# Report.
# ---------------------------------------------------------------------------

if (output_json) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    message("--json requested but jsonlite not available; falling back to text")
    output_json <- FALSE
  }
}

if (output_json) {
  payload <- list(
    repo_root = repo_root,
    target_generics = target_generics,
    plain_classes = plain_classes,
    amatrix_classes = amatrix_classes,
    total_required = length(findings),
    total_covered = sum(vapply(findings, `[[`, logical(1), "covered")),
    setMethod_calls_indexed = length(all_methods),
    missing = lapply(missing_entries, function(m) {
      list(generic = m$generic, lhs = m$lhs, rhs = m$rhs)
    })
  )
  cat(jsonlite::toJSON(payload, pretty = TRUE, auto_unbox = TRUE), "\n", sep = "")
} else {
  cat("tools/audit-dispatch.R — dispatch hermeticity audit\n")
  cat(sprintf("Repo:                       %s\n", repo_root))
  cat(sprintf("Generics audited:           %s\n", paste(target_generics, collapse = ", ")))
  cat(sprintf("setMethod calls indexed:    %d\n", length(all_methods)))
  cat(sprintf("Required (mixed) pairs:     %d\n", length(findings)))
  cat(sprintf("Covered pairs:              %d\n",
              sum(vapply(findings, `[[`, logical(1), "covered"))))
  cat(sprintf("Missing pairs:              %d\n", length(missing_entries)))
  cat("\n")

  if (verbose) {
    cat("All required pairs (including covered):\n")
    for (f in findings) {
      mark <- if (f$covered) "  [ok]" else "  [MISSING]"
      cat(sprintf("%s %s(%s, %s)\n", mark, f$generic, f$lhs, f$rhs))
    }
    cat("\n")
  }

  if (length(missing_entries) > 0L) {
    cat("MISSING dispatch methods (PR-blocking):\n\n")
    by_gen <- split(
      missing_entries,
      vapply(missing_entries, `[[`, character(1), "generic")
    )
    for (g in names(by_gen)) {
      cat(sprintf("  %s:\n", g))
      for (entry in by_gen[[g]]) {
        cat(sprintf("    setMethod(\"%s\", signature(x = \"%s\", y = \"%s\"), ...)\n",
                    g, entry$lhs, entry$rhs))
      }
      cat("\n")
    }
    cat("Every missing pair allows S4 dispatch to fall through to base or Matrix\n")
    cat("generics, coercing the amatrix operand and silently losing GPU residency.\n")
    cat("Add the missing methods in R/dispatch-hardening.R. See also\n")
    cat("planning_docs/quality-tracking.md \u00a77 rule 6.\n")
  } else {
    cat("All required dispatch pairs are covered. Dispatch hermeticity OK.\n")
  }
}

if (length(missing_entries) > 0L) {
  quit(status = 1L)
}
quit(status = 0L)
