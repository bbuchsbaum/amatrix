# amatrix_explain() - human-readable dispatch diagnostic.
#
# Shows exactly which backend was chosen for a given operation, why each
# candidate was accepted or rejected, and actionable suggestions for getting
# the best performance.

amatrix_explain <- function(x, op, y = NULL) {
  stopifnot(inherits(x, "aMatrix"), is.character(op), length(op) == 1L)

  plan  <- amatrix_backend_plan(x, op = op, y = y)
  cal   <- amatrix_calibration_info()
  resid <- amatrix_residency_info(x)

  .amatrix_explain_print(x, op, plan, cal, resid)
  invisible(plan)
}

# -- Formatting helpers --------------------------------------------------------

.amatrix_explain_print <- function(x, op, plan, cal, resid) {
  width <- 72L
  rule  <- function(title = "") {
    pad <- width - nchar(title) - 4L
    paste0("\u2500\u2500 ", title, " ", strrep("\u2500", max(pad, 0L)))
  }

  cat(rule(paste("amatrix dispatch:", op)), "\n")

  # -- Object summary ----------------------------------------------------------
  dims <- paste0("[", nrow(x), " \u00d7 ", ncol(x), "]")
  cat(sprintf("  object:    %s %s  precision=%s  preferred=%s\n",
    class(x)[[1L]], dims, x@precision, x@preferred_backend))

  live_be <- .amatrix_live_resident_backend(x)
  if (!is.null(live_be)) {
    cat(sprintf("  residency: GPU-resident on '%s'\n", live_be))
  } else {
    cat("  residency: host (not GPU-resident)\n")
  }

  cat("\n")

  # -- Per-candidate table -----------------------------------------------------
  cat(rule("candidates"), "\n")
  for (cand in plan$candidates) {
    .amatrix_explain_candidate(cand)
  }
  cat("\n")

  # -- Chosen summary ----------------------------------------------------------
  cat(rule("result"), "\n")
  chosen_path_label <- switch(plan$chosen_path,
    cold     = "cold path (upload + compute)",
    resident = "resident path (already on device)",
    plan$chosen_path
  )
  cat(sprintf("  chosen: %s  via %s\n\n", plan$chosen, chosen_path_label))

  # -- Suggestions -------------------------------------------------------------
  suggestions <- .amatrix_explain_suggestions(x, op, plan, cal, resid)
  if (length(suggestions) > 0L) {
    cat(rule("suggestions"), "\n")
    for (s in suggestions) cat(s, "\n")
    cat("\n")
  }

  cat(strrep("\u2500", width), "\n")
  invisible(NULL)
}

.amatrix_explain_candidate <- function(cand) {
  chosen_marker <- if (cand$chosen) "\u25ba CHOSEN " else "  ......  "

  # Build compact flag string
  flag <- function(val, label) {
    if (isTRUE(val)) label else paste0("NO-", label)
  }
  flags <- paste(
    flag(cand$registered,          "reg"),
    flag(cand$available,           "avail"),
    flag(cand$precision_compatible,"prec"),
    flag(cand$supported_cold,      "cold"),
    flag(cand$calibration_ok,      "calib"),
    if (cand$resident_active) "RESIDENT" else NULL
  )

  path_label <- if (is.na(cand$chosen_path)) "" else
    sprintf("  [%s]", cand$chosen_path)

  cat(sprintf("  %s  %-14s %s%s\n",
    chosen_marker, cand$name, flags, path_label))
}

.amatrix_explain_suggestions <- function(x, op, plan, cal, resid) {
  s <- character()
  add <- function(...) { s <<- c(s, paste0(...)) }

  chosen      <- plan$chosen
  chosen_path <- plan$chosen_path
  preferred   <- plan$preferred[[1L]]

  # -- Calibration missing -----------------------------------------------------
  if (is.null(cal) && preferred != "cpu") {
    add("* No calibration data found. GPU dispatch uses static capability")
    add("  checks only -- GPU may be chosen even when CPU is faster.")
    add("  \u2192 Run amatrix_calibrate() once to tune thresholds for this machine.")
  }

  # -- CPU fallback analysis ---------------------------------------------------
  if (chosen == "cpu" && preferred != "cpu") {
    known_backends <- c("cpu", names(.amatrix_optional_backend_specs()))
    gpu_cands <- Filter(
      function(c) c$name != "cpu" && c$name %in% known_backends,
      plan$candidates
    )

    for (cand in gpu_cands) {
      be <- cand$name
      if (!cand$registered) {
        add(sprintf("* '%s' not registered (package not installed?).", be))
        add(sprintf("  \u2192 Install amatrix.%s and retry.", be))
        next
      }
      if (!cand$available) {
        add(sprintf("* '%s' registered but reports available=FALSE.", be))
        add(sprintf("  \u2192 Check backend status: amatrix_backend_status('%s')", be))
        next
      }
      if (!cand$precision_compatible) {
        add(sprintf(
          "* '%s' requires precision 'fast' but object uses '%s'.",
          be, x@precision))
        add("  \u2192 Reconstruct with mode=\"fast\":")
        add(sprintf("      adgeMatrix(as.matrix(x), mode=\"fast\", backend=\"%s\")", be))
        next
      }
      if (!isTRUE(cand$calibration_ok)) {
        n_elem   <- nrow(x) * ncol(x)
        thresh   <- cal$thresholds[[be]][[op]]
        thresh_s <- if (!is.null(thresh) && !is.infinite(thresh))
          formatC(thresh, format = "d", big.mark = ",") else "unknown"
        add(sprintf(
          "* '%s' skipped by calibration: %s elements < threshold %s.",
          be, formatC(n_elem, format="d", big.mark=","), thresh_s))
        add("  \u2192 Use a larger matrix, or re-run amatrix_calibrate() if the")
        add("    threshold seems wrong for your hardware.")
        next
      }
      if (!cand$supported_cold) {
        add(sprintf("* '%s' does not support op '%s'.", be, op))
        next
      }
    }
  }

  # -- Cold path on GPU --------------------------------------------------------
  if (chosen != "cpu" && identical(chosen_path, "cold")) {
    add(sprintf(
      "* GPU cold path: matrix will be uploaded to '%s' each call.", chosen))
    add("  \u2192 For repeated ops on the same matrix, bind it resident first:")
    add(sprintf(
      "      X <- amatrix_bind_resident(X, \"%s\")", chosen))
    add("  \u2192 For the many-Y QR workflow, use:")
    add("      many_lm(X, Y, method=\"qr\", cache=TRUE)")
  }

  # -- Resident path: already optimal ------------------------------------------
  if (identical(chosen_path, "resident")) {
    add(sprintf(
      "* Optimal: '%s' resident path active. No upload cost.", chosen))
  }

  # -- Warm-up reminder --------------------------------------------------------
  if (chosen != "cpu" && is.null(.amatrix_state$warm_done[[chosen]])) {
    add(sprintf(
      "* First call to '%s' may be slow (JIT kernel compilation).", chosen))
    add(sprintf("  \u2192 Run amatrix_warm(\"%s\") before timed work.", chosen))
  }

  s
}
