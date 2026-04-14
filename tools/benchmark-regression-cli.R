#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x
r_string_literal <- function(x) encodeString(x, quote = "\"")
r_character_vector_literal <- function(x) {
  if (length(x) == 0L) {
    return("character()")
  }

  sprintf("c(%s)", paste(vapply(x, r_string_literal, character(1)), collapse = ", "))
}

raw_args <- commandArgs(trailingOnly = FALSE)
direct_file_paths <- sub("^--file=", "", grep("^--file=", raw_args, value = TRUE))
if (length(direct_file_paths) == 0L) {
  stop("benchmark-regression-cli.R must be invoked as a file", call. = FALSE)
}

script_path <- normalizePath(direct_file_paths[[1L]], winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(dirname(dirname(script_path)), winslash = "/", mustWork = TRUE)

setwd(repo_root)
options(amatrix.benchmark_regression.script_path = normalizePath(file.path("tools", "benchmark-regression.R"), winslash = "/", mustWork = TRUE))

script_target <- normalizePath(file.path("tools", "benchmark-regression.R"), winslash = "/", mustWork = TRUE)
dispatch_args <- commandArgs(trailingOnly = TRUE)
expr <- sprintf(
  paste(
    "setwd(%s); source(%s, local = globalenv());",
    "args <- parse_args(%s);",
    "initialize_regression_benchmark_context();",
    "if (isTRUE(args$worker)) run_worker(args) else run_master(args)"
  ),
  r_string_literal(repo_root),
  r_string_literal(script_target),
  r_character_vector_literal(dispatch_args)
)

warned_status <- NULL
launch_output <- withCallingHandlers(
  system2(
    file.path(R.home("bin"), "Rscript"),
    vapply(c("-e", expr), shQuote, character(1), USE.NAMES = FALSE),
    stdout = TRUE,
    stderr = TRUE
  ),
  warning = function(w) {
    warned_status <<- attr(w, "status") %||% warned_status
    invokeRestart("muffleWarning")
  }
)
launch_status <- attr(launch_output, "status") %||% warned_status %||% 0L

if (length(launch_output) > 0L) {
  cat(paste(launch_output, collapse = "\n"), sep = "\n")
  if (!grepl("\n$", paste(launch_output, collapse = "\n"))) {
    cat("\n")
  }
}

quit(save = "no", status = launch_status)
