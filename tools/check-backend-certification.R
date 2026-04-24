#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

r_string_literal <- function(x) encodeString(x, quote = "\"")

r_character_vector_literal <- function(x) {
  if (length(x) == 0L) {
    return("character()")
  }
  sprintf("c(%s)", paste(vapply(x, r_string_literal, character(1)), collapse = ", "))
}

backend_certification_bootstrap_direct_entry <- function() {
  raw_args <- commandArgs(trailingOnly = FALSE)
  direct_file_paths <- sub("^--file=", "", grep("^--file=", raw_args, value = TRUE))
  direct_file_paths <- direct_file_paths[basename(direct_file_paths) == "check-backend-certification.R"]
  if (length(direct_file_paths) == 0L ||
      "--safe-main" %in% commandArgs(trailingOnly = TRUE)) {
    return(invisible(FALSE))
  }

  script_path <- normalizePath(direct_file_paths[[1L]], winslash = "/", mustWork = TRUE)
  repo_root <- normalizePath(dirname(dirname(script_path)), winslash = "/", mustWork = TRUE)
  dispatch_args <- c("--safe-main", commandArgs(trailingOnly = TRUE))
  expr <- sprintf(
    paste(
      "setwd(%s);",
      "Sys.setenv(AMATRIX_BACKEND_CERTIFICATION_AUTORUN = \"0\");",
      "source(%s, local = globalenv());",
      "backend_certification_main(%s)"
    ),
    r_string_literal(repo_root),
    r_string_literal(script_path),
    r_character_vector_literal(dispatch_args)
  )

  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    vapply(c("-e", expr), shQuote, character(1), USE.NAMES = FALSE)
  )
  quit(save = "no", status = status)
}

backend_certification_bootstrap_direct_entry()

backend_certification_usage <- function() {
  paste(
    "amatrix backend certification gate",
    "",
    "Usage:",
    "  Rscript -e 'source(\"tools/check-backend-certification.R\")'",
    "  Rscript tools/check-backend-certification.R [options]",
    "",
    "Options:",
    "  --help                 Show this help and exit",
    "  --backends=a,b         Backend gates to run (default: all)",
    "                         Valid: mlx,arrayfire,opencl,metal,all",
    "  --allow-skips          Do not fail when a selected gate has skipped tests",
    "  --summary=PATH         Write a CSV summary",
    "",
    "Notes:",
    "  - The script relaunches direct file-entry invocations through Rscript -e",
    "    so MLX startup follows the safe Apple Silicon probe path.",
    "  - Release runs should use the default no-skip policy.",
    sep = "\n"
  )
}

backend_certification_gates <- function() {
  list(
    mlx = list(
      filter = "backend-certification-mlx",
      env = character()
    ),
    arrayfire = list(
      filter = "cross-backend-conformance|arrayfire-matmul-layout|regression-arrayfire-worker-crash",
      env = character()
    ),
    opencl = list(
      filter = "cross-backend-conformance|opencl-model-core|opencl-eigen|benchmark-harness",
      env = c(AMATRIX_OPENCL_PROBE_GPU = "1")
    ),
    metal = list(
      filter = "backend-certification-metal|sparse-backend|sparse-product-pathway|sparse-linalg",
      env = c(AMATRIX_METAL_PROBE_GPU = "1")
    )
  )
}

parse_backend_certification_args <- function(args) {
  out <- list(
    help = any(args %in% c("--help", "-h")),
    safe_main = "--safe-main" %in% args,
    backends = names(backend_certification_gates()),
    allow_skips = "--allow-skips" %in% args,
    summary = NULL
  )

  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) {
      next
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- pieces[[1L]]
    value <- paste(pieces[-1L], collapse = "=")
    if (identical(key, "backends")) {
      requested <- unique(strsplit(value, ",", fixed = TRUE)[[1L]])
      requested <- requested[nzchar(requested)]
      out$backends <- if ("all" %in% requested) names(backend_certification_gates()) else requested
    }
    if (identical(key, "summary")) {
      out$summary <- value
    }
  }

  valid <- names(backend_certification_gates())
  unknown <- setdiff(out$backends, valid)
  if (length(unknown) > 0L) {
    stop(sprintf("unknown backend gate(s): %s", paste(unknown, collapse = ", ")), call. = FALSE)
  }

  out
}

with_backend_certification_env <- function(env, code) {
  if (length(env) == 0L) {
    return(force(code))
  }

  old <- Sys.getenv(names(env), unset = NA_character_)
  on.exit({
    for (idx in seq_along(env)) {
      name <- names(env)[[idx]]
      if (is.na(old[[idx]])) {
        Sys.unsetenv(name)
      } else {
        do.call(Sys.setenv, stats::setNames(as.list(old[[idx]]), name))
      }
    }
  }, add = TRUE)

  for (idx in seq_along(env)) {
    do.call(Sys.setenv, stats::setNames(as.list(env[[idx]]), names(env)[[idx]]))
  }

  force(code)
}

run_backend_certification_gate <- function(name, gate) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools is required for backend certification gates", call. = FALSE)
  }

  start <- Sys.time()
  result <- with_backend_certification_env(gate$env, {
    devtools::test(filter = gate$filter, reporter = "summary", stop_on_failure = FALSE)
  })
  df <- as.data.frame(result)

  data.frame(
    backend = name,
    filter = gate$filter,
    contexts = nrow(df),
    failed = sum(df$failed %||% 0L),
    errors = sum(df$error %||% 0L),
    skipped = sum(df$skipped %||% 0L),
    warnings = sum(df$warning %||% 0L),
    duration_sec = round(as.numeric(difftime(Sys.time(), start, units = "secs")), 3),
    stringsAsFactors = FALSE
  )
}

backend_certification_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  parsed <- parse_backend_certification_args(args)
  if (isTRUE(parsed$help)) {
    cat(backend_certification_usage(), "\n")
    return(invisible(0L))
  }

  gates <- backend_certification_gates()[parsed$backends]
  summaries <- lapply(names(gates), function(name) {
    cat(sprintf("\n== Backend certification: %s ==\n", name))
    run_backend_certification_gate(name, gates[[name]])
  })
  summary <- do.call(rbind, summaries)

  cat("\n== Backend certification summary ==\n")
  print(summary, row.names = FALSE)

  if (!is.null(parsed$summary) && nzchar(parsed$summary)) {
    dir.create(dirname(parsed$summary), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(summary, parsed$summary, row.names = FALSE)
  }

  failing <- summary$failed > 0L | summary$errors > 0L
  skipped <- summary$skipped > 0L
  if (any(failing) || (!isTRUE(parsed$allow_skips) && any(skipped))) {
    if (any(failing)) {
      bad <- summary$backend[failing]
      cat(sprintf("Backend certification failures: %s\n", paste(bad, collapse = ", ")))
    }
    if (!isTRUE(parsed$allow_skips) && any(skipped)) {
      bad <- summary$backend[skipped]
      cat(sprintf("Backend certification skips are not allowed: %s\n", paste(bad, collapse = ", ")))
    }
    quit(save = "no", status = 1L)
  }

  invisible(0L)
}

if (identical(Sys.getenv("AMATRIX_BACKEND_CERTIFICATION_AUTORUN", unset = "1"), "1")) {
  backend_certification_main()
}
