#!/usr/bin/env Rscript

r_string_literal <- function(x) encodeString(x, quote = "\"")

usage <- function() {
  cat(
    paste(
      "ArrayFire runtime diagnostic",
      "",
      "Usage:",
      "  Rscript tools/debug-arrayfire-gpu.R [--runtime=cpu|opencl|cuda|oneapi] [--n=128]",
      "",
      "Examples:",
      "  Rscript tools/debug-arrayfire-gpu.R --runtime=cpu",
      "  Rscript tools/debug-arrayfire-gpu.R --runtime=opencl --n=64",
      sep = "\n"
    ),
    "\n"
  )
}

parse_args <- function(args) {
  out <- list(runtime = "", n = 128L, help = FALSE)
  for (arg in args) {
    if (arg %in% c("--help", "-h")) {
      out$help <- TRUE
      next
    }
    if (startsWith(arg, "--runtime=")) {
      out$runtime <- sub("^--runtime=", "", arg)
      next
    }
    if (startsWith(arg, "--n=")) {
      out$n <- as.integer(sub("^--n=", "", arg))
    }
  }
  out
}

locate_repo_root <- function() {
  raw_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", raw_args, value = TRUE)
  if (length(file_arg) == 0L) {
    return(normalizePath(".", winslash = "/", mustWork = TRUE))
  }

  script_path <- sub("^--file=", "", file_arg[[1L]])
  normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
}

launch_child <- function(repo_root, runtime, n) {
  runtime_clause <- if (nzchar(runtime)) {
    sprintf("Sys.setenv(AMATRIX_ARRAYFIRE_BACKEND=%s);", r_string_literal(runtime))
  } else {
    "Sys.unsetenv('AMATRIX_ARRAYFIRE_BACKEND');"
  }

  child_code <- paste0(
    "setwd(", r_string_literal(repo_root), ");",
    "source('tools/benchmark-helpers.R', local=globalenv());",
    "Sys.setenv(AMATRIX_BENCHMARK_ARRAYFIRE='1', AMATRIX_BENCHMARK_ARRAYFIRE_UNSAFE='1', AMATRIX_ARRAYFIRE_PROBE_GPU='1');",
    runtime_clause,
    "spec <- .benchmark_optional_backend_specs(include_arrayfire=TRUE)[['arrayfire']];",
    "stopifnot(!is.null(spec));",
    "stopifnot(.benchmark_enable_backend(spec));",
    "ns <- ensure_optional_backend_namespace('amatrix.arrayfire', repo_dir='backends/amatrix.arrayfire');",
    "diag <- get('amatrix_arrayfire_diagnostics', envir=ns, inherits=FALSE)();",
    "af_matmul <- get('amatrix_arrayfire_matmul', envir=ns, inherits=FALSE);",
    "af_crossprod <- get('amatrix_arrayfire_crossprod', envir=ns, inherits=FALSE);",
    "af_svd <- get('amatrix_arrayfire_svd', envir=ns, inherits=FALSE);",
    "cat(sprintf('ACTIVE_BACKEND=%s\\n', diag$active_backend));",
    "cat(sprintf('DEVICE_COUNT=%s\\n', diag$device_count));",
    "cat(sprintf('AVAILABLE_BACKENDS=%s\\n', diag$available_backends));",
    "cat(sprintf('LAPACK_AVAILABLE=%s\\n', diag$lapack_available));",
    sprintf("set.seed(1); x <- matrix(rnorm(%dL * %dL), nrow=%dL, ncol=%dL);", n, n, n, n),
    "y <- x;",
    "mm <- af_matmul(x, y);",
    "cat(sprintf('DIRECT_MATMUL=%s\\n', paste(dim(mm), collapse='x')));",
    "cp <- af_crossprod(x);",
    "cat(sprintf('DIRECT_CROSSPROD=%s\\n', paste(dim(cp), collapse='x')));",
    "sv <- af_svd(x, nu=5L, nv=5L);",
    "cat(sprintf('DIRECT_SVD_D=%s\\n', length(sv$d)));",
    "load_benchmark_amatrix();",
    "ax <- amatrix::adgeMatrix(x, preferred_backend='arrayfire', precision='fast');",
    "ay <- amatrix::adgeMatrix(y, preferred_backend='arrayfire', precision='fast');",
    "plan <- amatrix::amatrix_backend_plan(ax, 'matmul', y=ay);",
    "cat(sprintf('DISPATCH=%s|%s\\n', plan$chosen, plan$chosen_path));",
    "prod <- ax %*% ay;",
    "cat(sprintf('DISPATCH_MATMUL=%s\\n', paste(dim(prod), collapse='x')));"
  )

  warned_status <- NULL
  quoted_args <- vapply(c("-e", child_code), shQuote, character(1), USE.NAMES = FALSE)
  output <- withCallingHandlers(
    system2(
      file.path(R.home("bin"), "Rscript"),
      quoted_args,
      stdout = TRUE,
      stderr = TRUE
    ),
    warning = function(w) {
      warned_status <<- attr(w, "status") %||% warned_status
      invokeRestart("muffleWarning")
    }
  )
  status <- attr(output, "status") %||% warned_status %||% 0L

  list(status = status, output = output)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage()
    return(invisible(0L))
  }

  repo_root <- locate_repo_root()
  launch <- launch_child(repo_root, args$runtime, args$n)

  cat("ArrayFire diagnostic\n")
  cat(sprintf("repo_root: %s\n", repo_root))
  cat(sprintf("requested_runtime: %s\n", if (nzchar(args$runtime)) args$runtime else "<default>"))
  cat(sprintf("child_status: %s\n\n", launch$status))
  if (length(launch$output) > 0L) {
    cat(paste(launch$output, collapse = "\n"), sep = "\n")
    if (!grepl("\n$", paste(launch$output, collapse = "\n"))) {
      cat("\n")
    }
  }

  invisible(launch$status)
}

status <- main()
quit(save = "no", status = status)
