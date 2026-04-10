args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

script_path <- if (length(file_arg) >= 1L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path("tools", "smoke-install-load.R"), winslash = "/", mustWork = TRUE)
}

pkg_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
r_bin <- file.path(R.home("bin"), "R")
rscript_bin <- file.path(R.home("bin"), "Rscript")
`%||%` <- function(x, y) if (is.null(x)) y else x

run_checked <- function(cmd, args, env = character(), label) {
  out <- tryCatch(
    suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = TRUE, env = env)),
    error = function(e) {
      stop(sprintf("%s failed to launch: %s", label, conditionMessage(e)), call. = FALSE)
    }
  )
  status <- attr(out, "status") %||% 0L
  if (!identical(status, 0L)) {
    if (length(out)) {
      cat(paste(out, collapse = "\n"), file = stderr())
      cat("\n", file = stderr())
    }
    stop(sprintf("%s failed with exit status %s", label, status), call. = FALSE)
  }
  invisible(out)
}

smoke_lib <- tempfile("amatrix-smoke-lib-")
dir.create(smoke_lib, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Installing amatrix into temporary library: %s\n", smoke_lib))
run_checked(
  r_bin,
  c("CMD", "INSTALL", "-l", smoke_lib, pkg_dir),
  label = "R CMD INSTALL"
)

smoke_file <- tempfile("amatrix-smoke-load-", fileext = ".R")
writeLines(
  c(
    "lib_dir <- Sys.getenv('AMATRIX_SMOKE_LIB')",
    "if (!nzchar(lib_dir)) stop('AMATRIX_SMOKE_LIB is not set', call. = FALSE)",
    ".libPaths(c(lib_dir, .libPaths()))",
    "suppressPackageStartupMessages(library(amatrix))",
    "",
    "exports <- getNamespaceExports('amatrix')",
    "required_exports <- c(",
    "  'adgeMatrix', 'adgCMatrix', 'as_adgeMatrix', 'resident_handle',",
    "  'rh_rowSums', 'rh_colSums', 'rowmeans', 'colmeans', 'matmul', 'am_qr',",
    "  'sinkhorn'",
    ")",
    "missing_exports <- setdiff(required_exports, exports)",
    "if (length(missing_exports)) {",
    "  stop(sprintf('Missing exports: %s', paste(missing_exports, collapse = ', ')), call. = FALSE)",
    "}",
    "",
    "required_s3 <- list(",
    "  c('as.matrix', 'resident_handle'),",
    "  c('dim', 'resident_handle'),",
    "  c('print', 'resident_handle')",
    ")",
    "for (sig in required_s3) {",
    "  method <- utils::getS3method(sig[[1L]], sig[[2L]], optional = TRUE)",
    "  if (is.null(method)) {",
    "    stop(sprintf('Missing S3 method: %s.%s', sig[[1L]], sig[[2L]]), call. = FALSE)",
    "  }",
    "}",
    "required_namespace_symbols <- c(",
    "  'resident_handle', 'as.matrix.resident_handle', 'dim.resident_handle',",
    "  'nrow.resident_handle', 'ncol.resident_handle', 'print.resident_handle'",
    ")",
    "missing_symbols <- required_namespace_symbols[!vapply(required_namespace_symbols, exists, logical(1), envir = asNamespace('amatrix'), inherits = FALSE)]",
    "if (length(missing_symbols)) {",
    "  stop(sprintf('Missing namespace symbols: %s', paste(missing_symbols, collapse = ', ')), call. = FALSE)",
    "}",
    "",
    "required_s4 <- list(",
    "  c('rowMeans', 'adgeMatrix'),",
    "  c('colMeans', 'adgeMatrix'),",
    "  c('rowMeans', 'adgCMatrix'),",
    "  c('colMeans', 'adgCMatrix'),",
    "  c('qr', 'adgeMatrix'),",
    "  c('qr', 'adgCMatrix')",
    ")",
    "for (sig in required_s4) {",
    "  if (!methods::hasMethod(sig[[1L]], sig[[2L]])) {",
    "    stop(sprintf('Missing S4 method: %s(%s)', sig[[1L]], sig[[2L]]), call. = FALSE)",
    "  }",
    "}",
    "",
    "x_host <- matrix(c(1, 2, 3, 4), nrow = 2L)",
    "dense <- adgeMatrix(x_host)",
    "sparse <- adgCMatrix(matrix(c(1, 0, 0, 1), nrow = 2L))",
    "",
    "stopifnot(inherits(dense, 'adgeMatrix'))",
    "stopifnot(inherits(sparse, 'adgCMatrix'))",
    "stopifnot(isTRUE(all.equal(rowMeans(dense), base::rowMeans(x_host))))",
    "stopifnot(isTRUE(all.equal(colMeans(dense), base::colMeans(x_host))))",
    "stopifnot(isTRUE(all.equal(rowMeans(sparse), base::rowMeans(as.matrix(sparse)))))",
    "stopifnot(isTRUE(all.equal(colMeans(sparse), base::colMeans(as.matrix(sparse)))))",
    "stopifnot(inherits(qr(dense), 'amQR'))",
    "stopifnot(inherits(qr(sparse), 'amQR'))",
    "",
    "cat('Install/load smoke OK\\n')"
  ),
  con = smoke_file
)

cat("Running installed-package smoke checks in a fresh R session\n")
run_checked(
  rscript_bin,
  c(smoke_file),
  env = c(paste0("AMATRIX_SMOKE_LIB=", smoke_lib)),
  label = "installed-package smoke"
)

cat("Smoke install/load check passed\n")
