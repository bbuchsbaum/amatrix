#!/usr/bin/env Rscript
# tools/benchmark-regression.R
#
# Canonical performance regression harness for amatrix.
# Covers all major op families at three standard sizes.
#
# Usage:
#   Rscript tools/benchmark-regression.R           # compare to baseline
#   Rscript tools/benchmark-regression.R --update  # write/overwrite baseline
#
# Requires: bench, amatrix (amatrix.mlx optional but recommended on Apple Silicon)
#
# Baseline file: tools/baseline.csv (machine-local, not committed to git).
# Regenerate after hardware changes or major refactors with --update.

suppressPackageStartupMessages({
  library(amatrix)
})

if (!requireNamespace("bench", quietly = TRUE)) {
  stop("bench package required: install.packages('bench')", call. = FALSE)
}
library(bench)

update_mode <- "--update" %in% commandArgs(trailingOnly = TRUE)

mlx_ok <- tryCatch(
  requireNamespace("amatrix.mlx", quietly = TRUE) &&
    amatrix.mlx::amatrix_mlx_is_available(),
  error = function(e) FALSE
)
BACKENDS <- c("cpu", if (mlx_ok) "mlx")

# Standard sizes: small / medium / large.
# p is kept at 32 for wide ops (dist, many_lm) to keep runtime manageable.
SIZES <- list(
  small  = list(n = 256L,  p = 32L),
  medium = list(n = 1024L, p = 128L),
  large  = list(n = 4096L, p = 128L)
)

run_size <- function(sz, backend) {
  set.seed(1L)
  X <- matrix(rnorm(sz$n * sz$p), sz$n, sz$p)
  set.seed(2L)
  Y <- matrix(rnorm(sz$n * 3L), sz$n, 3L)

  aX <- adgeMatrix(X, preferred_backend = backend)

  # Cap dist/tcrossprod input to 512 rows to stay under GPU OOM at large sizes.
  n_dist <- min(sz$n, 512L)
  Xs     <- X[seq_len(n_dist), , drop = FALSE]

  res <- bench::mark(
    matmul     = { invisible(aX %*% matrix(rnorm(sz$p * sz$p), sz$p)) },
    crossprod  = { invisible(crossprod(aX)) },
    covariance = { invisible(am_covariance(aX)) },
    dist       = { invisible(am_dist(Xs)) },
    many_lm    = { invisible(am_many_lm(aX, Y, method = "qr", cache = FALSE)) },
    rsvd       = { invisible(am_rsvd(aX, k = 10L)) },
    iterations = 10L,
    check      = FALSE,
    memory     = FALSE,
    time_unit  = "ms"
  )

  res$backend <- backend
  res$size    <- paste0(sz$n, "x", sz$p)
  res
}

message("Running benchmark suite (", paste(BACKENDS, collapse = ", "),
        " × ", length(SIZES), " sizes × 10 reps) ...")

results <- do.call(rbind, lapply(SIZES, function(sz) {
  do.call(rbind, lapply(BACKENDS, function(b) {
    message("  ", b, " @ ", sz$n, "x", sz$p, " ...")
    run_size(sz, b)
  }))
}))

# Build tidy output frame
out <- data.frame(
  op        = as.character(results$expression),
  size      = results$size,
  backend   = results$backend,
  median_ms = as.numeric(results$median) * 1e3,
  stringsAsFactors = FALSE
)

# Compute speedup_vs_cpu
cpu_med <- out[out$backend == "cpu", c("op", "size", "median_ms")]
names(cpu_med)[3] <- "cpu_ms"
out <- merge(out, cpu_med, by = c("op", "size"), all.x = TRUE)
out$speedup_vs_cpu <- round(out$cpu_ms / out$median_ms, 2L)
out$cpu_ms <- NULL

out <- out[order(out$size, out$op, out$backend), ]

BASELINE <- file.path(dirname(sys.frame(0)$ofile %||% "."), "baseline.csv")
if (is.na(BASELINE) || !nzchar(BASELINE)) BASELINE <- "tools/baseline.csv"

if (update_mode || !file.exists(BASELINE)) {
  # Write a small header comment then the CSV
  hdr <- sprintf(
    "# amatrix benchmark baseline — %s\n# R %s on %s\n",
    format(Sys.time(), "%Y-%m-%d %H:%M"),
    paste0(R.version$major, ".", R.version$minor),
    Sys.info()[["nodename"]]
  )
  cat(hdr, file = BASELINE)
  write.table(out[, c("op", "size", "backend", "median_ms", "speedup_vs_cpu")],
              BASELINE, sep = ",", col.names = TRUE, row.names = FALSE,
              append = TRUE)
  message("Baseline written to ", BASELINE)
  message("\nResults:")
  print(out[, c("op", "size", "backend", "median_ms", "speedup_vs_cpu")],
        row.names = FALSE)
} else {
  # Read baseline (skip comment lines starting with #)
  base_lines <- readLines(BASELINE)
  base_data  <- read.csv(textConnection(
    paste(base_lines[!grepl("^#", base_lines)], collapse = "\n")))
  names(base_data)[names(base_data) == "median_ms"] <- "base_ms"

  comp <- merge(out, base_data[, c("op", "size", "backend", "base_ms")],
                by = c("op", "size", "backend"), all.x = TRUE)
  comp$ratio <- round(comp$median_ms / comp$base_ms, 3L)

  regs <- comp[!is.na(comp$ratio) & comp$ratio > 1.2, ]
  if (nrow(regs) > 0L) {
    message("\nREGRESSIONS DETECTED (>20% slower than baseline):")
    print(regs[, c("op", "size", "backend", "median_ms", "base_ms", "ratio")],
          row.names = FALSE)
  } else {
    message("\nOK — no regressions vs baseline.")
  }

  message("\nFull comparison:")
  print(comp[order(comp$size, comp$op, comp$backend),
             c("op", "size", "backend", "median_ms", "speedup_vs_cpu", "ratio",
               "base_ms")],
        row.names = FALSE)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
