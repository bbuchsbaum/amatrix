#!/usr/bin/env Rscript
# tools/benchmark-regression.R
#
# Canonical performance regression harness for amatrix.
# Covers all major op families at three standard sizes.
#
# Two timing variants per (backend × size):
#   cold  — adgeMatrix wraps host data; first op triggers lazy GPU upload.
#            Represents real-world first-call latency.
#   warm  — data is pre-resident on the GPU before timing starts.
#            Shows pure compute throughput, removing PCIe/Metal transfer cost.
#
# Usage:
#   Rscript tools/benchmark-regression.R             # compare to baseline
#   Rscript tools/benchmark-regression.R --update    # write/overwrite baseline
#   Rscript tools/benchmark-regression.R --warm-only # show warm results only
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

args        <- commandArgs(trailingOnly = TRUE)
update_mode <- "--update"    %in% args
warm_only   <- "--warm-only" %in% args

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

  # Cap dist input to 512 rows to stay under GPU OOM at large sizes.
  n_dist <- min(sz$n, 512L)
  Xs     <- X[seq_len(n_dist), , drop = FALSE]

  # Pre-generate RHS so the same matrix is reused across iterations.
  set.seed(3L)
  B_host <- matrix(rnorm(sz$p * sz$p), sz$p, sz$p)

  # ── Cold-start timing ──────────────────────────────────────────────────────
  # adgeMatrix is freshly constructed; first operation triggers lazy GPU upload.
  aX_cold <- adgeMatrix(X, preferred_backend = backend)

  cold <- bench::mark(
    matmul     = { invisible(aX_cold %*% B_host) },
    crossprod  = { invisible(crossprod(aX_cold)) },
    covariance = { invisible(am_covariance(aX_cold)) },
    dist       = { invisible(am_dist(Xs)) },
    many_lm    = { invisible(am_many_lm(aX_cold, Y, method = "qr", cache = FALSE)) },
    rsvd       = { invisible(am_rsvd(aX_cold, k = 10L)) },
    iterations = 10L,
    check      = FALSE,
    memory     = FALSE
  )
  cold$variant <- "cold"

  # ── Warm-start timing ──────────────────────────────────────────────────────
  # Pre-upload aX and the RHS to GPU; dry-run each op to prime caches.
  aX_warm <- adgeMatrix(X, preferred_backend = backend)
  aB_warm <- adgeMatrix(B_host, preferred_backend = backend)
  if (backend != "cpu") {
    invisible(aX_warm %*% aB_warm)  # trigger X upload + warm metal kernel
    invisible(crossprod(aX_warm))   # ensure X resident
  }

  warm <- bench::mark(
    matmul     = { invisible(aX_warm %*% aB_warm) },
    crossprod  = { invisible(crossprod(aX_warm)) },
    covariance = { invisible(am_covariance(aX_warm)) },
    dist       = { invisible(am_dist(Xs)) },
    many_lm    = { invisible(am_many_lm(aX_warm, Y, method = "qr", cache = TRUE)) },
    rsvd       = { invisible(am_rsvd(aX_warm, k = 10L)) },
    iterations = 10L,
    check      = FALSE,
    memory     = FALSE
  )
  warm$variant <- "warm"

  res          <- rbind(cold, warm)
  res$backend  <- backend
  res$size     <- paste0(sz$n, "x", sz$p)
  res
}

message("Running benchmark suite (", paste(BACKENDS, collapse = ", "),
        " \u00d7 ", length(SIZES), " sizes \u00d7 2 variants \u00d7 10 reps) ...")

results <- do.call(rbind, lapply(SIZES, function(sz) {
  do.call(rbind, lapply(BACKENDS, function(b) {
    message("  ", b, " @ ", sz$n, "x", sz$p, " ...")
    run_size(sz, b)
  }))
}))

# Build tidy output frame.
# as.numeric(bench_time) returns seconds; * 1e3 → milliseconds.
out <- data.frame(
  op        = as.character(results$expression),
  size      = results$size,
  backend   = results$backend,
  variant   = results$variant,
  median_ms = as.numeric(results$median) * 1e3,
  stringsAsFactors = FALSE
)

# Compute speedup_vs_cpu within each (variant × size × op).
cpu_med <- out[out$backend == "cpu", c("op", "size", "variant", "median_ms")]
names(cpu_med)[4] <- "cpu_ms"
out <- merge(out, cpu_med, by = c("op", "size", "variant"), all.x = TRUE)
out$speedup_vs_cpu <- round(out$cpu_ms / out$median_ms, 2L)
out$cpu_ms <- NULL

out <- out[order(out$variant, out$size, out$op, out$backend), ]

BASELINE <- "tools/baseline.csv"

if (update_mode || !file.exists(BASELINE)) {
  cols <- c("op", "size", "backend", "variant", "median_ms", "speedup_vs_cpu")
  write.csv(out[, cols], BASELINE, row.names = FALSE)
  message("Baseline written to ", BASELINE)
  message("\nResults:")
  print(out[, cols], row.names = FALSE)
} else {
  base_data <- read.csv(BASELINE)
  names(base_data)[names(base_data) == "median_ms"] <- "base_ms"

  merge_keys <- intersect(c("op", "size", "backend", "variant"),
                          names(base_data))
  comp <- merge(out, base_data[, c(merge_keys, "base_ms")],
                by = merge_keys, all.x = TRUE)
  comp$ratio <- round(comp$median_ms / comp$base_ms, 3L)

  # Regression check on warm-path (representative of steady-state perf).
  warm_comp <- comp[comp$variant == "warm", ]
  regs <- warm_comp[!is.na(warm_comp$ratio) & warm_comp$ratio > 1.2, ]
  if (nrow(regs) > 0L) {
    message("\nREGRESSIONS DETECTED (warm path >20% slower than baseline):")
    print(regs[, c("op", "size", "backend", "median_ms", "base_ms", "ratio")],
          row.names = FALSE)
  } else {
    message("\nOK \u2014 no regressions vs baseline (warm path).")
  }

  disp <- if (warm_only) comp[comp$variant == "warm", ] else comp
  message("\n", if (warm_only) "Warm" else "Full", " comparison:")
  print(disp[order(disp$variant, disp$size, disp$op, disp$backend),
             c("op", "size", "backend", "variant", "median_ms",
               "speedup_vs_cpu", "ratio", "base_ms")],
        row.names = FALSE)
}
