# Certify that MLX Metal probing and GPU ops are stable under direct
# `Rscript file.R` launch — the entry mode historically guarded against
# because Metal device init threw NSRangeException on older MLX
# (see planning_docs/mlx-spectral-benchmark-instability.md, upstream
# ml-explore/mlx#2691).
#
# Run this file DIRECTLY as `Rscript tools/certify-mlx-file-entry.R`
# (that is the point). Drive N iterations via the wrapper loop:
#   for i in $(seq 1 20); do Rscript tools/certify-mlx-file-entry.R || break; done
#
# Exit status 0 with final line "CERTIFY-OK" = this process survived
# Metal init plus matmul/crossprod/svd on the MLX path.

entry <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(entry) == 0L) {
  stop("certify-mlx-file-entry.R must be launched as `Rscript file.R` — ",
       "the direct file-entry mode is what is being certified")
}

Sys.setenv(AMATRIX_MLX_PROBE_GPU = "1")
if (!requireNamespace("amatrix.mlx", quietly = TRUE)) {
  cat("CERTIFY-SKIP: amatrix.mlx not installed\n")
  quit(status = 0L)
}

avail <- .Call("amatrix_mlx_native_available_bridge", PACKAGE = "amatrix.mlx")
if (!isTRUE(avail)) {
  cat("CERTIFY-SKIP: MLX native bridge reports unavailable\n")
  quit(status = 0L)
}

suppressMessages({
  library(amatrix)
  library(amatrix.mlx)
})

set.seed(42)
x <- matrix(rnorm(512 * 256), 512, 256)
A <- adgeMatrix(x, mode = "fast")

for (i in seq_len(3L)) {
  B <- A %*% t(A)
  cp <- crossprod(A)
  s <- svd(A)
  ref_cp <- crossprod(x)
  stopifnot(
    is.finite(sum(as.matrix(B))),
    max(abs(as.matrix(cp) - ref_cp)) < 1e-3 * max(abs(ref_cp)),
    all(is.finite(s$d)),
    abs(s$d[1] - svd(x, nu = 0, nv = 0)$d[1]) < 1e-2
  )
}

cat("CERTIFY-OK\n")
