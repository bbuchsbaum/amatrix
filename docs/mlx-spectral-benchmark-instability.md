# MLX Spectral Benchmark Instability

Date: 2026-04-09

## Summary

MLX spectral benchmarks are not stable under every scripted launch pattern on this
Apple Silicon setup.

Two distinct issues were observed:

- Direct `Rscript file.R` launch is a known unsafe mode for MLX Metal probing.
- Even under safer `Rscript -e '...; source(...)'` launch, combined multi-step
  spectral benchmark workers can still abort with an uncaught `NSException`.

## Evidence

- [backends/amatrix.mlx/R/backend.R](/Users/bbuchsbaum/code/amatrix/backends/amatrix.mlx/R/backend.R) documents that the Metal GPU probe is safe in `Rscript -e` / interactive / testthat contexts, but not in plain `Rscript file.R`.
- [backends/amatrix.mlx/src/amatrix_mlx_matmul.c](/Users/bbuchsbaum/code/amatrix/backends/amatrix.mlx/src/amatrix_mlx_matmul.c) documents the same `NSRangeException` risk during Metal device initialization.
- [tools/benchmark-qr-runtime.R](/Users/bbuchsbaum/code/amatrix/tools/benchmark-qr-runtime.R) and [tools/print-svd-factor-calibration.R](/Users/bbuchsbaum/code/amatrix/tools/print-svd-factor-calibration.R) already route MLX benchmarking through `Rscript -e` or one-shot commands for this reason.

## Current Working Assumption

Single-operation MLX spectral calls in a fresh safe-mode process are the most
reliable path right now.

Observed on this machine:

- one-shot MLX `svd()` in `Rscript -e`: works
- one-shot MLX `rsvd()` in `Rscript -e`: works
- combined benchmark-worker session for MLX spectral timing: can abort with
  `NSException`

## Harness Policy

The SVD backend harness should therefore:

- avoid direct `Rscript file.R` MLX workers
- run each MLX spectral cell in its own fresh `Rscript -e` process
- keep crash isolation at the row level so failed MLX cells become incident
  rows instead of taking down the full benchmark run

This is a runtime/workflow workaround, not a fix for the underlying MLX issue.
