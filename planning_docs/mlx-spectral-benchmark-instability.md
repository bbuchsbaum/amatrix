# MLX Spectral Benchmark Instability

Date: 2026-04-09

## RESOLVED 2026-07-01 (empirically, on current MLX)

The direct `Rscript file.R` NSException crash no longer reproduces on this
machine with Homebrew mlx 0.31.1 / mlx-c 0.6.0 (originally observed on an
older MLX; upstream ml-explore/mlx#2691 was closed as not-reproducible, so
there is no explicit fix commit to cite).

Evidence: 26 consecutive direct `Rscript file.R` runs — each performing
Metal device init via `amatrix_mlx_native_available_bridge` (the
historically crashing call) followed by GPU matmul/crossprod/svd loops —
completed with exit 0 and correct results. The repeatable certification
harness lives at [tools/certify-mlx-file-entry.R](/Users/bbuchsbaum/code/amatrix/tools/certify-mlx-file-entry.R)
(20/20 pass on 2026-07-01) and is wired into the nightly stress workflow.

Consequences:
- The file-entry probe guard (`.amatrix_mlx_direct_file_entry`) and the
  env-var probe gate are retired; MLX probing is default-on with
  `AMATRIX_MLX_PROBE_GPU=0` / `options(amatrix.auto_probe = FALSE)` as
  opt-outs.
- The first Metal probe of a session runs in a disposable child process
  (crash containment), so machines where the old crash still exists fall
  back to CPU with an informative reason instead of aborting the session.
- The `Rscript -e` benchmark-worker policy below is retained only as
  historical context for the *multi-step spectral worker* aborts, which
  were a separate (also unreproduced since) failure shape; the nightly
  certification run guards against regression.

The remainder of this document is retained as historical context.

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

Single-operation MLX spectral calls in a fresh top-level `Rscript -e` process
are the most reliable path right now. Nested subprocesses launched from the
general SVD harness are still unsafe for native MLX spectral work.

Observed on this machine:

- one-shot MLX `svd()` in `Rscript -e`: works
- one-shot MLX `rsvd()` in `Rscript -e`: works
- combined benchmark-worker session for MLX spectral timing: can abort with
  `NSException`
- sourcing the full SVD harness before native MLX resident initialization can
  also abort with `NSException`
- the standalone native RSVD runner at
  [tools/benchmark-mlx-native-rsvd.R](/Users/bbuchsbaum/code/amatrix/tools/benchmark-mlx-native-rsvd.R)
  avoids that initialization shape and records usable native MLX RSVD timings

## Harness Policy

The SVD backend harness should therefore:

- avoid direct `Rscript file.R` MLX workers
- keep exact MLX `svd` on safe CPU fallback because `mlx_linalg_svd` is
  CPU-stream-only in the current MLX bridge
- keep safe CPU fallback as the default for MLX `rsvd` inside the general SVD
  harness
- use `--mlx-native-spectral` only as a crash-probe mode for isolated MLX RSVD
  workers
- use `Rscript -e 'setwd("..."); source("tools/benchmark-mlx-native-rsvd.R")'`
  for actual native MLX RSVD timing
- keep crash isolation at the row level so failed MLX cells become incident
  rows instead of taking down the full benchmark run

This is a runtime/workflow workaround, not a fix for the underlying MLX issue.
