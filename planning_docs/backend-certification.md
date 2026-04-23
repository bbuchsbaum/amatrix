# Backend Certification Snapshot

This is the beta-readiness evidence ledger for optional accelerator backends.
It records what the current tests prove, not what the package aspires to prove.

Date: 2026-04-23

## Current Claim Matrix

| Backend | Current tier | Claimed beta surface | Evidence | Caveats |
|---|---|---|---|---|
| `cpu` | Authoritative | Full reference path, strict precision, Matrix/base-R oracle. | `devtools::test()` core conformance and coverage table in `planning_docs/quality-tracking.md`. | CPU failures are stop-ship defects. |
| `mlx` | Provisional supported | Apple Silicon fast path; `mode = "fast"` auto-selects MLX when healthy; dense matmul at calibrated sizes; fast Cholesky solve; broad dense/model conformance. | `tests/testthat/test-backend-certification-mlx.R`; `tests/testthat/test-cross-backend-conformance.R`; `tests/testthat/test-chol-factor.R`; `tests/testthat/test-backend-integration.R`. | Avoid plain `Rscript file.R` for MLX probe workers; use top-level `Rscript -e 'source(...)'`. Exact/strict paths may route to CPU by design. |
| `arrayfire` | Provisional supported, explicit probe | Dense fast backend when explicitly probed; cross-backend dense/model conformance; dense matmul layout regression fixed. | `tests/testthat/test-regression-arrayfire-matmul-layout.R`; ArrayFire blocks in `tests/testthat/test-cross-backend-conformance.R`; `tests/testthat/test-regression-arrayfire-worker-crash.R`. | Probe is explicit; Apple Silicon runtime defaults to ArrayFire CPU backend unless overridden. GPU/OpenCL ArrayFire paths remain crash-sensitive and are not the default claim. |
| `opencl` | Experimental | Explicit opt-in backend; dense matmul supported subset; rsvd benchmark accuracy fixed when OpenCL is installed and probed. | `tools/benchmark-regression.R`; `tests/testthat/test-opencl-model-core.R`; `tests/testthat/test-opencl-eigen.R`; `amatrix-x3c.1.1`. | Requires `AMATRIX_OPENCL_PROBE_GPU=1`. Not installed in the latest local conformance run, so beta tier remains experimental until fixed-hardware evidence is current. |
| `metal` | Experimental | Sparse/product-path work only. | Tracked under `amatrix-x3c.1`. | Not a general dense fast path for beta. |

## Latest Local Gates

```bash
Rscript --vanilla -e 'devtools::test(filter = "backend-certification-mlx|cross-backend-conformance|arrayfire-matmul-layout|backend-health|backend-integration", stop_on_failure = FALSE)'
```

Result on 2026-04-23: `281` pass, `0` fail, `0` skip.

```bash
Rscript --vanilla -e 'devtools::test(filter = "arrayfire", stop_on_failure = FALSE)'
```

Result on 2026-04-23: `11` pass, `0` fail, `0` skip.

## Release Gate Rule

Do not promote a backend tier from this file into `README.md` unless the
corresponding test evidence is green on hardware where that backend is actually
available. A backend that needs an explicit probe must say so in user-facing
docs and must route away cleanly when unprobed or unhealthy.
