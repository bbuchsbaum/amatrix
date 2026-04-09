# amatrix Quality Tracking: Accuracy & Performance

This document is the authoritative methodology reference for verifying that
amatrix is numerically correct and fast. Read it before adding a new op, before
interpreting a test failure, and before comparing benchmark numbers.

---

## 1. Two Pillars

| Pillar | Mechanism | Runs in CI? |
|--------|-----------|-------------|
| **Accuracy** | testthat conformance suite | Yes — every `devtools::test()` |
| **Performance** | `tools/benchmark-regression.R` | No — run manually by developer |

Every exported op must appear in both pillars. The coverage table in §4 tracks
current status.

---

## 2. Accuracy Methodology

### 2.1 Principle

Every exported op has a named base-R reference function. The test asserts:

```
max|am_result − ref| / max(|ref|, 1e-12)  <  tol
```

`expect_equal(..., tolerance = tol)` in testthat implements this as a relative
comparison.

### 2.2 Tolerances

| Context | Tolerance | Rationale |
|---------|-----------|-----------|
| GPU backends (MLX, ArrayFire) | `1e-4` | Float32 arithmetic; a GEMM over p columns accumulates O(p) ULPs. 1e-4 is conservative but not loose. |
| CPU backend | `1e-10` | Native float64; no precision gap. |
| Factorisation output (eigenvectors, singular vectors) | residual-based | Vectors are sign/rotation-ambiguous. Test `‖A·V − V·Λ‖_F / ‖A‖_F` instead. Threshold: `sqrt(tol)`. |
| Statistical ops (`am_covariance`) | same as above | Float32 covariance deviates ≈1e-6 vs `cov()` on typical data; 1e-4 GPU tol comfortably covers this. |
| `am_dist` vs `dist()` | `max(tol, 1e-6)` | `am_dist` uses the GEMM identity (‖x‖²+‖y‖²−2x·y); `dist()` uses direct element-wise sum. Both are float64-correct but differ by ≈1e-6, so CPU tolerance is floored at 1e-6 rather than 1e-10. |
| Algorithm convergence (am_rsvd, block Lanczos) | reconstruction error | `‖A − Â‖_F / ‖A‖_F ≤ 10 × optimal + 0.15` |

### 2.3 Test Organisation

```
tests/testthat/
  test-cross-backend-conformance.R   # primary per-op × per-backend harness
  test-conformance.R                 # basic dense/sparse alignment (1e-10)
  test-chol-factor.R                 # Cholesky factor lifecycle
  test-svd-factor.R                  # SVD factor, PCA, projection
  test-block-lanczos.R               # Block Lanczos convergence
  test-backend-integration.R         # QR decomposition across backends
  test-dispatch-primitives.R         # S4 dispatch correctness
  test-residency-lifecycle.R         # GPU residency management
  test-gemm-dist.R                   # .am_gemm, tiled am_dist
```

The cross-backend harness runs `.run_backend_conformance(backend, tol)` for every
available backend. GPU backends are auto-skipped when not installed.

### 2.4 Adding a New Op

1. Add a check inside `.run_backend_conformance()` with a clearly named
   `label = tag("<op_name>")`.
2. Name the base-R reference in the comment immediately above the check.
3. Update the coverage table in §4.
4. If the op is algorithm-level (iterative, non-deterministic), add a dedicated
   `test_that(...)` block and test reconstruction/residual, not raw values.

### 2.5 Interpreting Failures

- **GPU failure, CPU pass:** likely a float32 precision issue. Check if tol=1e-4
  is actually tight enough for this op at this problem size.
- **CPU failure:** genuine bug. The CPU backend is pure float64; numerical
  disagreement with base R at 1e-10 is always a defect.
- **Both fail:** API mismatch (wrong reference function, wrong `dimnames`, etc.)
  or a constructor/dispatch regression.

---

## 3. Performance Methodology

### 3.1 Protocol

| Parameter | Value |
|-----------|-------|
| Reps | 10 timed iterations (bench::mark default warmup) |
| Metric | `median_ms` (median wall time in milliseconds) |
| Comparison | `speedup_vs_cpu = cpu_ms / backend_ms` |
| Tool | `bench::mark(iterations=10, check=FALSE, memory=FALSE)` |

**Standard suites:**

- `dense`: `matmul`, `crossprod`, `covariance`, `dist`, `chol`, `solve_rhs`,
  `eigen_sym`, `many_lm`, `rsvd`, `sinkhorn`
- `sparse`: `spmv`, `spmm`, `block_lanczos`, `svd_factor_subspace`

Backends are probed systematically from the checkout:
- `cpu`
- `mlx`
- `opencl`
- `arrayfire` when explicitly enabled for benchmark runs
- `metal` for sparse products only

Each backend/op family runs in its own `Rscript` worker. That isolates hard
backend failures: a segfault or native abort is recorded as a `crash` incident
in the output instead of taking down the entire benchmark session.

For MLX spectral benchmarks, prefer fresh one-shot `Rscript -e` cells per
operation rather than combined file-entry workers. See
`docs/mlx-spectral-benchmark-instability.md`.

**Standard sizes:**

| Label | n | p | Notes |
|-------|---|---|-------|
| small | 256 | 32 | Cold-start / threshold regime |
| medium | 1024 | 128 | Main workload |
| large | 4096 | 128 | Memory-pressure regime |

Op-specific caps:
- `am_dist` input is capped at 512 rows to prevent OOM even at `large` size.
- `am_sinkhorn` uses square positive inputs at `128x128`, `512x512`, and
  `1024x1024` with a fixed 25-iteration loop so the resident iterative path is
  represented without dominating total runtime.

`tools/benchmark-regression.R` also auto-discovers repo-local optional backend
builds from `.tmp/opencl-lib`, `.tmp/lib`, `.tmp/backends-lib`, and
`.tmp/metal-lib`, so MLX, OpenCL, and Metal can be benchmarked directly from
this checkout when those backend packages are present.

### 3.2 Baseline File

```
tools/baseline.csv                     <- machine-local baseline snapshot
tools/benchmark-results/<timestamp>/   <- per-run artifacts
```

- The baseline stores canonical key columns plus `median_ms` for successful
  rows only.
- Each run writes:
  - `raw-results.csv` with one row per benchmark cell
  - `summary.csv` with CPU-relative speedups and baseline ratios
  - `incidents.csv` with `error`, `crash`, and `unavailable` rows
  - `metadata.rds` with host and runtime metadata
- **Baseline is machine-specific.** Numbers from a MacBook M3 are not comparable
  to numbers from an AWS g4dn instance. Regenerate after hardware changes.

### 3.3 Workflow

```bash
# First run on a machine (or after hardware change / major refactor):
Rscript tools/benchmark-regression.R --update

# Subsequent runs — compare to saved baseline:
Rscript tools/benchmark-regression.R

# Focus a run while debugging:
Rscript tools/benchmark-regression.R --backends=cpu,mlx --suites=sparse
```

The script prints a regression summary, an incident summary, and the artifact
paths for the run. Regressions are reported with
`ratio_vs_baseline = current_ms / baseline_ms`.

### 3.4 Regression Threshold

**>20% slower** than baseline triggers a manual review. This is intentionally
loose to absorb OS scheduling noise. A genuine regression (algorithm change,
extra allocation, routing regression) typically appears as >50%.

Do not update the baseline to mask a regression. Fix the code first.

Cold-path numbers include wrap/materialization cost by constructing fresh
`adgeMatrix`/`adgCMatrix` inputs inside the timed closure. Warm dense rows reuse
resident dense inputs; warm sparse product rows use the resident sparse product
path and are recorded as `resident` variants. Sparse iterative rows
(`block_lanczos`, `svd_factor_subspace`) use `warm` to mean the sparse lhs is
pre-bound before the timed loop so the compiled product-plan path can reuse the
resident operand across repeated products inside the algorithm.

### 3.5 Interpreting Speedup

`speedup_vs_cpu` is the ratio of CPU median time to this backend's median time.
Values < 1.0 mean the backend is *slower* than CPU (common at small sizes due to
PCIe/Metal overhead). This is expected and not a bug — the routing thresholds in
`docs/backend-benchmarks.md` exist precisely to avoid GPU dispatch below the
crossover point.

---

## 4. Coverage Table

Legend: **tested** = in cross-backend harness or dedicated suite | **partial** = tested but not vs canonical base R | **—** = not yet benchmarked

| Op | Base-R Reference | Accuracy | Benchmark |
|----|-----------------|----------|-----------|
| `%*%` / `am_matmul` | `base::%*%` | tested | tested |
| sparse `%*%` (`SpMV` / `SpMM` via `adgCMatrix`) | `Matrix::%*%` | tested (dispatch + calibration + sparse suites) | tested |
| `crossprod` / `am_crossprod` | `base::crossprod` | tested | tested |
| `tcrossprod` / `am_tcrossprod` | `base::tcrossprod` | tested | tested |
| `+`, `*`, scalar ewise | `+`, `*` | tested | tested |
| `rowSums` / `am_row_sums` | `base::rowSums` | tested | tested |
| `colSums` / `am_col_sums` | `base::colSums` | tested | tested |
| `chol` / `am_chol` | `base::chol` | tested (conformance + factor suites) | tested |
| `solve` / `am_solve` | `base::solve` | tested (conformance + backend-local smoke) | tested |
| `am_rsvd` | truncated SVD | tested (recon + ortho) | tested |
| `am_dist` | `as.matrix(dist())` | tested | tested |
| `am_kernel` linear/rbf/poly | manual formula | tested | — |
| `am_covariance` | `cov()` | tested | tested |
| `am_eigen` | `eigen()` / `eigh()` | tested (values + residual) | tested |
| `svd()` dispatch | `base::svd` | tested (values only) | tested |
| `am_many_lm` | `lm.fit()` loop | tested | tested |
| `am_block_lanczos` / `am_block_svd` | `La.svd` | own suite | tested |
| `am_sinkhorn` | plain Sinkhorn loop | own suite | tested |
| `am_irlba` | `am_rsvd` | own suite | — |
| `am_chol_solve` / `am_chol_factor` | `solve(chol(...))` | own suite | — |
| `am_svd_factor` + project/reconstruct | reconstruction | own suite | partial |
| `am_qr` + qr.* methods | `base::qr.*` | own suite | tested |
| `am_ridge_fit` | `lm.fit` + penalty | partial | — |
| `am_wls_fit` | weighted `lm.fit` | partial | — |
| `am_pca_coef` | `prcomp` loadings | partial | — |
| `am_diag` | `base::diag` | partial | — |
| `t()` / `am_transpose` | `base::t` | tested (dispatch) | — |
| `am_ewise` | `+`, `-`, `*`, `/` | tested | — |
| `amatrix_bind_resident` | resident key present + backend reuse | tested | — |
| `amatrix_resident_backend_for` | backend choice for resident hot path | tested | — |
| `amatrix_prepare_operands` | auto-bind repeated product operands | tested | — |
| `amatrix_compile_product` | reusable sparse/dense product plan | tested | — |

---

## 5. How to Run

### Run full accuracy suite
```r
devtools::test()
# or just the cross-backend harness:
devtools::test(filter = "cross-backend-conformance")
```

### Capture a fresh performance baseline
```bash
Rscript tools/benchmark-regression.R --update
```

### Check for performance regressions
```bash
Rscript tools/benchmark-regression.R
```

### Run the install/load smoke gate
Use this after touching `NAMESPACE`, exports, S3/S4 registrations, or package
startup/registration code.
```bash
Rscript tools/smoke-install-load.R
```

### Run a specific subsystem benchmark (deep-dive)
```bash
Rscript tools/benchmark-dense-products.R
Rscript tools/benchmark-opencl-qr.R
Rscript tools/benchmark-sinkhorn.R
Rscript tools/benchmark-flagship-many-y.R
# ... see tools/ for the full list
```

---

## 6. Keeping This Document Current

- When a new op is added to NAMESPACE: add a row to §4 with status `missing`.
- When a test is added: update §4 to `tested` or `partial`.
- When a benchmark is added: update §4 benchmark column.
- When tolerances change: update §2.2 and add a one-line rationale.

The doc is prose, not code — it does not auto-update. Treat a stale coverage
table as a test failure.
