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
- `opencl`
- `metal` for sparse products only
- `arrayfire` only when explicitly enabled for benchmark runs and not hard-gated by platform safety policy

`mlx` is intentionally excluded from the canonical per-group worker harness by
default because isolated worker launch is still unstable on Apple Silicon.
Benchmark MLX with dedicated top-level runners such as
`tools/benchmark-mlx-native-rsvd.R`, or explicitly re-enable worker-mode MLX
only for crash-probing via `AMATRIX_BENCHMARK_MLX_WORKERS=1`.

Each backend/op family runs in its own `Rscript` worker. That isolates hard
backend failures: a segfault or native abort is recorded as a `crash` incident
in the output instead of taking down the entire benchmark session.

For MLX spectral benchmarks, prefer fresh one-shot `Rscript -e` cells per
operation rather than combined file-entry workers. See
`planning_docs/mlx-spectral-benchmark-instability.md`.

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
`planning_docs/backend-benchmarks.md` exist precisely to avoid GPU dispatch below
the crossover point.

---

## 4. Coverage Matrix

Every exported *operation* (excluding management / introspection APIs listed
in §4.1) must appear in this matrix. For release, every row must have all four
test-type columns populated and a benchmark row. This matrix is the authoritative
schema; the PR-gate test `tests/testthat/test-coverage-table.R` reads `NAMESPACE`
and fails if an op is missing a row, and — when `AMATRIX_COVERAGE_STRICT=1` — if
any row has a gap.

**Legend:**

- `✓` present and verified
- `○` partial (present but not at full rigor — e.g. Oracle test exists but uses a non-canonical reference, or only a subset of the input space is covered)
- `—` gap: will be closed by Tracks 3 (tests) / 4 (benchmarks) / 5 (tiers)

**Status columns:**

| Column | Meaning |
|---|---|
| Oracle | Differential test against a named base-R / Matrix / reference function (§2.1) |
| Metamorphic | Algebraic-invariant / property test (`chol(A) %*% t(chol(A)) == A`, etc.). See Track 3 task 2. |
| Adversarial | Evil-input coverage: empty, 1×1, rank-deficient, `Inf`/`NaN`, extreme scale/conditioning. See Track 3 task 3. |
| Regression | Dedicated minimal repro for any bug ever filed against this op. See Track 6 task 1. |
| Benchmark | Row in `tools/benchmark-regression.R`. |
| Tiers | Per-backend status. Format: `C<status> M<status> A<status> O<status> X<status>` where C=CPU, M=MLX, A=ArrayFire, O=OpenCL, X=Metal. See §8 for tier definitions and current backend tier. |

> **Current honest state (2026-04-13):** Most ops have Oracle coverage via the
> cross-backend conformance harness and a benchmark row, but Metamorphic /
> Adversarial / Regression are almost entirely gaps — that is Track 3's mandate.
> Per-backend tier columns use `?` where the backend is registered but its tier
> under the op is not yet proven by the gate. See §8.

### 4.1 Excluded from coverage (management / introspection / constructors)

These exports are not linear-algebra operations. They are covered by dedicated
test files and are intentionally excluded from the coverage matrix. New exports
matching these patterns are exempt by the coverage-table test:

- `amatrix_*` — policy, calibration, residency, dispatch, memory, warmup, execution, explain, register, backend-query, cache (`test-backend-integration.R`, `test-residency-lifecycle.R`, `test-dispatch-*.R`)
- `with_amatrix` — scope helper
- `adgeMatrix`, `adgCMatrix`, `as_adgeMatrix`, `as_adgCMatrix` — constructors / coercions (`test-constructors.R`)
- `resident_handle`, `as.matrix.resident_handle`, `ncol.resident_handle`, `nrow.resident_handle`, `rh_rowSums`, `rh_colSums` — resident-handle helpers (`test-bind-resident.R`, `test-residency-lifecycle.R`)
- `kron_matrix`, `as.matrix.KronMatrix` — KronMatrix wrappers (the op `kron` is covered below)

### 4.2 Operations matrix

#### 4.2.1 Dense products and GEMM family

| Export | Oracle | Metamorphic | Adversarial | Regression | Benchmark | Tiers |
|---|---|---|---|---|---|---|
| `matmul` / `%*%` method | ✓ | ✓ | ✓ | — | ✓ | C✓ M? A? O? X— |
| `gemm` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `addmm` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `dot` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `crossprod` method | ✓ | ✓ | ✓ | — | ✓ | C✓ M? A? O? X— |
| `tcrossprod` method | ✓ | ✓ | ✓ | — | ✓ | C✓ M? A? O? X— |
| `crossprod_add_diag` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `crossprod_weighted` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `tcrossprod_weighted` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `xty_weighted` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `quad_form` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `kron` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `sparse %*%` (adgCMatrix method) | ✓ | — | — | — | ✓ | C✓ M? A? O? X— |

#### 4.2.2 Reductions, scaling, element-wise

| Export | Oracle | Metamorphic | Adversarial | Regression | Benchmark | Tiers |
|---|---|---|---|---|---|---|
| `rowsums` / `rowSums` method | ✓ | ✓ | ✓ | — | ✓ | C✓ M? A? O? X— |
| `colsums` / `colSums` method | ✓ | ✓ | ✓ | — | ✓ | C✓ M? A? O? X— |
| `rowmeans` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `colmeans` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `rowscale` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `colscale` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `sym` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `trace` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `trace_estim` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `ewise` | ✓ | ✓ | — | — | ✓ | C✓ M? A? O? X— |
| `am_ewise_inplace` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `am_sweep` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `am_sweep_inplace` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `segment_sum` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `segment_mean` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `am_scatter_mean` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `am_rowargmax` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `am_rowargmin` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `am_colargmax` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `am_colargmin` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `t()` method | ○ | ✓ | — | — | — | C✓ M? A? O? X— |

#### 4.2.3 Solvers and factorizations

| Export | Oracle | Metamorphic | Adversarial | Regression | Benchmark | Tiers |
|---|---|---|---|---|---|---|
| `solve` method | ✓ | ○ | ○ | — | ✓ | C✓ M? A? O? X— |
| `solve_triangular` | ○ | ✓ | — | — | — | C✓ M? A? O? X— |
| `chol` method | ✓ | — | — | — | ✓ | C✓ M? A? O? X— |
| `chol_factor` | ✓ | ✓ | ✓ | ✓ | — | C✓ M? A? O? X— |
| `chol_solve` | ✓ | — | — | — | — | C✓ M? A? O? X— |
| `chol_solve_batches` | ✓ | — | — | — | — | C✓ M? A? O? X— |
| `chol_diag` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `chol_logdet` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `batch_chol` | ○ | ✓ | — | — | — | C✓ M? A? O? X— |
| `batch_solve` | ○ | ✓ | — | — | — | C✓ M? A? O? X— |
| `batch_crossprod` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `lu_factor` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `lu_solve` | ○ | ✓ | — | — | — | C✓ M? A? O? X— |
| `am_qr` | ✓ | ✓ | — | ✓ | ✓ | C✓ M? A? O? X— |
| `qr_info` | ○ | ✓ | — | — | — | C✓ M? A? O? X— |
| `qr_downdate` | ✓ | — | ✓ | — | — | C✓ M? A? O? X— |

#### 4.2.4 Spectral and iterative

| Export | Oracle | Metamorphic | Adversarial | Regression | Benchmark | Tiers |
|---|---|---|---|---|---|---|
| `svd` method | ✓ | — | — | — | ✓ | C✓ M? A? O? X— |
| `svd_factor` | ✓ | — | — | — | ○ | C✓ M? A? O? X— |
| `svd_project` | ✓ | — | — | — | — | C✓ M? A? O? X— |
| `svd_reconstruct` | ✓ | — | — | — | — | C✓ M? A? O? X— |
| `eigh` | ✓ | — | — | — | ✓ | C✓ M? A? O? X— |
| `rsvd` | ✓ | — | — | — | ✓ | C✓ M? A? O? X— |
| `irlba` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `irlba_native` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `block_lanczos` | ✓ | ✓ | ✓ | — | ✓ | C✓ M? A? O? X— |
| `block_svd` | ✓ | — | — | — | ✓ | C✓ M? A? O? X— |
| `mat_sqrt` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `mat_pow` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `mat_log` | ○ | — | — | — | — | C✓ M? A? O? X— |

#### 4.2.5 Statistical models

| Export | Oracle | Metamorphic | Adversarial | Regression | Benchmark | Tiers |
|---|---|---|---|---|---|---|
| `lm_fit` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `many_lm` | ✓ | ✓ | — | — | ✓ | C✓ M? A? O? X— |
| `array_lm` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `lm_loo_cv` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `wls_fit` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `ridge_fit` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `ridge_path` | ○ | ✓ | ✓ | — | — | C✓ M? A? O? X— |
| `pca_coef` | ○ | — | — | — | — | C✓ M? A? O? X— |

### 4.2.7 Composed workflow tests

Workflow tests chain multiple exported ops in a user-facing pattern that
matches the "plain R" idiom. They are the highest-value cross-backend
conformance tests because they verify that a shared factor or resident
tensor survives intact across multiple consumer ops. The ops are already
covered individually in §4.2.x; these rows document the end-to-end
chain test.

| Workflow | Ops chained | Test location | Backends covered |
|---|---|---|---|
| Gaussian Process regression (fit + predict + log marginal likelihood) | `kernel_matrix` → `chol_factor` → `chol_solve` (reused) → `matmul` → `chol_logdet` + `quad_form` | `tests/testthat/test-cross-backend-conformance.R::.run_gp_conformance` (and per-backend `test_that` blocks) | cpu (1e-10), mlx (1e-4), arrayfire (1e-4), opencl (1e-4 when installed) |

The GP conformance runner asserts: posterior mean vs base-R reference,
posterior variance vs base-R reference, log marginal likelihood vs
reference, `K α = y` metamorphic invariant, and non-negative predictive
variance. The reference is a 15-line plain-R GP via
`chol` / `backsolve` / `forwardsolve`.

#### 4.2.6 Distances, kernels, covariance, sinkhorn, Woodbury

| Export | Oracle | Metamorphic | Adversarial | Regression | Benchmark | Tiers |
|---|---|---|---|---|---|---|
| `dist_matrix` | ✓ | ✓ | — | — | ✓ | C✓ M? A? O? X— |
| `pairwise_sqdist_argmin` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `kernel_matrix` | ✓ | ✓ | ✓ | — | — | C✓ M? A? O? X— |
| `covariance` | ✓ | ✓ | — | — | ✓ | C✓ M? A? O? X— |
| `correlation` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `sinkhorn` | ✓ | — | — | — | ✓ | C✓ M? A? O? X— |
| `woodbury_solve` | ○ | — | — | — | — | C✓ M? A? O? X— |
| `woodbury_logdet` | ○ | — | — | — | — | C✓ M? A? O? X— |

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

- When a new op is added to NAMESPACE: add a row to §4 with gaps marked `—`.
  The PR-gate test `tests/testthat/test-coverage-table.R` will fail until the
  row exists.
- When a test is added: update the relevant cell in §4 from `—` / `○` to `✓`.
- When a benchmark is added: update §4 Benchmark column.
- When a backend earns or loses a tier for an op: update the Tiers column for
  that op. Tier labels in §8 are also updated.
- When tolerances change: update §2.2 and add a one-line rationale.

The doc is prose, not code — it does not auto-update. Treat a stale coverage
table as a test failure.

---

## 7. Stop-Ship Rules

The release gate (Track 2) blocks shipping if any of the following are true.
These are hard rules, not guidelines. Overriding requires a written, signed-off
justification recorded in `cran-comments.md`.

1. **CPU correctness failure.** Any test in the cross-backend conformance suite
   that fails on the `cpu` backend. CPU is the authoritative reference (§8);
   a CPU failure is always a genuine numerical defect.
2. **Supported-tier backend crash.** Any backend listed as `supported` in §8
   segfaults, aborts, or produces a non-finite result on any op at any of the
   three size classes (small / main / xlarge). Experimental-tier backends are
   exempt from this rule but still must not leak state into later tests.
3. **Unexplained performance regression.** Any op × backend × size row in
   `tools/benchmark-regression.R` whose `median_ms` is > 20% slower than
   `tools/baseline.csv` on release hardware **without** a written explanation
   landed in the same commit as a baseline update. Baseline updates that silently
   absorb a regression are a stop-ship rule violation.
4. **Missing coverage row.** Any entry in `NAMESPACE` that qualifies as an
   operation (per §4.1 exclusions) and is absent from the §4 matrix. This is
   mechanically enforced by `test-coverage-table.R` on every PR.
5. **Orphan repro outside the active suite.** Any file under
   `tests/testthat/_problems/` at release time. Every known failing case must
   either be fixed and promoted into the active suite (as `test-regression-*.R`)
   or explicitly deleted with a commit message explaining the decision.
6. **Residency tripwire trip.** Any host-copy of a GPU-resident tensor during a
   conformance run with `options(amatrix.residency.tripwire = TRUE)`. A trip
   means S4 dispatch fell through to a base-R generic and silently copied the
   tensor home — which is a correctness-visible performance bug.
7. **Non-empty fallback log after a clean conformance run.** `amatrix_fallback_log()`
   must be empty after `devtools::test()` completes. A non-empty log means a
   backend claimed to support an op it cannot actually execute.
8. **Backend tier label drift.** The tier labels in §8 and `README.md` must be
   generated from gate evidence (Track 5 task 6), not hand-edited. A PR that
   hand-edits the tier section without corresponding gate evidence is rejected.

---

## 8. Backend Tiers (Honest Assessment)

The release makes explicit support claims per backend. A backend only claims
a tier it has earned, and the label is generated from gate evidence (7 green
nightlies + benchmark net-benefit on its intended workloads), not from
aspiration.

### 8.1 Tier definitions

| Tier | User-visible claim | Gate requirement |
|---|---|---|
| **Authoritative** | Reference of record for correctness. Always available. | Passes Correctness Contract (§2) at the tightest tolerance (`1e-10`) for every op. CPU is the only authoritative backend and is the numerical oracle for all others. |
| **Supported** | First-class: listed in README as fast path. Conformance-clean, crash-free, benchmarked as net-beneficial. | Green on every Correctness cell for the op × size classes it claims, zero crashes across 7 consecutive nightly runs, benchmark evidence that the backend is ≥ 0.9× CPU (§3.4 floor) on its intended workloads, auto-fallback (Track 5 task 4) empty after a clean conformance run. |
| **Experimental** | Available, opt-in, warns on first use. Not polished, not blocked. | Registered and loadable. README explicitly labels as experimental. Health probe (Track 5 task 1) may mark it unhealthy; in that case the dispatcher routes away silently. |

### 8.2 Current tier (2026-04-13, pre-release)

These tiers are **honest provisional assessments**, not final release claims.
The release gate (Track 2) will re-evaluate against real nightly evidence and
the README tier matrix will be regenerated from that evidence (Track 5 task 6).

| Backend | Tier | Basis |
|---|---|---|
| `cpu` | **Authoritative** | Base R / Matrix BLAS. Pure float64. CPU failures are genuine defects. Reference for all Oracle tests. |
| `mlx` | **Provisional Supported** | Primary fast path on Apple Silicon. Conformance-green on macOS. **Caveat:** multi-worker isolation is unstable (see `planning_docs/mlx-spectral-benchmark-instability.md`); default `amatrix.mlx.workers = 1`. Benchmark fast-path status needs re-measurement under Track 4 gates before the tier is confirmed. |
| `arrayfire` | **Provisional Supported** | Portable GPU option. Has known backend-specific gotchas (non-resident `tcrossprod` requires `nrow(X) == nrow(Y)`; see `feedback_af_tcrossprod` memory). Needs crash-free 7-nightly run before the tier is confirmed. |
| `opencl` | **Experimental** | Registered via `amatrix.opencl`. Limited op coverage, limited hardware coverage in CI. Opt-in via `getOption("amatrix.optional_backends", TRUE)`. |
| `metal` | **Experimental** | Sparse-product-only path via `amatrix.metal`. Limited scope; not a general fast path. |

> **Aspiration vs enforcement.** The aspiration is that all five backends
> reach Supported by release. The enforcement is that a backend ships in the
> tier the evidence proves at gate time. Anything labelled Provisional above
> will either be promoted to Supported or demoted to Experimental by the
> release gate — without apology.

### 8.3 What the tier buys you

- **Authoritative backend**: never skipped, always run, always the oracle.
- **Supported backend**: tested at every gate (PR + nightly), documented in
  README, appears in the fast-path routing table, subject to the speed floor
  (§3.4 and Track 4 task 2). A regression here is stop-ship (§7).
- **Experimental backend**: exists, may work, may not. Health probe gates
  routing. README explicitly says so. Not subject to §7 crash rule. Users who
  care about it opt in explicitly.

---

## 9. Machinery

The following artifacts enforce the contracts in this document:

- `tests/testthat/test-coverage-table.R` — reads `NAMESPACE` and §4, enforces
  row existence. With `AMATRIX_COVERAGE_STRICT=1`, also enforces no gaps.
- `tools/audit-dispatch.R` (Track 3 task 7) — enforces dispatch hermeticity.
- `tools/benchmark-regression.R` — the performance harness; the §3.4 threshold
  and §3 baseline discipline are the teeth.
- `amatrix_backend_status()`, `amatrix_fallback_log()` (Track 5) — runtime
  telemetry that the conformance suite asserts on.
- `.github/workflows/R-CMD-check.yaml` (PR gate), `nightly-stress.yaml` (nightly
  gate), `release-gate.yaml` (release gate) — the CI lanes (Track 2).
