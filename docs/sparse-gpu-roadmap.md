# Sparse GPU Roadmap For amatrix

## Goal

Turn sparse GPU support into a narrow, benchmark-backed extension of the
existing backend contract rather than a second product center.

The current repo already has the basic shape for this:

- `adgCMatrix` routes `%*%` through `matmul()`
- backend features already include `sparse_spmm`
- MLX and ArrayFire already expose sparse-dense multiply hooks
- sparse calibration already distinguishes `spmv` from `spmm`
- CPU and `Matrix` remain the semantic fallback

What is missing is not a greenfield sparse design. What is missing is a clear
product stance about which sparse GPU paths are worth hardening and which ones
should remain CPU-first.

## Current State

The sparse posture in the repo today is:

- Public sparse objects are still `dgCMatrix`-compatible `adgCMatrix` values.
- Sparse semantics are still CPU-first in the backend contract.
- Optional backends may accelerate sparse-dense products opportunistically.
- Sparse factorizations, indexing, and structural sparse operations still
  primarily belong to CPU `Matrix` paths.

In practical terms, the codebase already contains:

- cold sparse-dense multiply hooks in MLX and ArrayFire
- sparse resident-store hooks for backends that can cache sparse operands
- resident sparse-dense multiply wrappers
- calibration buckets for `spmv` and `spmm`
- CPU fallback coverage for sparse-sparse multiply and sparse factorizations

That means the roadmap should start from hardening `SpMM`, not from searching
for a general-purpose sparse library first.

## Cross-Backend Pathway That Already Exists

Yes, but it is intentionally narrow.

There is already a backend-spanning pathway for sparse-dense product execution:

1. `adgCMatrix` methods route `%*%`, `crossprod()`, and `tcrossprod()` through
   the ordinary `matmul` / `crossprod` / `tcrossprod` wrappers.
2. Backend planning uses the same `supports(op, x, y)` contract and the same
   calibration machinery that dense products use, including separate `spmv`
   versus `spmm` thresholds.
3. Backends that want sparse acceleration can opt in by:
   - accepting `adgCMatrix` in `supports()`
   - advertising `sparse_spmm`
   - implementing sparse-dense product entry points
   - optionally implementing sparse resident storage plus `spmm_resident`
4. If no optional backend qualifies, CPU `Matrix` remains the fallback.

In other words, the pathway that spans backends today is:

`adgCMatrix` object model -> common product wrappers -> common backend planner
and calibration -> optional backend-specific sparse-dense kernel ->
common CPU fallback.

That is a real cross-backend pathway.

What does not yet span backends as a coherent product surface:

- sparse factorizations
- sparse spectral decompositions
- sparse-sparse GPU multiply as a benchmark-backed default
- a fully documented sparse resident chain across every optional backend

So the right statement is:

> We already have a cross-backend sparse product pathway. We do not yet have a
> full cross-backend sparse linear algebra product.

## Product Rule

Preserve this rule:

> Sparse is CPU-native first; GPU sparse is opportunistic and must justify
> itself with benchmark-backed wins on real workloads.

Corollaries:

- Do not widen the public sparse contract before the existing `SpMM` path is
  trustworthy.
- Do not promise sparse GPU factorization in v1.
- Do not treat sparse-sparse multiply as a default GPU path.
- Do not let portable-but-mediocre sparse GPU coverage outrank platform-native
  high-performance backends where those exist.

## Priority Order

### 1. Harden `SpMM` first

Target surface:

- `adgCMatrix %*% matrix`
- `adgCMatrix %*% adgeMatrix`
- `crossprod(adgCMatrix, dense)`
- `tcrossprod(adgCMatrix, dense)` when it can be implemented as sparse-dense
  multiply against a transposed dense RHS

Why first:

- This is already partially implemented.
- It matches the repo's dense-first plus many-RHS philosophy.
- Upload cost amortizes better with multi-column dense RHS than with vectors.
- It unlocks the most plausible real wins for iterative/model workloads.

What "done" means:

- backend planning chooses GPU only when sparse structure and RHS width justify it
- MLX and ArrayFire parity is covered by tests that compare against `Matrix`
- benchmarks show a real crossover against CPU on representative sparse shapes
- resident sparse reuse works reliably for repeated products

### 2. Productize `SpMV` second

Target surface:

- `adgCMatrix %*% numeric`
- `adgCMatrix %*% matrix(ncol = 1)`

Why second:

- The codebase already calibrates `spmv` separately, so the planner can make
  this distinction cleanly.
- `SpMV` is critical for Krylov, PageRank, and graph-style workloads.
- GPUs often underperform on small or poorly reused `SpMV`; this must stay
  benchmark-gated.

What "done" means:

- sparse vector products use the same planner and sparse residency machinery as
  `SpMM`
- heuristics prefer resident sparse operands or large enough `nnz`
- at least one real iterative workflow shows a credible benefit

### 3. Treat `SpGEMM` as optional until a workload forces it

Target surface:

- `adgCMatrix %*% adgCMatrix`

Why not earlier:

- Memory growth is harder to control.
- Performance depends heavily on sparsity pattern and output fill-in.
- Correctness and structural expectations are harder than sparse-dense multiply.
- The current CPU fallback already preserves semantics.

What would justify it:

- a benchmarked graph or sparse algebra workflow in this repo genuinely needs it
- a backend offers a mature sparse-sparse multiply path
- output sparsity semantics are specified and tested clearly

Until then, sparse-sparse multiply should remain a correctness-first CPU path.

### 4. Keep sparse factorizations on CPU for now

Remain CPU-first:

- `chol(adgCMatrix)`
- `qr(adgCMatrix)`
- `solve(adgCMatrix, ...)`
- sparse `svd` / `eigen`

Reason:

- `Matrix`/SuiteSparse already define the stable semantic baseline.
- Factorization surfaces are much larger than multiply surfaces.
- A partial GPU factorization story is harder to test and easier to oversell.

The first sparse GPU win should be products, not sparse direct solvers.

## Backend Strategy

### Preferred hierarchy

If the goal is maximal performance rather than maximum checkbox coverage, the
backend order should be:

1. platform-native high-performance sparse backend where available
2. existing MLX / ArrayFire sparse-dense paths where they win
3. CPU `Matrix` fallback everywhere else

More concretely:

- NVIDIA: a future `cuSPARSE`-style backend is the strongest candidate for peak
  sparse GPU performance
- AMD: a future `rocSPARSE`-style backend is the analogous path
- Apple: keep sparse GPU opportunistic and benchmark-led; do not force the
  dense MLX story to become a broad sparse promise without evidence
- portable OpenCL-style options: acceptable for experiments, not the default
  strategic answer for "maximal power and speed"

### Contract shape

No new public sparse class hierarchy is needed for the first pass.

The existing contract is already sufficient if backends do the following:

- advertise `sparse_spmm` in `features()`
- answer `supports(op = "matmul", x = adgCMatrix, y = dense)` honestly
- optionally support sparse resident storage plus `spmm_resident`
- reject unsupported sparse shapes predictably so CPU wins as fallback

## Benchmark Requirements

Sparse GPU work should only advance when all three of these are true:

1. Numeric parity is stable.
2. Planner thresholds are stable.
3. The benchmark win appears on more than one sparsity family.

Minimum benchmark matrix families:

- uniform random sparsity
- banded sparsity
- block-diagonal or clustered sparsity
- graph-like power-law sparsity

Minimum benchmark slices:

- `spmv` with one RHS column
- `spmm` with several RHS widths
- cold path versus resident sparse reuse
- density buckets that show where CPU should still win

## Implementation Checklist

### Phase A. Freeze the sparse GPU product boundary

- Keep the backend contract language aligned with the current sparse-first CPU
  rule.
- Document existing sparse GPU coverage as opportunistic `SpMM`, not general
  sparse acceleration.
- Make sure README and planning docs do not imply sparse GPU factorization.

### Phase B. Harden current `SpMM`

- Audit MLX and ArrayFire sparse-dense paths against the same shape matrix.
- Add focused parity tests for:
  - `adgCMatrix %*% dense`
  - `crossprod(adgCMatrix, dense)`
  - `tcrossprod(adgCMatrix, dense)`
  - resident sparse reuse when available
- Confirm planner behavior around:
  - low `nnz`
  - semi-dense sparse matrices
  - narrow versus wide RHS

### Phase C. Promote `SpMV` deliberately

- Add explicit benchmark slices that justify separate `spmv` thresholds.
- Verify that vector and one-column matrix inputs take the same intended path.
- Validate at least one iterative workload that reuses a resident sparse
  operand.

### Phase D. Defer `SpGEMM` unless needed

- Keep sparse-sparse multiply correctness tests.
- Do not introduce backend-specific sparse-sparse kernels until a workload
  requires them.
- When revisited, define expected output class and structural semantics first.

### Phase E. Revisit backend inventory only after Phase B/C data exists

- If MLX or ArrayFire already deliver the needed sparse-dense win, keep scope
  local and harden what exists.
- If they do not, evaluate whether a dedicated sparse backend package is
  justified.
- Do not broaden the backend matrix before the benchmark case is clear.

## Acceptance Criteria

The sparse GPU roadmap is succeeding when all of the following are true:

- `SpMM` is benchmark-backed on at least one real workload and several sparse
  matrix families.
- `SpMV` is separately calibrated and only selected when it actually wins.
- sparse-sparse multiply remains correctness-first unless and until a workload
  proves otherwise.
- docs, tests, and backend planning all describe the same sparse posture.
- CPU `Matrix` remains the semantic safety net for unsupported sparse work.

## Non-Goals

- Do not turn `amatrix` into a general sparse GPU framework.
- Do not promise backend-neutral sparse factorization parity in v1.
- Do not treat portable sparse GPU libraries as inherently preferable to
  platform-native high-performance backends.
- Do not count kernel existence as product readiness without calibration and
  workload-level benchmarks.
