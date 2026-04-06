# amatrix: Next Three Major Features

**Plan author:** planner
**Date:** 2026-04-05
**Scope:** amatrix core package + amatrix.mlx + amatrix.arrayfire backends
**Order:** A (Conformance Harness) -> B (GPU Cholesky + `am_lm_fit` GPU path) -> C (Hardened S3/S4 dispatch + ALTREP research stub)

---

## Context

amatrix is a GPU-accelerated R matrix library with a pluggable backend system. Core package (`amatrix`) defines `adgeMatrix` (dense) and `adgCMatrix` (sparse) S4 classes, a residency cache, and a dispatch pipeline (`amatrix_dispatch_op` + `_try_resident_*`). Two backends currently ship:

- `amatrix.arrayfire` — portable GPU via ArrayFire (f32 internal, row-major bridge).
- `amatrix.mlx` — Apple Silicon via MLX (f32 internal, row-major bridge). Preferred on macOS (39x vs 9x speedup over AF in benchmarks).

Existing public API (`NAMESPACE`) already exports: `am_matmul`, `am_crossprod`, `am_tcrossprod`, `am_solve`, `am_chol`, `am_qr`, `am_svd`, `am_rsvd`, `am_dist`, `am_kernel`, `am_lm_fit`, `am_ridge_fit`, `am_wls_fit`, `am_many_lm`, `am_array_lm`, `am_irlba`, `am_pca_coef`, `am_correlation`, `am_covariance`, etc.

**Known pain points that motivate this plan:**

1. `helper-conformance.R` has per-op scaffolding but no cross-backend matrix that runs the same operation on `cpu`, `mlx`, `arrayfire` and asserts numerical agreement with base R. `am_rsvd`, `am_dist`, `am_kernel`, `tcrossprod` are not covered.
2. `amatrix_mlx_solve_resident` and `amatrix_mlx_chol_resident` (in `backends/amatrix.mlx/R/backend.R:623,635`) are stubs that call `base::solve`/`base::chol` on the host. `am_ridge_fit`/`am_lm_fit` therefore do a host round-trip for the p×p solve — fine when p is small, catastrophic when the caller wanted true GPU residency for a million-voxel workflow.
3. Existing S4 dispatch on `adgeMatrix` covers the `Matrix` generics, but `base::%*%` is a primitive that bypasses S4 in certain LHS=`matrix`/RHS=`adgeMatrix` combinations and in packages that call `.Primitive("%*%")` directly. The manual `numeric %*% adgeMatrix` and `matrix %*% adgeMatrix` methods in `methods-dense.R:19-26` already document this friction (irlba's `mult(v, A)` pattern).

---

## Work Objectives

1. **Feature A:** Ship a cross-backend conformance harness that runs every amatrix op on every available backend and diffs against base R, plus a benchmark helper for the same matrix.
2. **Feature B:** Implement real C-level GPU Cholesky factor + solve in both MLX and AF backends, wire them into `solve_resident`/`chol_resident`, add `am_chol_solve(A, B)` public wrapper, and make `am_lm_fit(X, Y, lambda=0)` use a GPU-resident normal-equations path end-to-end (no host round-trip on the p×p factor/solve).
3. **Feature C:** Harden the S3/S4 dispatch surface so that `base::%*%`, `base::crossprod`, `base::tcrossprod` on `adgeMatrix` always reach the GPU path regardless of LHS/RHS class, and land a documentation-grade ALTREP research note describing what a true ALTREP implementation would require (scoped as a future major release item, not built now).

---

## Guardrails

### Must Have
- **Numerical agreement:** Conformance tolerance is `1e-5` for f32 backends (MLX/AF), `1e-10` for CPU/base. These thresholds are asserted per-op, per-backend.
- **Backend opt-in:** Tests must `skip_if_backend_package_missing` for MLX/AF; CI without GPUs must still pass.
- **No API breakage:** `am_lm_fit`'s existing signature (`intercept`, `include_fitted`, `include_residuals`, `cache`, `method`) is preserved. `lambda = 0` is added as a new parameter with a default that makes the existing OLS path unchanged.
- **Residency cache respected:** New `am_chol_solve` must use `_try_resident_*` helpers so a cached `adgeMatrix` LHS stays on-device across calls.
- **GPU Cholesky correctness:** Both MLX and AF bridges must detect non-SPD input and `Rf_error` with a clear message matching `base::chol`'s "not positive definite" wording.
- **f32 precision modes:** If caller uses `precision = "strict"`, Cholesky solve must refuse the MLX/AF path and fall through to CPU. `precision = "fast"` enables the f32 GPU path.

### Must NOT Have
- No new public API surface beyond `am_chol_solve` and the `lambda` argument.
- No changes to the `adgeMatrix`/`adgCMatrix` class definitions or slots.
- No true ALTREP C code in this plan — only a research note + hardened dispatch.
- No dependency on Matrix-internal symbols that aren't already imported in NAMESPACE.
- No changes to `am_rsvd` math (that was just fixed in commit 37a2918).
- No regeneration of existing `.so`/`.o` files as part of review — implementers build locally.

---

## Task Flow

```
Feature A (Conformance Harness)
  |
  v
Feature B (GPU Cholesky + am_lm_fit GPU path)
  |   ^-- Feature A is a dependency: new bridges need conformance coverage
  v
Feature C (Hardened dispatch + ALTREP research note)
      ^-- Feature A is a dependency: new dispatch paths need conformance coverage
```

---

## Detailed TODOs

### Feature A: Cross-Backend Conformance Harness

**Files to create:**
- `tests/testthat/helper-conformance-matrix.R` — new helper that enumerates `(backend, op, input_shape)` cases.
- `tests/testthat/test-conformance-backends.R` — new test file that iterates the matrix.
- `tests/testthat/test-conformance-numerics.R` — new test file for `am_dist`, `am_kernel`, `am_rsvd` numerical agreement.
- `inst/benchmarks/conformance-bench.R` — benchmark helper that runs the same matrix and prints timings (not a test).

**Files to modify:**
- `tests/testthat/helper-conformance.R` — add `available_backend_names()` utility, add `conformance_tolerance(backend)` that returns `1e-5` for mlx/arrayfire and `1e-10` for cpu.
- `DESCRIPTION` — add `bench` or `microbenchmark` to Suggests (prefer `bench` for memory tracking).

**Key implementation decisions:**

1. **Backend matrix shape:** A list of lists. Each row is `list(backend = "mlx", op = "matmul", shapes = list(c(64,64), c(512,32)), tol = 1e-5)`. `test-conformance-backends.R` walks this list and calls `testthat::test_that(sprintf("%s %s [%dx%d]", backend, op, nr, nc), ...)`.

2. **Op coverage:** `matmul`, `crossprod`, `tcrossprod` (both with `y=NULL` and `y=non-null`), `ewise` (for `+`, `-`, `*`, `/`, scalar and matrix RHS), `rowSums`, `colSums`, `solve`, `chol`, `svd`, `am_rsvd`, `am_dist` (three methods), `am_kernel` (five kernels).

3. **Input shapes per op:** small (2x2, 4x3) for correctness, one medium (128x64) to exercise the GEMM path without blowing up CI time. SVD/am_rsvd use tall (256x32) and wide (32x256) shapes. Non-SPD inputs excluded from `solve`/`chol`; use `z <- matrix(rnorm(n*n), n); crossprod(z) + diag(n)*0.5` per `test-conformance.R:50`.

4. **Expected-value computation:** Compute the ground truth **once** per shape on the CPU using base R or `Matrix`, cache it, then diff each backend's output against it. Avoids N backends x M ops x K shapes redundant CPU computations.

5. **Skip logic:** `skip_if_backend_package_missing(spec)` already exists. Wrap each backend block in `for (spec in optional_backend_specs()) { ... skip if missing ... }`. CPU is always included.

6. **Tolerance selection:** `conformance_tolerance(backend)` returns `1e-5` for f32 backends, `1e-10` for f64 CPU, `1e-4` for `am_rsvd` (randomized — needs looser bound), `1e-3` for `am_kernel` rbf/laplacian (exp() amplifies f32 error).

7. **Benchmark helper separation:** `inst/benchmarks/conformance-bench.R` is a **runnable script**, not a test. It uses `bench::mark()` to time each `(backend, op, shape)` cell and prints a tibble. Invoked via `Rscript inst/benchmarks/conformance-bench.R` in CI nightly, not on every PR.

**Acceptance criteria:**
- `devtools::test(filter = "conformance-backends")` passes on a dev machine with either MLX or AF available, and skips cleanly on a machine with neither.
- Every exported op in the table above has at least one shape tested per available backend.
- `Rscript inst/benchmarks/conformance-bench.R` produces a tibble with columns `backend`, `op`, `shape`, `median`, `mem_alloc`, `n_itr` and exits 0.
- When the harness is pointed at a deliberately-broken backend (e.g., using `make_recording_backend` with a corrupted `matmul` wrapper), it fails with a clear message identifying which `(backend, op, shape)` cell mismatched.

---

### Feature B: GPU Cholesky + `am_lm_fit` GPU Path

**Files to create:**
- `backends/amatrix.mlx/src/amatrix_mlx_cholesky.c` — new C file with MLX Cholesky + triangular solve kernels.
- `backends/amatrix.arrayfire/src/arrayfire_cholesky.c` — new C file (or section inside `arrayfire_bridge.c`) with AF Cholesky + solve kernels.
- `R/chol-solve.R` — new R file exporting `am_chol_solve(A, B, upper = TRUE)`.
- `tests/testthat/test-am-chol-solve.R` — unit tests for the new wrapper.
- `tests/testthat/test-am-lm-fit-gpu.R` — end-to-end test that `am_lm_fit(X, Y, lambda=0)` stays GPU-resident.

**Files to modify:**
- `backends/amatrix.mlx/src/init.c` — register `amatrix_mlx_chol_bridge`, `amatrix_mlx_chol_solve_bridge` Call entries.
- `backends/amatrix.mlx/R/backend.R` — replace `amatrix_mlx_solve_resident` and `amatrix_mlx_chol_resident` stubs (lines 623-638) with real bridge calls; add `amatrix_mlx_chol_solve_bridge` R wrapper; advertise `"chol"` and `"solve"` in `capabilities()`.
- `backends/amatrix.mlx/src/amatrix.mlx.so` — rebuilt by implementer.
- `backends/amatrix.arrayfire/src/init.c` — register AF chol/solve Call entries.
- `backends/amatrix.arrayfire/R/backend.R` — replace `amatrix_arrayfire_solve_resident` / `amatrix_arrayfire_chol_resident` (lines 196,204) with real bridge calls.
- `R/wrappers.R` — add `am_chol_solve()` dispatch function following the `am_solve` pattern (lines 415-448), and add `.amatrix_try_resident_chol_solve` helper after `.amatrix_try_resident_chol` (line 270-279).
- `R/models-lm.R` — `am_lm_fit` signature (line 609) gets a new `lambda = 0` parameter; internally route through `.amatrix_lm_core_gpu` when `lambda == 0` and `inherits(X, "adgeMatrix")`, else keep existing behavior. Add `.amatrix_lm_core_gpu` helper that uses `am_chol_solve(crossprod(X) + lambda*I, crossprod(X, Y))`.
- `NAMESPACE` — add `export(am_chol_solve)`.
- `R/backend-registry.R` — add `"chol_solve"` to the list of recognized capability names (documentation only, no enforcement).

**API design:**

```r
# New public function
am_chol_solve <- function(A, B, upper = TRUE) {
  # A: adgeMatrix or matrix (SPD, p x p)
  # B: adgeMatrix, matrix, or numeric vector (p x q or length p)
  # upper: whether to return the upper-triangular Cholesky factor convention
  #
  # Returns X such that A %*% X == B, via GPU Cholesky when backend supports it.
  # Dispatches: preferred_backend of A, falling through to CPU on unsupported precision.
}

# Extended signature (lambda is new, defaults to 0 = unchanged behavior)
am_lm_fit <- function(
  X, Y,
  lambda = 0,              # NEW: L2 penalty; 0 = OLS (no penalty)
  intercept = FALSE,
  include_fitted = TRUE,
  include_residuals = TRUE,
  cache = TRUE,
  method = c("normal", "qr")
)
```

**Key implementation decisions:**

1. **MLX Cholesky bridge strategy:** MLX exposes `mlx_linalg_cholesky(a, upper)` and `mlx_linalg_cholesky_inv` (check exact symbol in `mlx/c/linalg.h`). If a full `cho_solve` is not exposed in MLX C API, implement as two triangular solves using the existing `amatrix_mlx_solve_triangular_bridge` (R/backend.R:464). The new `amatrix_mlx_chol_solve_bridge` C function signature:
   ```c
   SEXP amatrix_mlx_chol_solve_bridge(SEXP A, SEXP B, SEXP upper);
   ```
   Flow: Cholesky factor A -> two triangular solves against B -> return f64 R matrix. f32 internal precision — the `precision = "strict"` check rejects this path.

2. **AF Cholesky bridge strategy:** ArrayFire exposes `af_cholesky` and `af_solve(A, B, AF_MAT_UPPER | AF_MAT_LOWER)`. Cleaner than MLX because the solve hop is native. Signature:
   ```c
   SEXP am_af_chol_solve_bridge(SEXP A, SEXP B, SEXP upper);
   ```
   Use existing row-major f32 conversion helpers (`copy_r_to_row_major_f32`, `copy_row_major_f32_to_r` at `arrayfire_bridge.c:10,18`).

3. **Resident path:** Add `chol_solve_resident(a_key, b_key, out_key)` to both backends. When `A` and `B` are both `adgeMatrix` with resident keys, this avoids any host transfer — a major win for the `(X'X + lambda*I)^-1 X'Y` path where `X'X + lambda*I` is p×p and `X'Y` is p×q with q in the millions. Both backends register the op via `capabilities = function() c(..., "chol_solve")`.

4. **SPD detection:** Both bridges check the Cholesky return status; on failure, `Rf_error("'a' is not positive definite")` (matches `base::chol` wording exactly so downstream `tryCatch` handlers keep working).

5. **`am_lm_fit` routing logic:**
   ```r
   # Pseudo-code for the new GPU path inside am_lm_fit
   if (inherits(X, "adgeMatrix") && X@precision == "fast") {
     XtX <- am_crossprod(X)                          # resident p x p
     XtY <- am_crossprod(X, Y)                       # resident p x q
     A <- if (lambda > 0) am_ewise("+", XtX, lambda * diag(p)) else XtX
     beta <- am_chol_solve(A, XtY)                   # resident p x q
   } else {
     # existing CPU / host path unchanged
   }
   ```
   The existing `.amatrix_lm_core` (line 314) stays intact for non-adge inputs and `precision = "strict"`.

6. **Cache integration:** Existing `.amatrix_lm_cache_value` (line 159) caches `XtX`. For `lambda != 0`, the penalized `XtX + lambda*I` is a different matrix, so cache the base `XtX` and apply `lambda*I` on each call. The existing `am_ridge_fit` already does this (line 651-652). No new cache machinery needed — re-use `.amatrix_penalty_matrix`.

7. **Dispatch priority for solve/chol:** Per user request, AF first, then MLX, then CPU. This matches the existing `.amatrix_backend_for` resolution order when both backends are registered — no code change needed, but verify with `amatrix_backend_plan(x, "chol_solve")` in a test.

8. **Why Cholesky, not QR:** X'X + lambda*I is SPD by construction when lambda >= 0. Cholesky is ~2x faster than QR and uses half the memory. QR on X directly is more numerically stable when X is rank-deficient, but `am_lm_fit(..., method = "qr")` already exists for that case. This plan only touches the `method = "normal"` path.

9. **Float precision loss analysis:** On f32, condition number kappa(X'X) = kappa(X)^2. For ill-conditioned X (kappa > 1e3), f32 Cholesky loses accuracy. Document in `?am_chol_solve` that `precision = "strict"` is recommended when kappa(X) > 1e3, and test with both well-conditioned and moderately-conditioned inputs.

**Acceptance criteria:**
- `am_chol_solve(A, B)` for small SPD A (8x8) agrees with `base::solve(A, B)` to `1e-5` (f32 backends) / `1e-10` (CPU).
- `am_lm_fit(X, Y, lambda=0)` returns identical `$coefficients` (within f32 tol) to `lm.fit(X, Y)$coefficients` for well-conditioned random X with shape 1000x50, Y shape 1000x100.
- `am_lm_fit(X, Y, lambda=0.1)` returns identical `$coefficients` to `am_ridge_fit(X, Y, lambda=0.1)`.
- When X is `adgeMatrix(..., precision="fast", preferred_backend="mlx")` or `"arrayfire"`, the solve does **not** materialize to host — verified by checking that a resident key is bound to the output `coefficients`.
- Non-SPD input to `am_chol_solve` errors with `"not positive definite"` message on all three backends.
- Conformance harness from Feature A includes `am_chol_solve` and passes on MLX and AF.
- `amatrix_backend_plan(X, "chol_solve")` returns AF first, MLX second, CPU third when all are registered.
- Build succeeds for both backend packages (`R CMD INSTALL backends/amatrix.mlx` and `.../amatrix.arrayfire`) on an Apple Silicon machine.

---

### Feature C: Hardened S3/S4 Dispatch + ALTREP Research Note

**Files to create:**
- `R/dispatch-hardening.R` — new R file with additional S4 methods and `.onLoad` hooks.
- `tests/testthat/test-dispatch-primitives.R` — regression tests for `base::%*%` pass-through.
- `docs/altrep-research.md` — research note (not user docs) on true ALTREP path.

**Files to modify:**
- `R/methods-dense.R` — add missing LHS/RHS combinations; ensure `Matrix %*% adgeMatrix` (dense Matrix) routes through GPU.
- `R/zzz.R` — add `.onLoad` hook that registers S3 methods for `%*%`, `crossprod`, `tcrossprod` via `.S3method()` for `adgeMatrix` class (complements the S4 methods for packages that use S3 dispatch).
- `NAMESPACE` — add `S3method("%*%", adgeMatrix)` etc. if the `.S3method` approach is chosen.

**Key implementation decisions:**

1. **Scope:** This feature is NOT true ALTREP. It is a hardening pass over existing S3/S4 dispatch to plug remaining holes where `base::%*%` and `base::crossprod` fall through to `as.matrix()` materialization. True ALTREP is deferred to a future major release (research note only).

2. **Gap audit first:** Before writing code, run a test matrix of `(LHS_class, RHS_class, op)` combinations — `adgeMatrix`, `matrix`, `dgeMatrix`, `dgCMatrix`, `numeric`, `integer`, `array` — against `%*%`, `crossprod`, `tcrossprod`. For each combination, assert that the result is still `adgeMatrix` (GPU-resident). Expected gaps based on the current code:
   - `array %*% adgeMatrix` — likely falls through (S4 doesn't know `array` class).
   - `adgeMatrix %*% numeric` — routed via `am_matmul` which promotes to column matrix (works, but verify).
   - `crossprod(adgeMatrix, matrix)` — verify it reaches `am_crossprod` and not the base method.
   - `base::%*%(matrix, adgeMatrix)` when called as a primitive (some packages use `.Primitive("%*%")`) — definitely falls through.

3. **`.Primitive("%*%")` trap:** R's `%*%` is a primitive that consults S4 dispatch only if either argument has S4 class. `adgeMatrix` does, so normally this works. The failure mode is packages that do `m <- as.matrix(x); m %*% v` — we can't fix that without ALTREP. Document it.

4. **S3 method registration:** R's `base::%*%` primitive uses `UseMethod` for S3 classes. Register `.S3method("%*%", "adgeMatrix", am_matmul)` in `.onLoad`. This catches cases where a caller does `m %*% v` with `m` having class attribute but S4 dispatch not firing (rare but happens).

5. **New S4 methods to add:**
   ```r
   setMethod("%*%", signature(x = "adgeMatrix", y = "array"), function(x, y) am_matmul(x, y))
   setMethod("%*%", signature(x = "array", y = "adgeMatrix"), function(x, y) am_transpose(am_crossprod(y, t(x))))
   setMethod("crossprod", signature(x = "adgeMatrix", y = "matrix"), function(x, y, ...) am_crossprod(x, y, ...))
   setMethod("crossprod", signature(x = "matrix", y = "adgeMatrix"), function(x, y, ...) {
     # t(x) %*% y where y is adgeMatrix: t(am_matmul(y_t, x_t)) = t(tcrossprod(y, x))
     am_transpose(am_tcrossprod(y, x))  # verify dims
   })
   setMethod("tcrossprod", signature(x = "adgeMatrix", y = "matrix"), function(x, y, ...) am_tcrossprod(x, y, ...))
   ```
   Each new method needs a round-trip test in `test-dispatch-primitives.R`.

6. **ALTREP research note content (`docs/altrep-research.md`):**
   - What ALTREP is: R 3.5+ mechanism for custom SEXP classes with lazy `DATAPTR`.
   - Why we'd want it: `x %*% v` where `x` is a base `matrix` stored GPU-side goes to `do_matprod` in `src/main/array.c`, which calls `REAL(x)` which calls `DATAPTR(x)`. ALTREP lets us intercept that (at the cost of materialization on any `DATAPTR` access).
   - Why it's hard: (a) `DATAPTR` is called by innumerable base R paths (subsetting, copying, `.Call` into other packages), so any call silently materializes and kills residency; (b) ALTREP classes need a companion of method hooks (`Elt`, `Get_region`, `Is_sorted`, etc.) that all need to handle the GPU-resident case; (c) we can't intercept `%*%` specifically — we'd need to intercept `DATAPTR` and rely on the C-level matmul noticing our class and routing to our kernel, which requires either a custom `%*%` primitive replacement (not possible in user code) or a `Matrix`-style SEXP wrapper that base R doesn't know about.
   - Best path forward: track R-core's ALTREP matrix proposal (pre-PEP discussions in r-devel); if it lands, revisit. Otherwise, continue S4 dispatch + targeted S3 method registration.
   - Estimated effort if pursued: 400-800 LoC of C, 2-4 weeks, high risk of silent residency loss bugs.

**Acceptance criteria:**
- Dispatch gap audit table committed to `tests/testthat/test-dispatch-primitives.R` with one test per `(LHS, RHS, op)` combination.
- For every `(LHS, RHS, op)` row where a GPU backend is available, the result is `adgeMatrix` or has a resident key bound.
- `irlba::irlba(x_adge, nv=5)` runs without falling back to `as.matrix()` on the hot path (check by registering a recording backend and asserting `matmul` count equals the number of Lanczos iterations).
- `docs/altrep-research.md` exists and covers: what, why, blockers, effort estimate, decision ("defer to future major release").
- No new C code in `amatrix` core package for Feature C.

---

## Success Criteria

- Conformance harness passes on CPU-only, CPU+MLX, CPU+AF, and CPU+MLX+AF machines.
- `am_lm_fit(X, Y, lambda=0)` on 1000x50 design with 1000x100,000 response returns coefficients agreeing with `lm.fit` to f32 tolerance, in under 500ms on Apple Silicon with MLX (benchmark target, not hard gate).
- `am_chol_solve` is exported and documented with roxygen.
- All existing tests (`devtools::test()` on amatrix core) continue to pass.
- Both backend packages (`R CMD INSTALL backends/amatrix.mlx`, `.../amatrix.arrayfire`) build and test cleanly.
- S3/S4 dispatch gap audit shows zero unexpected host materialization for the supported combinations.

---

## Dependencies Between Features

```
A (Conformance Harness) -- independent, must ship first
    |
    +-- provides test coverage for -->
    |
    v
B (GPU Cholesky + am_lm_fit)  -- depends on A for regression guarantees
    |
    +-- new chol_solve op -->  added to conformance matrix in A
    |
    v
C (Hardened Dispatch)  -- depends on A for dispatch-gap regression tests
    |
    +-- exercises am_chol_solve via Matrix/base generics, so benefits from B being in place
```

**Hard dependencies:** B cannot land without A's conformance harness (otherwise the new C bridges are untested). C cannot land without A's dispatch-coverage harness.

**Soft dependencies:** C's dispatch hardening exercises `am_chol_solve` through `solve(crossprod(adge))` style calls, which is a nice validation path for B but not strictly required.

---

## Estimated Complexity

| Feature | LoC (approx) | Files touched | C code | Risk | Timeline |
|---|---|---|---|---|---|
| A | 400-600 R | 4 new + 2 modified | none | low | 2-3 days |
| B | 300-500 C + 200 R | 5 new + 7 modified | yes (MLX+AF) | medium | 4-6 days |
| C | 150-250 R + docs | 3 new + 3 modified | none | low-medium | 2-3 days |

**Total:** ~8-12 days of focused work, assuming backend packages build cleanly on the dev machine.

---

## Open Questions

Items to revisit during or after execution (will be written to `.omc/plans/open-questions.md`):

- [ ] MLX C API: confirm whether `mlx_linalg_cholesky` exists or whether we need to compose triangular solves manually.
- [ ] Should `am_chol_solve(A, B)` accept `B` as a numeric vector (length p) and auto-promote to a p x 1 matrix? (likely yes, matching `am_solve` convention at wrappers.R:258).
- [ ] Does `capabilities()` string list need a formal enum in `backend-registry.R`, or is the current free-form character vector fine?
- [ ] For Feature C gap audit, should we add a `dispatch_trace` debug mode to `amatrix_dispatch_op` that logs every op + chosen backend, to make "unexpected host materialization" visible at runtime?
- [ ] ALTREP research note: should we track the r-devel matrix-ALTREP discussion in a beads issue so we get notified if it moves?
