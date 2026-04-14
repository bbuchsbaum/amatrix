# Round 3 Bug Hunt: GPU Cross-Backend Conformance
# Hunter 02 — Execution-Based Verification
# Date: 2026-04-14

## Environment
- R 4.5.1 / macOS 23.3.0
- Backends loaded: amatrix (CPU), amatrix.mlx, amatrix.arrayfire
- Method: live Rscript execution, no source modification
- Conformance suite: 94 PASS / 0 FAIL / 2 SKIP (opencl not installed)

---

## (a) Verification of Round-2 Inferred Bugs

### amatrix-dev — host_cache_valid unconditionally TRUE (03-resident-coherence.md V1)
**Status: REFUTED (for the normal GPU matmul path)**

Probe: created adgeMatrix on mlx, ran `A %*% A`, inspected `host_cache_valid` on result.

```
Before GPU push - host_cache_valid: FALSE
After GPU matmul - host_cache_valid on result: FALSE
After GPU matmul - host_deferred on result: FALSE
Max diff from expected (m %*% m): 0
RESULT: PASS - materialize returns correct value
```

The current implementation does NOT set `host_cache_valid = TRUE` on matmul results. The R2 inferred hazard (bind-then-serve-stale-@x) was not triggered by this path. The stale-read scenario from residency.R:97 requires a code path that calls `.amatrix_bind_resident` on an object whose `@x` slot predates the GPU op — this requires internal API access not reachable from normal user-facing ops in the current build. **The normal GPU result path correctly leaves `host_cache_valid = FALSE`.**

---

### amatrix-qm2 — mixed-backend registry corruption (04-composition-chains.md F1)
**Status: INCONCLUSIVE (lazy promotion, no pre-push)**

Probe: `A_mlx %*% B_af` where A was preferred_backend="mlx", B was preferred_backend="arrayfire".

```
Mixed matmul succeeded, result backend: mlx
B has no residency entry (may never have been pushed)
Numerical correctness max diff: 0
```

B had no residency entry at the time of the op — neither object was pre-pushed to GPU before the matmul. The op used lazy promotion, so there was no ArrayFire buffer to corrupt. The F1 corruption scenario requires B to already be GPU-resident on AF before the mixed op. This is not triggered by the normal constructor path. The structural risk in `.amatrix_bind_resident` overwriting without prior `_release` (residency.R:84) remains a code-level concern but **could not be reproduced with normal API usage**. Verdict: inferred risk not refuted, but not executable via public API in current build.

---

### amatrix-cth — deferred adgeMatrix dangling after release (04-composition-chains.md G2)
**Status: INCONCLUSIVE — `amatrix_release_resident` not exported**

The public API function `amatrix_release_resident` does not exist. The internal function `.amatrix_release_resident` exists but is not exported. `new_adgeMatrix_deferred` exists internally.

When the release was performed via `.amatrix_release_resident(B)` where B had a valid `@x` host copy (length 9, non-NaN), the post-release op succeeded correctly (diff = 0). This is the safe path — B had a real host copy so fallback to `@x` worked.

The dangerous G2 scenario (deferred object with NaN `@x` + release + plan fallback) requires `new_adgeMatrix_deferred` to produce an object with sentinel `@x`, which is not reachable without internal API access. **Inconclusive due to missing public API surface for the dangerous path.**

New finding: `amatrix_release_resident` is documented in round-2 sketches but is NOT exported from the package. User-facing GPU memory management is therefore not possible without `:::`.

---

### amatrix-fqh — cache keys miss backend identity (06-cache-invalidation.md B3/B4)
**Status: CONFIRMED (execution evidence)**

Probe: computed `chol_factor` for `adgeMatrix(m, preferred_backend="cpu")`. Inspected cache keys.

```
Cache keys after CPU chol: chol:20260414081624...:am:1
Backend name in cache key: FALSE
CONFIRMED BUG: cache keys have no backend component
```

Additional finding: when `chol_factor` is called on an `adgeMatrix` with `preferred_backend="mlx"`, the resulting `amChol` object has `@backend = "cpu"` — the actual dispatch fell back to CPU. The cache key `chol:<object_id>` records no backend. So:

1. If a different backend later computes chol for the same `object_id` (impossible in practice since object_id is unique per object), the stale factor would be returned.
2. More concretely: `amChol@backend = "cpu"` even when computed on an mlx-preferred object — the `@backend` slot silently misreports the computation backend. This makes the cached factor's provenance unverifiable.

**Min repro:**
```r
library(amatrix); library(amatrix.mlx)
m <- crossprod(matrix(rnorm(16),4,4)) + 5*diag(4)
A <- adgeMatrix(m, preferred_backend = "mlx")
f <- chol_factor(A)
stopifnot(f@backend == "mlx")  # FAILS: f@backend == "cpu"
# Cache key: chol:<object_id>  -- no backend component
```

---

### amatrix-dmy — LOO CV mixed-precision (04-composition-chains.md E1)
**Status: REFUTED for current test case (may be input-dependent)**

Probe: `lm_loo_cv(X, y)` with both `precision="strict"` and `precision="fast"` on a 20x4 design matrix.

```
LOO CV strict (float64) residuals == reference: max diff = 0
LOO CV fast  (float32) residuals == reference: max diff = 0
```

Both precision modes returned identical results with zero diff against pure-R reference. This suggests either: (a) MLX is computing at full float64 precision regardless of `precision="fast"` setting for this size, (b) the mixed-precision path in `qr-downdate.R:143,145` is not triggered for this input, or (c) the bug was fixed between R2 and R3. **Cannot confirm the mixed-precision correctness bug for the tested input size and seed.**

---

## (b) New Chain-Bug Findings

### NEW-01 — Cache leak confirmed under mutation loop [CONFIRMED]
**Severity: MEDIUM | File: R/chol-factor.R:154**

10 sequential mutations (`A[1,1] <- A[1,1] + 0.001`) each producing a new `object_id` leave 10 orphaned `chol:<old_id>` entries in `.amatrix_state$model_cache`. Entries are never evicted with default `cache_max_size = Inf`.

```
Cache entries before mutation loop: 1
Cache entries after 10 mutations: 11
CONFIRMED BUG (V3): cache leaking - 10 orphaned entries
```

**Min repro:**
```r
library(amatrix)
m <- crossprod(matrix(rnorm(16),4,4)) + 5*diag(4)
A <- adgeMatrix(m)
chol_factor(A)
for (i in 1:10) { A[1,1] <- A[1,1] + 0.001; chol_factor(A) }
length(ls(amatrix:::.amatrix_state$model_cache))  # returns 11, not 1
```

In a LOO-CV loop with n=1000 this accumulates 1000 leaked chol entries.

---

### NEW-02 — amChol@backend slot silently misreports CPU dispatch for MLX-preferred objects [CONFIRMED]
**Severity: MEDIUM | File: R/chol-factor.R**

When `chol_factor` is called on an `adgeMatrix` with `preferred_backend="mlx"`, the returned `amChol` has `@backend = "cpu"`. This means:
- The `@backend` slot cannot be trusted to identify which backend computed the factor.
- The existing cache-key design (no backend component) compounds this: a downstream check comparing `f@backend` against the current preferred backend will always see "cpu" and cannot detect a backend mismatch.
- Any future logic that gates cache reuse on `f@backend == preferred_backend` will silently bypass the gate for CPU-fallback factors on GPU-preferred objects.

```
f_mlx@backend: cpu   (expected: mlx)
f_cpu@backend: cpu
```

This is distinct from amatrix-fqh (which concerns the cache key string); this is a slot-level metadata integrity bug on the returned factor object.

---

### NEW-03 — `amatrix_release_resident` not exported; GPU memory management not user-accessible [CONFIRMED]
**Severity: LOW (API gap) | File: NAMESPACE / R/residency.R**

The function `amatrix_release_resident` (called in round-2 sketches and the G2 scenario) does not exist as an exported symbol. The internal `.amatrix_release_resident` exists. Users cannot explicitly free GPU buffers from the public API — they rely entirely on GC finalizers. This:
- Makes it impossible to reproduce or avoid the G2 (deferred + plan + release) scenario from user code.
- Means GPU memory pressure cannot be managed proactively in long-running sessions.
- All round-2 sketches that call `amatrix_release_resident` are non-executable as written.

---

### NEW-04 — Chain 3: chol_factor result accessor mismatch (`$` vs `@`) [CONFIRMED]
**Severity: LOW (usability/doc) | File: R/chol-factor.R**

`chol_factor` returns an `amChol` S4 object. Round-2 sketches and user-facing documentation use `$L` (list-style accessor) but the actual slot is `@factor` (plain R matrix, not adgeMatrix). `amatrix_materialize_dense(f$L)` fails; `f@factor` works and is a plain matrix.

```
chol_factor slotNames: factor, factor_obj, source_id, precision, backend
chol_factor names: (empty - not a list)
```

The `$` operator is not defined for `amChol`, causing a hard error in any code that follows the round-2 sketch pattern.

---

### NEW-05 — H1 sparse->dense coerce: backend IS preserved (round-2 refuted)
**Status: REFUTED**

```
adgCMatrix backend: mlx
Coerced to adgeMatrix, backend: mlx
Backend preserved: mlx
```

`as(adgCMatrix_mlx, "adgeMatrix")` correctly preserves `preferred_backend = "mlx"`. The H1 inferred bug from round-2 does not reproduce. A `setAs("adgCMatrix","adgeMatrix")` path apparently exists and transfers metadata correctly.

---

### NEW-06 — Post-GPU-release fallback to host `@x` works (G2 safe path confirmed)
**Status: PASS (safe path)**

After `.amatrix_release_resident(B)` where B had a populated `@x`, subsequent `B %*% A` correctly fell back to the host `@x` data with zero numerical error. The dangerous G2 path (NaN `@x` + release) requires the internal deferred constructor and is not exercisable from the public API.

---

## (c) Test Suite Gaps Under Real GPU

### Conformance suite results
```
FAIL 0 | WARN 0 | SKIP 2 | PASS 94
Duration: 1.1s
Skipped: amatrix.opencl backend not installed (lines 444, 464)
```

All 94 assertions pass. The suite is fast (1.1s) suggesting most tests use small matrices and may not stress GPU code paths.

### Identified gaps

**Gap 1 — No test for cache key backend identity.**
No test verifies that `chol_factor` or `am_svd` results computed on different backends for the same numerical data produce cache entries distinguishable by backend. The `amChol@backend` misreport (NEW-02) and the backend-free cache key (amatrix-fqh) are both untested.

**Gap 2 — No test for cache growth under repeated mutation.**
No test checks that `length(ls(.amatrix_state$model_cache))` stays bounded after N mutations. The V3 cache leak (NEW-01) is therefore untested.

**Gap 3 — No test for `amChol@backend` slot correctness.**
No assertion verifies `chol_factor(mlx_obj)@backend == "mlx"`. The slot silently reports "cpu" for MLX-preferred objects (NEW-02).

**Gap 4 — No mixed-backend op tests with pre-pushed operands.**
All conformance tests appear to use objects without explicit prior GPU pushes. The F1 mixed-backend registry corruption scenario (amatrix-qm2) is not covered because it requires both operands to be resident on different backends before the op.

**Gap 5 — No test for `amatrix_release_resident` (or its internal equivalent) + subsequent op.**
The G2 deferred+release+plan scenario is untested. The internal `.amatrix_release_resident` fallback-to-host path is untested in the suite.

**Gap 6 — LOO CV precision test.**
No test covers `lm_loo_cv` with `precision="fast"` vs `precision="strict"` against a pure-R reference, which is the only way to detect the E1 mixed-precision bug if it exists for certain input shapes.

---

## Summary Table

| Bug ID | R2 Label | Status | Severity | Evidence |
|--------|----------|--------|----------|---------|
| amatrix-dev | host_cache_valid unconditional TRUE | REFUTED (normal path) | — | host_cache_valid=FALSE on all GPU results tested |
| amatrix-qm2 | mixed-backend registry corruption | INCONCLUSIVE | HIGH | No pre-push in normal API; structural risk unverified |
| amatrix-cth | deferred dangling after release | INCONCLUSIVE | HIGH | `amatrix_release_resident` not exported; safe path works |
| amatrix-fqh | cache key no backend identity | CONFIRMED | MEDIUM | Direct execution: `grep("mlx\|cpu", keys) == FALSE` |
| amatrix-dmy | LOO CV mixed precision | REFUTED (tested inputs) | — | max diff = 0 vs pure-R reference |
| NEW-01 | Cache leak under mutation loop | CONFIRMED | MEDIUM | 10 mutations → 11 orphaned cache entries |
| NEW-02 | amChol@backend slot misreports CPU for MLX obj | CONFIRMED | MEDIUM | f@backend="cpu" when preferred_backend="mlx" |
| NEW-03 | amatrix_release_resident not exported | CONFIRMED | LOW | `exists("amatrix_release_resident")` = FALSE |
| NEW-04 | amChol `$` accessor fails (S4 not list) | CONFIRMED | LOW | `$ operator not defined for this S4 class` |
| H1 | adgCMatrix→adgeMatrix loses backend | REFUTED | — | backend preserved in actual coerce |
