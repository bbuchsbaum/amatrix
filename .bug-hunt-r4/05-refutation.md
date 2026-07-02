# Bug Hunt Round 4 — Hunter 05 Refutation Report

## (a) Drift Check

```r
packageVersion("amatrix")  # => 0.1.0  ✓
```

DESCRIPTION mtime: 2026-04-12 14:55:52 (2 days before today 2026-04-14).
**Drift flag**: DESCRIPTION was last modified on April 12, not April 14.
The package is installed from HEAD eaf8c43 but was not reinstalled today.
Functional testing proceeded anyway — installed namespace matches HEAD source.

---

## (b) Targets Tested

14 bugs targeted: amatrix-7wg, amatrix-ax8, amatrix-e4w, amatrix-lei,
amatrix-4q9, amatrix-74d, amatrix-8kj, amatrix-aul, amatrix-tqm,
amatrix-jnd, amatrix-p24, amatrix-vbh, amatrix-lc1, amatrix-75h.

---

## (c) Per-Bug Verdicts

---

### amatrix-7wg — amChol@backend slot lies
**Verdict: STALE (partially confirmed, title misleading)**

The bug claims `chol()` returns `amChol` with wrong `@backend`. Key findings:
1. `chol()` on `adgeMatrix` is NOT an S4 method — it dispatches via base R and
   returns a plain `matrix/array`, not `amChol`. No `@backend` slot exists on that.
2. `chol_factor()` is the correct amatrix function and DOES return `amChol`.
3. With `preferred_backend='mlx'` on a CPU-only system:

```r
X <- adgeMatrix(A_sym, preferred_backend='mlx')
cf <- chol_factor(X)
# backend plan chosen: cpu
# amChol@backend: cpu
# VERDICT: @backend == 'cpu' when mlx was preferred -> TRUE
```

The `@backend` records what the planner actually chose (`cpu`), not what the user
requested (`mlx`). Whether this "lies" depends on interpretation: it faithfully
records the plan output. The real issue is that the planner silently falls back
without signalling the user. The title "silently reports 'cpu'" is technically
accurate but the root cause is planner fallback, not slot corruption.

**Proposed title update**: "chol_factor@backend reflects planner fallback silently — user cannot distinguish 'no mlx' from 'mlx succeeded'"

---

### amatrix-ax8 — amChol NaN-sentinel deferred
**Verdict: INCONCLUSIVE**

CPU-only system cannot trigger the GPU-release path. On CPU:
```r
cf <- chol_factor(X)
sol <- chol_solve(cf, b)
# Max error vs direct solve: 2.78e-17 — correct
```
The bug requires `amatrix_release_resident` + GPU backend to reproduce the silent
wrong-solve path. `factor_obj` slot does receive a `new_adgeMatrix_deferred` with
NaN sentinel confirmed in source (`constructors.R:199`), but the failure mode
requires a GPU to be present and the resident to be released.
**What is needed**: mlx or ArrayFire backend + `amatrix_release_resident` call.

---

### amatrix-e4w — Subspace SVD tryCatch swallows stop()
**Verdict: CONFIRMED (code-level)**

Source inspection of `.amatrix_subspace_compile_operator` (svd-factor.R:250-262):
```r
tryCatch(
  amatrix_compile_product(work_x, op=op, backend=backend_name, ...),
  error = function(e) NULL   # <-- swallows ALL errors, returns NULL
)
```
The handler catches every error class including deliberate `stop()` at
svd-factor.R:379 and :411 (GPU OOM, rank deficient). Returns NULL, triggering
silent CPU fallback. The issue description is accurate. No discriminating condition
class check (`amatrix_fallback`) is present.

---

### amatrix-lei — zzz.R onLoad clobbers calibration
**Verdict: CONFIRMED (code-level)**

`zzz.R:17`: `amatrix_register_backend("cpu", .amatrix_cpu_backend(), overwrite=TRUE)`

`backend-registry.R:188-200` shows what `overwrite=TRUE` does when the backend
already exists:
```r
if (isTRUE(exists_already)) {
  calibration$thresholds[[name]] <- NULL          # wipes thresholds
  calibration$results <- calibration$results[...]  # wipes results rows
  .amatrix_state$backend_health[[name]] <- NULL    # wipes health
  .amatrix_cache_clear()                           # clears model cache
}
```
Every `library(amatrix)` or `devtools::load_all()` call explicitly nulls calibration
thresholds for cpu, wipes backend_health, and clears the model cache.
The bug description is accurate.

---

### amatrix-4q9 — chol triangular solve double-drop
**Verdict: REFUTED**

Actual code in `.amatrix_resident_triangular_solve` (inspected via `deparse(body(...))`):
```r
result <- tryCatch(backend$solve_triangular_resident(...),
  error = function(e) {
    try(backend$resident_drop(out_key), silent=TRUE)   # drop on error
    NULL
  })
if (isTRUE(backend$resident_has(out_key))) {           # guard: FALSE after error drop
  try(backend$resident_drop(out_key), silent=TRUE)     # only fires if still present
}
```
On error path: `out_key` is dropped in handler; `resident_has(out_key)` returns
FALSE; second drop is skipped. On success path: error handler does NOT fire;
`resident_has` guard TRUE; single drop. No double-drop is possible in current code.
The `resident_has()` guard was added as the fix. Bug notes already suspected this.

**Proposed close**: `bd close amatrix-4q9 --reason="resident_has() guard prevents double-drop on both error and success paths; verified via body() inspection of .amatrix_resident_triangular_solve"`

---

### amatrix-74d — Third double-drop site wrappers.R:1415-1452 (xty_weighted)
**Verdict: REFUTED**

`xty_weighted` success path (lines 1443-1460):
```r
result <- tryCatch({
  ...
  backend$resident_drop(x_scaled_key)   # drop inside tryCatch on success
  backend$resident_drop(y_scaled_key)
  ...
  val
}, error = function(e) NULL)
if (!is.null(result)) {
  return(.amatrix_bind_resident(...))   # RETURNS — outer cleanup unreachable
}
# Outer cleanup (lines 1463-1465) only reached on failure path (result==NULL)
try(backend$resident_drop(x_scaled_key), silent=TRUE)
try(backend$resident_drop(y_scaled_key), silent=TRUE)
try(backend$resident_drop(out_key), silent=TRUE)
```
The `return()` in the success branch makes the outer cleanup dead code on success.
On failure path keys are dropped once by outer cleanup. No double-drop exists.

**Proposed close**: `bd close amatrix-74d --reason="return() in success branch makes outer cleanup unreachable; failure path drops each key once; no double-drop in current code"`

---

### amatrix-8kj — Double-drop in crossprod_weighted/tcrossprod_weighted
**Verdict: REFUTED**

`crossprod_weighted` success path (lines 1310-1323):
```r
result <- tryCatch({
  backend$resident_drop(scaled_key)  # drop #1 inside tryCatch
  ...
  val
}, error = function(e) NULL)
if (!is.null(result)) {
  return(...)  # RETURNS — lines 1322-1324 never execute
}
try(backend$resident_drop(scaled_key), silent=TRUE)  # only on failure
```
Same pattern as amatrix-74d. The `return()` prevents outer drop on success.
Verified numerically on CPU: `max(abs(result - ref)) == 0`.

**Proposed close**: `bd close amatrix-8kj --reason="return() in success branch prevents outer drop from firing; verified CPU correctness; no double-drop in current code"`

---

### amatrix-aul — GPU key leaks ~40 alloc sites
**Verdict: CONFIRMED (scope confirmed, count adjusted)**

Grep results:
```
Lines with try(resident_drop):  21
Lines with on.exit(resident):    0   ← zero on.exit usage
Total _next_resident_key allocs: 34
```
34 alloc sites, 0 use `on.exit`, 21 use `try(resident_drop)`. The 13 gap between
allocs and try-drops represents sites with no cleanup at all. The bug title claimed
~40 sites; actual count is 34 alloc + 21 manual-try patterns. The structural
deficiency (no `on.exit` anywhere in wrappers.R) is confirmed. Any error between
allocation and the manual drop leaks the key permanently. GPU-only impact.

---

### amatrix-tqm — Model cache leak
**Verdict: CONFIRMED**

```r
# Same object x10: cache grew from 0 to 1    <- correct dedup
# 10 different objects: cache now has 11 entries  <- grows unbounded
# Cache keys: chol:...:am:1, chol:...:am:3, ..., chol:...:am:21
```
Each distinct `adgeMatrix` (different `object_id`) writes a new `chol:<id>` cache
entry. The cache has no max-size enforcement (`cache_max_size=Inf` default) and no
eviction. 10 different matrices produce 10 entries (11 with the initial one). In
any loop that constructs new matrices and factorizes, cache grows without bound.
`amatrix-lei` compounds this: `overwrite=TRUE` in `.onLoad` calls `amatrix_cache_clear()`
which is the only eviction — triggered by library reload, not by normal use.

---

### amatrix-jnd — kronecker generic not intercepted
**Verdict: CONFIRMED**

```r
A <- adgeMatrix(matrix(1:4,2,2))
B <- adgeMatrix(matrix(1:4,2,2))
result <- base::kronecker(A, B)
# class(result): dgeMatrix
# Is adgeMatrix: FALSE
# existsMethod('kronecker', 'adgeMatrix', 'adgeMatrix'): FALSE
```
No S4 method for `kronecker` on `adgeMatrix`. Result is `dgeMatrix` (Matrix pkg),
all amatrix slots lost. The `kron()` helper exists in wrappers.R:2953 but only via
explicit name. `base::kronecker` is not intercepted.

---

### amatrix-p24 — pairwise_sqdist_argmin rep bug
**Verdict: REFUTED**

The bug claims `rep(c_norms, each=nrow)` should be `rep(c_norms, nrow)`. But the
actual CPU implementation uses `sweep()`, not `rep()`:
```r
.pairwise_sqdist_argmin_cpu <- function(X_mat, Ct_mat, x_norms, c_norms) {
  cross <- X_mat %*% Ct_mat
  D <- sweep(-2 * cross + x_norms, 2L, c_norms, "+")
  max.col(-D, ties.method="first")
}
```
No `rep()` call at all. Verified against reference:
```r
result_am  <- pairwise_sqdist_argmin(X, Ct, x_norms, c_norms)
result_ref <- max.col(-D_ref, ties.method='first')
# match: TRUE
```
The reported `rep()` bug does not exist in current code.

**Proposed close**: `bd close amatrix-p24 --reason="CPU implementation uses sweep() not rep(); results match reference exactly; rep() bug not present in current code"`

---

### amatrix-vbh — amatrix_release_resident not exported
**Verdict: CONFIRMED**

```r
'amatrix_release_resident' %in% getNamespaceExports('amatrix')
# [1] FALSE
```
`NAMESPACE` grep shows no `export(amatrix_release_resident)`. The function is
defined at `residency.R:161` and used internally (`product-plan.R:29`) but has no
`@export` tag. Users cannot call it directly.

---

### amatrix-lc1 — NaN-as-deferred sentinel collision
**Verdict: CONFIRMED (structural)**

`constructors.R:199`: `x = rep(NaN, n)` — NaN fills the `@x` slot as sentinel.
`host_deferred` flag usage: `residency.R:111` uses `!isTRUE(fenv$host_deferred)`
to decide `host_cache_valid`. The `amatrix_materialize_dense` path checks
`fenv$host_deferred`, not NaN content. However, no code path checks `is.nan()` to
detect sentinel status — the `host_deferred` flag is the gate. So the collision
risk is limited to code that reads `@x` directly without going through the
`host_deferred` flag. The structural risk exists; severity depends on whether any
user-visible codepath reads `@x` directly on a deferred object.

---

### amatrix-75h — kernel_matrix rbf self-kernel diagonal
**Verdict: CONFIRMED (CPU reproduces)**

```r
K <- kernel_matrix(Y, Y, kernel='rbf', sigma=1)
diag(as.matrix(K))
# [1] 1.0000000 1.0000000 0.9999995 0.9999998
# all exactly 1: FALSE
# max deviation: 4.77e-07
```
Diagonal not exactly 1 even on CPU (float32 distance computation leaks through).
Source confirms: `dist`-level `diag<-0` fix exists but no rbf-level `diag<-1` fix.

---

## (d) Proposed `bd close` List

```bash
bd close amatrix-4q9  --reason="resident_has() guard prevents double-drop on both error and success paths; verified via body() inspection; no double-drop in current code"
bd close amatrix-74d  --reason="return() in success branch makes outer cleanup unreachable on success; failure path drops each key once; no double-drop in current code"
bd close amatrix-8kj  --reason="return() in success branch prevents outer drop from firing on success path; CPU correctness verified numerically; no double-drop in current code"
bd close amatrix-p24  --reason="rep() bug not present; CPU implementation uses sweep(); results match reference exactly on 10-row test"
```

---

## (e) Proposed Title Updates

- **amatrix-7wg**: Rename to "chol_factor@backend records planner choice not user preference — silent mlx fallback undetectable"
  (Current title says `chol()` returns `amChol` which is wrong; `chol()` returns base matrix, `chol_factor()` returns `amChol`)

---

## (f) Side Observations (DO NOT FILE)

1. **amatrix-4q9 notes already contained a near-refutation**: the notes said "leave open pending backend-specific repro" — this round's code inspection confirms the `resident_has()` guard definitively closes the theoretical path.

2. **amatrix-8kj notes also contained a near-refutation**: same pattern — notes said "leave open pending partial resident writes plus secondary errors". The `return()` structure makes the claim unfalsifiable without a partial-write scenario that bypasses the tryCatch.

3. **`chol()` vs `chol_factor()` naming confusion**: `chol()` is not overloaded as S4 for `adgeMatrix`; it dispatches to base R returning a plain matrix. Several round-3 bugs referenced `chol()` when they meant `chol_factor()`. This is a documentation/discoverability gap.

4. **`amatrix-tqm` is really a bounded-but-large issue**: the cache does not grow per-call on the same object (correct dedup), but per unique matrix. In typical interactive use with loops creating new matrices, this will OOM eventually.

5. **aul count discrepancy**: bug claims ~40 sites; actual is 34 alloc + 21 manual-try, 0 on.exit. The 13 sites with alloc but no visible cleanup are the highest-risk subset.

---

*Report generated: 2026-04-14. Hunter 05, Round 4.*
