# Cache Invalidation Bug Hunt — Round 2

## 1. Cache Inventory

| Cache | Storage Location | Keys | Values | Readers | Writers | Clearers |
|---|---|---|---|---|---|---|
| **backend registry** | `.amatrix_state$backends` (env) | backend name (string) | backend list | `amatrix_backend_plan`, `amatrix_dispatch_op`, `amatrix_backend_for`, `amatrix_backend_status` | `amatrix_register_backend` | nothing (no unregister) |
| **calibration** | `.amatrix_state$calibration` (list) | `cal$thresholds[[backend]][[op]]` | numeric threshold | `.amatrix_calibration_ok` (called from `amatrix_backend_plan`) | `amatrix_calibrate` (session + disk), `.amatrix_load_calibration` (lazy disk load) | nothing on re-register |
| **disk calibration** | `~/.cache/R/amatrix/calibration.rds` | backend + op nested keys | same as above | `.amatrix_load_calibration` | `amatrix_calibrate(persist=TRUE)` | nothing |
| **backend health** | `.amatrix_state$backend_health` (env-as-list) | backend name | `list(status, reason, timestamp)` | `amatrix_backend_status` | `amatrix_backend_health_probe`, `amatrix_dispatch_op` (on error), `amatrix_backend_health_mark` | nothing on re-register |
| **fallback log** | `.amatrix_state$fallback_log` (list) | integer index | event record | `amatrix_fallback_log` | `amatrix_dispatch_op` on fallback | `amatrix_fallback_log_reset` only |
| **model cache** (LRU) | `.amatrix_state$model_cache` (env) + `cache_atime` (env) | string key e.g. `"chol:<object_id>"` | amChol / amQR / amSVD factor objects | `chol-factor.R`, `svd-factor.R`, `.amatrix_cache_get` | `.amatrix_cache_set` | `amatrix_gc(cache=TRUE)`, LRU eviction |
| **residency registry** | `.amatrix_state$residency` (env) | `"obj:<object_id>"` | `list(backend, resident_key, sparse, finalizer_env)` | `.amatrix_live_resident_backend`, `.amatrix_resident_entry`, `amatrix_dispatch_op` | `.amatrix_bind_resident` | finalizer, `.amatrix_release_resident`, `.amatrix_drop_resident_binding` |
| **am_product_plan closure** | closure environment of `plan_fun` in `product-plan.R` | N/A (captured variable `backend_name`) | string backend name + `lhs_bound` object | every `plan_fun(y)` call | `amatrix_compile_product` | no invalidation — plan is fire-and-forget |
| **balanced-deprecation flag** | `.amatrix_state$balanced_deprecation_warned` | — | logical | `.amatrix_balanced_deprecation_warned` | `.amatrix_warn_balanced_deprecation_once` | never |
| **default policy / precision** | `.amatrix_state$default_policy`, `$default_precision` | — | string | `amatrix_default_policy`, `amatrix_backend_preference` | `amatrix_set_default_policy/precision` | never (by design) |

**Total caches: 10**

---

## 2. Registry Mutation Inventory

| Mutation | Location | Description |
|---|---|---|
| `amatrix_register_backend(name, backend, overwrite=FALSE/TRUE)` | `backend-registry.R:118` | Adds or replaces a backend in `$backends` |
| `.amatrix_try_register_optional_backend(name)` | `backend-registry.R:38` | Auto-registers from package namespace on first use |
| CPU re-register at load | `zzz.R:17` | `amatrix_register_backend("cpu", …, overwrite=TRUE)` every `onLoad` |
| `amatrix_set_default_policy(policy)` | `policy.R:152` | Changes `$default_policy` |
| `amatrix_set_default_precision(precision)` | `policy.R:246` | Changes `$default_precision` |

There is **no `amatrix_unregister_backend`** function. The only removal path is `overwrite=TRUE` on re-registration, which silently replaces the backend list object in `$backends`.

**Total mutations: 5** (3 on the registry, 2 on policy/precision)

---

## 3. Cache-Mutation Dependency Graph (missing edges)

```
Mutation                          Should invalidate               Does it?
─────────────────────────────────────────────────────────────────────────────
register_backend(name, overwrite) → calibration[name]             NO  ← BUG 1
register_backend(name, overwrite) → backend_health[name]          NO  ← BUG 2
register_backend(name, overwrite) → residency entries for name    NO  ← BUG 3 (partial)
register_backend(name, overwrite) → am_product_plan closures      NO  ← BUG 4
set_default_policy(Y)             → am_product_plan closures      NO  (by design; user responsibility)
set_default_policy(Y)             → calibration                   N/A (calibration not keyed on policy)
─────────────────────────────────────────────────────────────────────────────
```

Missing edges: **4** confirmed, **1** by-design gap.

---

## 4. Scenario Analysis

### S1: register A → plan op X uses A → overwrite-register A (different caps) → plan op X again

**Code path**: `amatrix_backend_plan` → `amatrix_backend_preference` → loop over preferred → `backend$supports(op, x)`.

`backend$supports` is called live from the freshly replaced backend object in `$backends`. The plan does NOT cache results across calls. Each call to `amatrix_backend_plan` re-fetches the backend via `.amatrix_get_backend(candidate_name)` which reads from `$backends` directly.

**Verdict: SAFE.** The planner is stateless — no plan memoization.

---

### S2: register A → calibrate (writes timing data) → overwrite-register A (new caps/impl) → plan uses old thresholds

**Code path**: `amatrix_backend_plan` → `.amatrix_calibration_ok(x, op, "A")` → `.amatrix_load_calibration()` → reads `$calibration$thresholds[["A"]][[op_key]]`.

`amatrix_register_backend` at `backend-registry.R:177` is a single `assign(name, backend, envir = .amatrix_state$backends)`. It does **nothing** to `$calibration`.

**Reproducer sketch**:
```r
amatrix_register_backend("mlx", fast_mlx_v1, overwrite = TRUE)
amatrix_calibrate(backend = "mlx", persist = FALSE)
# $calibration$thresholds$mlx$gemm == 16384  (v1 was slow below this)

amatrix_register_backend("mlx", fast_mlx_v2, overwrite = TRUE)
# v2 is much faster; threshold should be 512
m <- adgeMatrix(matrix(rnorm(1000), 10, 100), preferred_backend = "mlx", precision = "fast")
plan <- amatrix_backend_plan(m, "matmul")
# plan$chosen == "cpu"  WRONG — old threshold from v1 still gates v2
```

**Verdict: CONFIRMED BUG (Severity: High)**. Old calibration data silently governs a new backend implementation. Any re-registration (including auto-register retry with a freshly-loaded package version) will serve stale thresholds.

**Fix location**: `backend-registry.R:172-178` — after `assign(name, backend, envir = $backends)`, add:
```r
if (!is.null(.amatrix_state$calibration)) {
  .amatrix_state$calibration$thresholds[[name]] <- NULL
}
```
The disk copy also needs to be cleaned; the safest fix also deletes the on-disk file for `name` or marks it stale.

---

### S3: register A → set_policy(X) → plan op uses X → set_policy(Y) → plan op again

**Code path**: `amatrix_backend_preference` (in `policy.R:298`) → reads `x@preferred_backend`, `x@policy`, `amatrix_default_policy()` live. `amatrix_default_policy()` reads `$default_policy` which was just updated.

**Verdict: SAFE.** Policy is read fresh on every plan call.

---

### S4: set_backend_capabilities for A (now supports f16) → stale plans

There is no `set_backend_capabilities` function. Capabilities are re-queried via `backend$capabilities()` live in `amatrix_backend_plan`. Any change to the backend object via `overwrite=TRUE` is reflected immediately.

**Verdict: SAFE** (capabilities are not separately cached).

---

### S5: warmup populates fast-path → backend overwrite → stale warmup state

`amatrix_warm` does not write any state into `.amatrix_state`. It only calls `backend$warm()` or runs dummy ops. There is no "fast-path table" populated by warmup.

**Verdict: SAFE** (no warmup cache).

---

### S6: model-cache.R caches per-(model, backend) artifacts → backend re-registration → key collision

Model cache keys are built from `object_id` (e.g. `"chol:<session_id>:am:<counter>"`), not backend names. Keys are session-unique and content-addressed to the matrix data, not to the backend that computed the factor. So re-registering a backend with `overwrite=TRUE` does not cause key collision.

**However**, a subtler issue exists: a QR or Cholesky factor computed by `mlx_v1` may be stored in the model cache under key `"chol:<object_id>"`. After `overwrite=TRUE` registers `mlx_v2` (different implementation, possibly different numerical precision), the stale factor from `mlx_v1` is returned from cache. The key encodes no backend identity.

**Reproducer sketch**:
```r
# mlx_v1 has a bug that produces a wrong Cholesky factor
X <- adgeMatrix(spd_matrix, preferred_backend = "mlx", precision = "fast")
f1 <- am_chol(X)  # cached as "chol:<X@object_id>"

amatrix_register_backend("mlx", mlx_v2_fixed, overwrite = TRUE)  # bug fixed in v2

f2 <- am_chol(X)  # returns f1 from cache — v2 never called
# f2 is wrong because cache key == "chol:<X@object_id>", no backend in key
```

**Verdict: CONFIRMED BUG (Severity: Medium)**. Cache key `chol-factor.R:130` is `paste0("chol:", X@object_id)` — no backend component. Same applies to `svd-factor.R:594` which uses `.amatrix_svd_cache_key`. After a backend overwrite the cached factor may have been computed by the old implementation.

**Fix location**: `chol-factor.R:130` and `svd-factor.R:472`. Cache key must incorporate `x@preferred_backend` (or the actually-chosen backend name at computation time).

---

### Additional: S5b — residency entries survive backend overwrite

When `amatrix_register_backend("mlx", new_impl, overwrite=TRUE)` is called, all existing residency entries in `.amatrix_state$residency` that reference `"mlx"` remain untouched. The next call to `.amatrix_live_resident_backend(x)` will attempt `backend$resident_has(entry$resident_key)` against the **new** backend object. If the new implementation has a fresh (empty) resident store, `resident_has` returns `FALSE` and the entry is treated as dead — so the matrix is transparently re-uploaded.

This is actually **safe** in practice (graceful degradation), but it means:
1. GPU memory tracked by the old backend implementation is never freed via `resident_drop` — the key is just silently abandoned.
2. `amatrix_memory_stats()` may report stale byte counts until `amatrix_gc()` is called.

**Verdict: SUSPECTED BUG (Severity: Low-Medium)** — GPU memory leak on backend overwrite if the underlying device store is shared between old and new backend implementations (same MLX process). No explicit `resident_drop` is called for old keys before the backend is replaced.

---

### Additional: S7 — backend_health not cleared on overwrite

`backend-registry.R:177` replaces the backend but leaves `$backend_health[["mlx"]]` intact. If the old backend was marked `"unhealthy"` (e.g. canary failure), the new backend inherits that mark. `amatrix_backend_status()` will show the new backend as `"unhealthy"` until a fresh `amatrix_backend_health_probe` is run.

**Verdict: CONFIRMED BUG (Severity: Medium)**. Not a data-correctness bug but causes `amatrix_backend_status()` to show misleading health state and could cause user confusion or automated health-check failures.

**Fix location**: `backend-registry.R:177` — after `assign(name, ...)`, add:
```r
if (!is.null(.amatrix_state$backend_health)) {
  .amatrix_state$backend_health[[name]] <- NULL
}
```

---

### Additional: S8 — am_product_plan captures backend_name in closure forever

`amatrix_compile_product` (product-plan.R:167) captures `backend_name` and `lhs_bound` in the `plan_fun` closure at compile time. If the chosen backend is later replaced via `overwrite=TRUE`, all existing compiled plans continue to dispatch to the old backend name. `.amatrix_get_backend(backend_name)` inside `plan_fun` will fetch the **new** backend object but with the old capability/precision assumptions baked into `lhs_bound`'s `@preferred_backend` and `@precision` slots.

More critically, if `lhs_bound` is GPU-resident on the old backend's store, and the new backend has a fresh empty store, the residency will be silently invalidated (safe), but the plan will then fall back to host — a silent performance regression with no warning.

**Verdict: SUSPECTED BUG (Severity: Low)** — no correctness failure (graceful residency invalidation) but silent performance regression and misleading metadata (`meta$backend` in the plan still shows old name).

---

## 5. Confirmed and Suspected Bugs

### Confirmed Bugs

| # | File | Line(s) | Description | Severity |
|---|---|---|---|---|
| B1 | `R/backend-registry.R` | 177 | `amatrix_register_backend(overwrite=TRUE)` does not clear `$calibration$thresholds[[name]]`. Stale per-op thresholds from old backend govern dispatch for the new one. | **High** |
| B2 | `R/backend-registry.R` | 177 | `amatrix_register_backend(overwrite=TRUE)` does not clear `$backend_health[[name]]`. Stale health status (e.g. "unhealthy") persists for a freshly replaced backend. | **Medium** |
| B3 | `R/chol-factor.R` | 130 | Model cache key `paste0("chol:", X@object_id)` has no backend component. Factor computed by old backend implementation is returned after backend overwrite. | **Medium** |
| B4 | `R/svd-factor.R` | 472–476 | `.amatrix_svd_cache_key` similarly encodes no backend identity. SVD factors from old backend survive backend overwrite. | **Medium** |

### Suspected Bugs

| # | File | Line(s) | Description | Severity |
|---|---|---|---|---|
| S1 | `R/backend-registry.R` | 177 | No `resident_drop` called for existing residency entries when backend is overwritten. Old GPU memory keys are abandoned without cleanup. | **Low-Medium** |
| S2 | `R/product-plan.R` | 193–217 | `am_product_plan` closures capture `backend_name` at compile time; silent performance regression after backend overwrite (resident store invalidated, falls back to host with no warning). | **Low** |

### Double-free / double-clear check (S6)

No double-free found. The finalizer in `residency.R:9-52` removes the entry from `$residency` before calling `backend$resident_drop`, so a second GC sweep or explicit `amatrix_gc()` would find no entry and skip. `amatrix_gc()` removes the entry before checking liveness, so the finalizer would find nothing to drop. This is safe.

---

## 6. Summary Statistics

- Caches inventoried: **10**
- Registry mutations: **5**
- Missing invalidation edges: **4 confirmed + 2 suspected**
- Confirmed bugs: **4**
- Suspected bugs: **2**
- Double-free / double-clear issues: **0**
