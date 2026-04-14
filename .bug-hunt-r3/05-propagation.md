# Round 3 Bug Hunt — Hunter 05: Same-Shape Propagation Pass

Generated: 2026-04-14  
Scope: R/*.R — propagation search for each round-2 bug signature in adjacent ops  
Method: grep + static read; all findings marked **inferred** unless noted

---

## Methodology

For each of the 17 round-2 bugs, I distilled a grep-able anti-pattern signature,
searched the full R/ tree, filtered out sites already flagged in round 2, and
cross-referenced against `bd list --status=open` to avoid double-counting.

---

## Bug-by-Bug Analysis

---

### BUG: amatrix-dev — `host_cache_valid` set TRUE unconditionally
**Original site:** `residency.R:97` (now line 112 after refactor to `cache_state$host_cache_valid`)  
**Signature:** `host_cache_valid\s*<-\s*TRUE` anywhere outside the bind function, OR `host_cache_valid` never reset to FALSE after an in-place op.

**Round-2 flagged sites:** `residency.R:112` (the one TRUE-setting site)

**Grep result:** Exactly one TRUE-setting site. Zero FALSE-setting sites exist anywhere. The `cache_state$host_cache_valid <- FALSE` path appears only inside `.amatrix_update_resident_aliases` (residency.R:137), which is called exclusively when a resident key is renamed (alias update). It is NOT called when:

**NEW PROPAGATION CANDIDATES:**

| File:Line | Context | Confidence |
|-----------|---------|------------|
| `R/wrappers.R` — all `_resident` op wrappers (e.g. matmul_resident, crossprod_resident, ewise_resident, solve_resident, chol_resident — ~12 sites) | Each of these allocates `out_key`, calls the backend op, then wraps the result via `.amatrix_bind_resident(value, backend_name, out_key)`. `.amatrix_bind_resident` sets `host_cache_valid=TRUE` on the returned object whose `@x` still holds the pre-op host data (stale). | **inferred** |
| `R/chol-factor.R:244` — `.amatrix_amchol_wrap_resident_result` calls `.amatrix_bind_resident` on a `new_adgeMatrix_deferred` with `@x = NaN` | Deferred object gets `host_cache_valid=TRUE` at bind. If `host_deferred` flag is not set, the fast path at `residency.R:452` returns NaN data. | **inferred** |
| `R/resident-handle.R` — `as_adgeMatrix.resident_handle` path when `defer_host=FALSE` downloads host data then calls `_bind_resident`; that sets `host_cache_valid=TRUE` correctly — but only if the download completed without error. If the download produced a partial matrix, the flag is set TRUE prematurely. | **inferred (low confidence)** |

**Total new candidates: ~14** (12 wrapper sites + chol deferred wrap + handle coerce edge)

---

### BUG: amatrix-4q9 — Double-drop on chol triangular solve
**Original site:** `chol-factor.R:175–197` — error handler drops `out_key` at line 191, then unconditional `resident_has` check at line 195 drops it again.  
**Signature:** `out_key <- _next_resident_key` → `tryCatch(... error = function(e) { try(backend$resident_drop(out_key), ...) ; NULL })` → `if (isTRUE(backend$resident_has(out_key))) { try(backend$resident_drop(out_key), ...) }` — two drop sites for the same key.

**Round-2 flagged sites:** `chol-factor.R:191,196`

**Grep result:** The double-drop pattern (drop in error handler + unconditional drop after tryCatch) appears **only** in `chol-factor.R`. All other resident op wrappers use the single-drop pattern: drop in error handler, then return NULL (no second `resident_has` guard). The chol path is uniquely broken here because the unconditional post-tryCatch drop was apparently added defensively but conflicts with the error-handler drop.

**NEW PROPAGATION CANDIDATES:** None. The adjacent ops (matmul, crossprod, tcrossprod, solve, ewise at wrappers.R:402–708) all use the safe single-error-handler-drop pattern: `error = function(e) { try(backend$resident_drop(out_key), ...) ; NULL }` with no second unconditional drop. The double-drop is isolated to `chol-factor.R`.

**Total new candidates: 0**

---

### BUG: amatrix-8kj — Double-drop in crossprod_weighted / tcrossprod_weighted
**Original site:** `wrappers.R:1283–1299` (crossprod_weighted) and `wrappers.R:1345–1359` (tcrossprod_weighted) — `scaled_key` is dropped inside the tryCatch on the success path and also dropped unconditionally in the failure cleanup block.  
**Signature:** Two-key allocation (`scaled_key` + `out_key`) where `scaled_key` is dropped inside tryCatch success, then also dropped in the cleanup block below.

**Round-2 flagged sites:** wrappers.R:1283–1307, 1345–1367

**Grep result for two-key allocation pattern:**

**NEW PROPAGATION CANDIDATES:**

| File:Line | Context | Confidence |
|-----------|---------|------------|
| `R/wrappers.R:1415–1452` — `tcrossprod_weighted_XY` (or equivalent three-key function) | Allocates `x_scaled_key`, `y_scaled_key`, AND `out_key` (three keys). Cleanup block at line 1446–1448 drops all three unconditionally. If the success path already consumed/dropped any of these, the cleanup is a double-drop. The pattern is structurally identical to amatrix-8kj. | **inferred** |
| `R/wrappers.R:627` — `ewise_resident` wrapper | Allocates `out_key`, error drops it, but then `.amatrix_cleanup_temp_resident(list(lhs, rhs), backend_name)` runs for the lhs/rhs — if lhs or rhs was a freshly promoted temp key AND out_key was the same counter value (impossible due to incrementing counter, but the cleanup ordering is fragile), race exists. Lower confidence. | **inferred (low)** |

**Total new candidates: 2**

---

### BUG: amatrix-aul — ~40 alloc sites use manual try-drop instead of on.exit
**Original site:** `wrappers.R` ~40 sites — `out_key` allocated before tryCatch, cleaned up manually after, with no `on.exit` guard.  
**Signature:** `out_key <- .amatrix_next_resident_key(...)` followed by `tryCatch(...)` followed by `try(backend$resident_drop(out_key), ...)` without an `on.exit` wrapping the entire block.

**Round-2 flagged sites:** wrappers.R representative sites ~40

**Grep result — additional files NOT in wrappers.R:**

**NEW PROPAGATION CANDIDATES:**

| File:Line | Context | Confidence |
|-----------|---------|------------|
| `R/backend-planning.R:378–400` | `out_key <- .amatrix_next_resident_key(resident_backend_name)` at line 378; tryCatch at line 379; cleanup `try(backend$resident_drop(out_key), ...)` at line 400. No `on.exit` guard. Identical pattern to the wrappers.R family. | **inferred** |
| `R/chol-factor.R:180–199` | `out_key <- .amatrix_next_resident_key(backend_name)` at 180; tryCatch at 181; no `on.exit`. Already flagged as amatrix-4q9 but also qualifies here as a missing `on.exit`. | **inferred** (overlaps amatrix-4q9) |
| `R/irlba.R:503–530` | `out_key <- .amatrix_next_resident_key(backend_name)` (around line 503); tryCatch `error = function(e) FALSE`; success check `if (!success) { if (isTRUE(backend$resident_has(out_key))) backend$resident_drop(out_key) }` — uses `resident_has` guard but NOT `on.exit`. If the tryCatch throws AND `resident_has` itself throws, key leaks. | **inferred** |

**Total new candidates: 3 (backend-planning.R:378 is NEW and most actionable)**

---

### BUG: amatrix-86l — S4 Math group missing for adgeMatrix/adgCMatrix
**Original site:** No `setMethod("Math", ...)` anywhere.  
**Signature:** `setMethod("Math", ...)` absent; also check `setMethod("Summary", ...)`.

**Round-2 flagged sites:** Both adgeMatrix and adgCMatrix in 01-s4-dispatch-grid.md

**Grep result:** Zero hits for `setMethod("Math"` and zero hits for `setMethod("Summary"` in all R files.

**NEW PROPAGATION CANDIDATES:**

| Missing method | Class | Consequence | Confidence |
|----------------|-------|-------------|------------|
| `setMethod("Summary", "adgeMatrix", ...)` | adgeMatrix | `sum(A)`, `max(A)`, `prod(A)`, `any(A)`, `all(A)` bypass GPU; go through inherited Matrix Summary method. Return type is scalar (correct) but residency is silently broken and GPU path never fires. `any(A > 0)` where `A > 0` returned an adgeMatrix will materialize unnecessarily. | **inferred** |
| `setMethod("Summary", "adgCMatrix", ...)` | adgCMatrix | Same issue. | **inferred** |
| `setMethod("Math", "adgeMatrix", ...)` | adgeMatrix | Already flagged amatrix-86l but propagation: any pipeline using `sqrt`, `exp`, `log`, `abs` on a GPU-resident adgeMatrix silently demotes to dgeMatrix. | **inferred** (amatrix-86l) |
| `setMethod("Math", "adgCMatrix", ...)` | adgCMatrix | Same. | **inferred** (amatrix-86l) |

**Distinct new candidate: Summary group missing for both classes — not filed in round 2.**  
**Total new candidates: 2 (Summary,adgeMatrix and Summary,adgCMatrix)**

---

### BUG: amatrix-j5a — diag<- replacement missing for adgeMatrix
**Original site:** No `setReplaceMethod("diag<-", ...)` anywhere.  
**Signature:** `setReplaceMethod("diag<-", ...)` absent for adgeMatrix and adgCMatrix.

**Round-2 flagged sites:** Both classes in 01-s4-dispatch-grid.md

**Grep result:** `setReplaceMethod` exists only for `[`, `[<-`, and `dimnames<-`. No `diag<-` replacement for either class.

**NEW PROPAGATION CANDIDATES:**

| Missing method | Class | Consequence | Confidence |
|----------------|-------|-------------|------------|
| `setReplaceMethod("diag<-", "adgCMatrix", ...)` | adgCMatrix | `diag(S) <- v` for a sparse adgCMatrix fires Matrix's `diag<-,Matrix-method`, returning a plain dgCMatrix. All amatrix slots lost. Round-2 report mentioned adgeMatrix; adgCMatrix was not explicitly confirmed as a separate issue. | **inferred** |

**Total new candidates: 1 (adgCMatrix diag<- — not explicitly filed as its own issue)**

---

### BUG: amatrix-jnd — kronecker not intercepted
**Original site:** No `setMethod("kronecker", ...)` anywhere.  
**Signature:** `setMethod("kronecker", ...)` absent.

**Round-2 flagged sites:** Both adgeMatrix and adgCMatrix in 01-s4-dispatch-grid.md

**Grep result:** Zero hits for `setMethod("kronecker"` in all R files.

**NEW PROPAGATION CANDIDATES:**

| Missing method | Class | Consequence | Confidence |
|----------------|-------|-------------|------------|
| `kronecker(adgCMatrix, ANY)` | adgCMatrix | Same as adgeMatrix — fires base::kronecker, returns plain matrix. The round-2 issue amatrix-jnd mentions both classes but only the dense case was the primary example. | **inferred** |
| `kron()` (wrappers.R:2860) exists as a custom function but does NOT call `setMethod("kronecker", ...)` | Both | Users calling the standard `kronecker()` generic get demoted to plain matrix. `kron()` is the workaround but it's not discoverable from `?kronecker`. | confirmed propagation |

**Total new candidates: 1 (explicit adgCMatrix gap confirmation)**

---

### BUG: amatrix-cth — Deferred adgeMatrix dangles after release
**Original site:** `product-plan.R:232–261` — `lhs_bound` in compiled plan closure holds deferred adgeMatrix; GPU release makes plan's fallback operate on NaN host data.  
**Signature:** `new_adgeMatrix_deferred` object passed to a function that stores it in a closure or long-lived structure WITHOUT checking liveness before use; fallback reads `@x` = NaN.

**Round-2 flagged sites:** product-plan.R:232–261

**Grep result for `new_adgeMatrix_deferred` usage:**

**NEW PROPAGATION CANDIDATES:**

| File:Line | Context | Confidence |
|-----------|---------|------------|
| `R/chol-factor.R:236–244` — `.amatrix_amchol_wrap_resident_result` creates `new_adgeMatrix_deferred` then calls `.amatrix_bind_resident` | The deferred object's `@x` = NaN sentinel. If the GPU buffer for `resident_key` is later freed (e.g. via `amatrix_gc()` or backend overwrite), the `factor_obj` slot in `amChol` now holds a dangling deferred object. `chol_solve` will call `.amatrix_resident_triangular_solve` → if resident key is gone, falls back to `factor@factor` (the host factor), but `factor@factor` is `matrix(numeric(0), 0, 0)` for GPU-only factors (chol-factor.R:147). Result: zero-row matrix used as triangular factor → wrong solve result silently. | **inferred** |
| `R/svd-factor.R` — SVD result objects store `ut_am` (an adgeMatrix, possibly deferred) as a slot in `amSVD` | If the GPU buffer backing `ut_am` is freed before `svd_project` or `svd_reconstruct` is called, the project/reconstruct will materialize NaN data. Same pattern as product-plan.R G2. | **inferred** |

**Total new candidates: 2 (chol-factor deferred factor_obj, svd-factor ut_am slot)**

---

### BUG: amatrix-fqh — Model cache keys miss backend identity
**Original site:** `chol-factor.R:136` key = `paste0("chol:", X@object_id)` and `svd-factor.R:472` `.amatrix_svd_cache_key` — no backend component.  
**Signature:** `paste0("chol:"` or `paste0("svd:"` or any model cache key construction that omits the backend name or precision.

**Round-2 flagged sites:** chol-factor.R:136, svd-factor.R:472

**Grep result:**

**NEW PROPAGATION CANDIDATES:**

| File:Line | Context | Confidence |
|-----------|---------|------------|
| `R/chol-factor.R:46` — `.amatrix_amchol_factor_matrix` builds `cache_key <- paste0("chol:", factor@source_id)` to UPDATE a cached amChol | This is a write-back path: it reads the cache, updates the `@factor` slot, and re-stores. The key uses `source_id` (the original matrix object_id), not backend. If a factor was computed by mlx_v1, cached, and then mlx_v2 is registered (backend overwrite), this path would overwrite the stale cache entry with new host data from the wrong factor. Separate from the primary miss because this is a cache UPDATE path, not lookup. | **inferred** |
| `R/models-lm.R` — QR/XtX model cache keys | `models-lm.R` uses model_cache for QR and XtX factors (grep confirms `model_cache` usage). The cache key construction in `models-lm.R` for QR factors likely follows the same pattern (object_id only, no backend). Grep shows no explicit `paste0("qr:"` or `paste0("lm:"` string — keys may be constructed differently, but the risk exists. | **inferred (low — needs targeted read)** |

**Total new candidates: 1 confirmed (chol-factor.R:46 update path)**

---

### BUG: amatrix-dmy — lm_loo_cv mixed-precision
**Original site:** `qr-downdate.R:143–145` — LOO loop uses original `X` (possibly float64) but QR was built from `X_am` (possibly float32 via precision="fast").  
**Signature:** Any function that builds a factor from an adgeMatrix `X_am` but then accesses numeric data from the original plain `X` for residual computation.

**Round-2 flagged sites:** qr-downdate.R:143,145

**Grep result:**

**NEW PROPAGATION CANDIDATES:**

| File:Line | Context | Confidence |
|-----------|---------|------------|
| `R/qr-downdate.R:131` — `qr_downdate.amQR` calls `am_qr(.amatrix_qr_arg(X_sub))` where `X_sub <- X[-row_idx, , drop=FALSE]` | If `X` is plain matrix, `X_sub` is plain matrix → `am_qr` wraps to adgeMatrix with defaults (may change precision). The downdate QR may use different precision than `qr_full`. Same mixed-precision hazard as the LOO loop. | **inferred** |
| `R/models-lm.R` — `wls_fit` uses `X_arg` (adgeMatrix, possibly float32) but `weights` (plain numeric, float64) and combines in residual computation | If `X_arg` is float32 and `weights` is float64, the weighted residual `(y - X %*% beta) * sqrt(w)` mixes precisions. Not as severe as LOO but the same class of precision-mixing bug. | **inferred** |

**Total new candidates: 2**

---

### BUG: amatrix-qm2 — Mixed-backend op corrupts shared object's registry binding
**Original site:** `residency.R:84–106` / `wrappers.R` — `.amatrix_bind_resident` overwrites the registry entry for an object without checking if it already has a binding on a DIFFERENT backend; the old backend's buffer is abandoned.  
**Signature:** `.amatrix_bind_resident(x, ...)` called on an object `x` that already has a live registry entry under a different backend name.

**Round-2 flagged sites:** residency.R:84–106 (the bind function itself)

**Callers of `.amatrix_bind_resident` that may call it on already-resident objects:**

**NEW PROPAGATION CANDIDATES:**

| File:Line | Context | Confidence |
|-----------|---------|------------|
| `R/wrappers.R:162–170` — `am_transpose` resident path: creates new object, binds it | This creates a fresh object each time (new object_id), so no overwrite of existing binding. SAFE. | — |
| `R/wrappers.R:374–378` — sparse matmul result bind | Result is always a fresh deferred object. SAFE. | — |
| `R/resident-handle.R` — `as_adgeMatrix.resident_handle` calls `_bind_resident` on a newly created object | Fresh object, no pre-existing binding. SAFE. | — |
| `R/bind-resident.R:~67` — `amatrix_bind_resident` (the exported function) | This IS the public API for explicitly binding a matrix to a backend. If called twice on the same object (once for mlx, once for arrayfire), the second call silently overwrites. No guard. Any user-facing code or internal code that calls this twice is the propagation site. | **inferred** |
| `R/backend-planning.R:~350–400` — `amatrix_dispatch_op` calls `_prepare_resident_arg` which may call `_bind_resident` on the RHS | If x is on mlx and y is on arrayfire, `_prepare_resident_arg(y, "mlx", promote_amatrix=TRUE)` uploads y to mlx and calls `_bind_resident(y, "mlx", new_key)` — overwriting y's arrayfire binding. This is the amatrix-qm2 propagation site confirmed by round-2 report (Pattern F1). Already in amatrix-qm2. | (already flagged) |

**Total new candidates: 1 (bind-resident.R public API double-bind, distinct from internal call)**

---

### BUG: amatrix-ubq — Backend overwrite leaves stale calibration thresholds
**Original site:** `backend-registry.R:184` — `assign(name, backend, ...)` with no calibration invalidation.  
**Signature:** `assign(name, backend, envir = .amatrix_state$backends)` without clearing `$calibration$thresholds[[name]]`.

**Round-2 flagged sites:** backend-registry.R:184

**Grep result:** Only ONE `assign(name, backend, envir = .amatrix_state$backends)` site exists (line 184). The `zzz.R:17` site calls `amatrix_register_backend("cpu", ..., overwrite=TRUE)` on every `.onLoad` — meaning EVERY package reload reregisters cpu with no calibration invalidation. If the user ran `amatrix_calibrate(backend="cpu")` before the reload, the stale cpu thresholds govern dispatch after reload.

**NEW PROPAGATION CANDIDATES:**

| File:Line | Context | Confidence |
|-----------|---------|------------|
| `R/zzz.R:17` — `amatrix_register_backend("cpu", .amatrix_cpu_backend(), overwrite=TRUE)` | Called every `.onLoad`. If calibration data exists from a prior session (loaded from disk via `.amatrix_load_calibration`), the cpu backend is replaced but its calibration thresholds remain. On a machine where the cpu backend was slow (triggering low thresholds), reloading the package with a new CPU backend implementation silently inherits wrong thresholds. | **inferred** |

**Total new candidates: 1**

---

### BUG: amatrix-2nh — Backend overwrite leaves stale backend_health
**Original site:** `backend-registry.R:184` — same assign site as amatrix-ubq; `$backend_health[[name]]` not cleared.  
**Signature:** Same as amatrix-ubq — no invalidation of `$backend_health` on overwrite.

**Round-2 flagged sites:** backend-registry.R:184

**NEW PROPAGATION CANDIDATES:**

| File:Line | Context | Confidence |
|-----------|---------|------------|
| `R/zzz.R:17` — same `.onLoad` reregistration of cpu | If cpu was marked "unhealthy" in a prior session (e.g. canary failure during testing), reloading the package leaves `$backend_health[["cpu"]] = "unhealthy"`. `amatrix_dispatch_op` won't skip cpu (it uses health for GPU backends only, per current logic), but `amatrix_backend_status("cpu")` would show stale "unhealthy" state. | **inferred** |
| `R/backend-registry.R:388–395` — `.amatrix_backend_health_mark` is called from `amatrix_dispatch_op` error handlers at lines 388 and 422 | These mark a backend "unhealthy" on any runtime error. There is no TTL or re-probe on the health mark — once unhealthy, only `amatrix_backend_health_probe()` can clear it. If a transient GPU error marks mlx "unhealthy", it stays unhealthy for the session. This is a stale-health propagation risk orthogonal to the overwrite path. | **inferred** |

**Total new candidates: 2**

---

### BUG: amatrix-3ka — `.amatrix_bind_resident` leaks prior key on rebind
**Original site:** `residency.R:90–116` — `_bind_resident` overwrites the registry entry without releasing the prior resident key.  
**Signature:** Any code path that calls `_bind_resident` on an object that already has a live entry in `.amatrix_state$residency`.

**Round-2 flagged sites:** residency.R:90–116

**NEW PROPAGATION CANDIDATES:**

| File:Line | Context | Confidence |
|-----------|---------|------------|
| `R/irlba.R:155,631,634,699,702` — irlba uses `on.exit` to clean up GPU keys | The `on.exit` drops keys via `backend$resident_drop` directly, not via `_release_resident`. If the irlba function previously called `_bind_resident` on an intermediate object, the registry entry for that object is never cleared — the direct `resident_drop` frees the device buffer but leaves a dangling entry in `$residency`. Next call to `_live_resident_backend` on that object will call `resident_has` (which returns FALSE since the buffer is gone), so the dangling entry is functionally inert but accumulates. | **inferred** |
| `R/bind-resident.R` — `amatrix_bind_resident` (exported) | If a user calls this on an object that was already bound (e.g. re-uploading after a precision change), the prior key is silently leaked. No user-visible warning. | **inferred** (same as amatrix-3ka, propagated to public API) |

**Total new candidates: 2**

---

### BUG: amatrix-uu2 — dispatch cold-path tryCatch swallows GPU error class on fallback re-dispatch
**Original site:** `backend-planning.R:412–430` — the cold-path `tryCatch` catches GPU errors, calls `fallback()` (re-dispatch to CPU), but the original error class is lost. `conditionMessage(e)` is logged but `class(e)` is discarded.  
**Signature:** `tryCatch(... error = function(e) { .amatrix_log_fallback(...) ; return(fallback()) })` — the original error object `e` is not re-signalled or wrapped.

**Round-2 flagged sites:** backend-planning.R:412–430

**NEW PROPAGATION CANDIDATES:**

| File:Line | Context | Confidence |
|-----------|---------|------------|
| `R/backend-planning.R:379–400` — the RESIDENT path tryCatch (lines 379–393) | Error handler calls `.amatrix_log_fallback` and `.amatrix_backend_health_mark` then returns NULL. The NULL then triggers a second dispatch path. This is the SAME swallow pattern — the GPU error class is discarded. | **inferred** |
| `R/svd-factor.R:250–262` — `.amatrix_subspace_compile_operator` wraps `amatrix_compile_product` in `tryCatch(... error = function(e) NULL)` | If `amatrix_compile_product` raises an `amatrix_error_backend_unavailable` (or similar typed condition — once amatrix-6m9 is fixed), that class is swallowed to NULL. The subspace SVD then silently falls back to a CPU operator with no logged event and no condition signalled. | **inferred** |
| `R/models-lm.R:579–591` — GPU broadcast_ewise tryCatch | Already flagged in round 2 (02-antipattern-greps.md) as a key-leak bug (amatrix-aul propagation), but ALSO swallows error class — dual propagation. | **inferred** |

**Total new candidates: 3** (resident path + svd compile operator + models-lm overlap)

---

### BUG: amatrix-833 — Subspace SVD stop() inside tryCatch silently swallowed
**Original site:** `svd-factor.R:379` and `svd-factor.R:411` — `stop("subspace SVD did not discover a usable range space")` and `stop("projected core is numerically rank-deficient")` are called INSIDE `.amatrix_subspace_svd`, which is called from `.amatrix_subspace_compile_operator` at svd-factor.R:252–261 inside `tryCatch(... error = function(e) NULL)`.  
**Signature:** `stop(...)` calls inside functions that are wrapped by callers in `tryCatch(... error = function(e) NULL)`.

**Round-2 flagged sites:** svd-factor.R:379,411

**Grep for other stop() calls inside functions called from error=NULL tryCatch:**

**NEW PROPAGATION CANDIDATES:**

| File:Line | Context | Confidence |
|-----------|---------|------------|
| `R/chol-factor.R` — `.amatrix_resident_triangular_solve` contains no `stop()` internally; its tryCatch at line 181–193 would swallow any backend-raised condition. But the backend function `solve_triangular_resident` may itself call `stop()`. If the backend raises a typed condition (future state post amatrix-6m9), the `tryCatch(error=function(e){ try(drop); NULL })` at line 181 swallows it. | **inferred (low)** |
| `R/irlba.R:303,367,400` — `.Call("amatrix_block_reorth_bridge", ...)` C-level errors swallowed | Already flagged in round 2 as a separate issue. The C code may call `Rf_error()` inside a section that is wrapped in `tryCatch(error=function(e) FALSE/NULL)`. Any typed condition class from C is discarded. | **inferred** (already noted in r2) |
| `R/qr.R:453` — `backend$resident_materialize(key)` inside `tryCatch(error=function(e) NULL)` | If materialization fails with a typed backend condition, caller receives NULL silently. | **inferred** (already noted in r2) |

**Total new candidates: 1 (chol triangular-solve backend error swallow)**

---

### BUG: amatrix-hjj — Resident-path tryCatch returns NULL silently, no amatrix_fallback condition class
**Original site:** Multiple resident-path wrappers return NULL on error with no condition signalled; `amatrix_dispatch_op` interprets NULL as "resident path unavailable" and transparently falls back to CPU.  
**Signature:** `result <- tryCatch(backend$..., error = function(e) NULL)` in a resident-path function where the NULL is silently interpreted as "try fallback" with no structured event or condition class emitted.

**Round-2 flagged sites:** Multiple wrappers.R resident functions

**NEW PROPAGATION CANDIDATES (sites where NULL fallback is indistinguishable from "not supported"):**

| File:Line | Context | Confidence |
|-----------|---------|------------|
| `R/wrappers.R:402–408` — `_try_resident_matmul` | `error = function(e) { try(drop); NULL }` — GPU OOM returns the same NULL as "backend doesn't support matmul_resident". Caller cannot distinguish. | **inferred** |
| `R/wrappers.R:506–514` — `_try_resident_crossprod` | Same pattern. | **inferred** |
| `R/wrappers.R:576–585` — `_try_resident_tcrossprod` | Same pattern. | **inferred** |
| `R/wrappers.R:670–691` — `_try_resident_solve` | Same pattern. | **inferred** |
| `R/wrappers.R:697–708` — `_try_resident_chol` | Same pattern. | **inferred** |
| `R/backend-planning.R:379–400` — dispatch resident path | `result` is NULL on error; this already logs via `_log_fallback` (partial mitigation), but the error class is still discarded. | **inferred** |

**Total new candidates: 5** (the five resident op wrappers in wrappers.R that were not in the original amatrix-hjj filing)

---

## Summary Table

| Parent Bug | Signature | Round-2 Sites | NEW Candidates | Files |
|------------|-----------|---------------|----------------|-------|
| amatrix-dev | `host_cache_valid=TRUE` never reset in op wrappers | residency.R:112 | ~14 | wrappers.R (all resident ops), chol-factor.R:244 |
| amatrix-4q9 | double-drop (error handler + post-tryCatch) | chol-factor.R:191,196 | 0 | — |
| amatrix-8kj | two-key alloc, success-path drop + cleanup drop | wrappers.R:1283,1345 | 2 | wrappers.R:1415–1452 |
| amatrix-aul | `out_key` alloc without `on.exit` | wrappers.R ~40 sites | 3 | backend-planning.R:378, irlba.R:503 |
| amatrix-86l | `setMethod("Math",...)` absent | both classes | 2 | Summary group — new for both classes |
| amatrix-j5a | `setReplaceMethod("diag<-",...)` absent | adgeMatrix | 1 | adgCMatrix gap |
| amatrix-jnd | `setMethod("kronecker",...)` absent | both classes | 1 | adgCMatrix explicit confirmation |
| amatrix-cth | deferred adgeMatrix in long-lived structure | product-plan.R:232 | 2 | chol-factor.R:236, svd-factor.R ut_am |
| amatrix-fqh | cache key omits backend name | chol-factor.R:136, svd-factor.R:472 | 1 | chol-factor.R:46 (update path) |
| amatrix-dmy | mixed-precision X vs X_am | qr-downdate.R:143 | 2 | qr-downdate.R:131, models-lm.R wls |
| amatrix-qm2 | `_bind_resident` overwrites without releasing prior | residency.R:90 | 1 | bind-resident.R public API |
| amatrix-ubq | no calibration clear on overwrite | backend-registry.R:184 | 1 | zzz.R:17 (.onLoad cpu re-register) |
| amatrix-2nh | no health clear on overwrite | backend-registry.R:184 | 2 | zzz.R:17, health mark no TTL |
| amatrix-3ka | prior key leaked on rebind | residency.R:90 | 2 | irlba.R direct drop, bind-resident.R API |
| amatrix-uu2 | GPU error class swallowed in fallback | backend-planning.R:412 | 3 | resident path, svd compile op, models-lm |
| amatrix-833 | stop() inside function called from error=NULL | svd-factor.R:379,411 | 1 | chol triangular-solve backend errors |
| amatrix-hjj | NULL fallback indistinguishable from unsupported | multiple wrappers | 5 | wrappers.R matmul/crossprod/tcrossprod/solve/chol |

**Total new propagation candidates: ~43**

---

## Top 5 Most Actionable New Candidates

### 1. `backend-planning.R:378` — Missing `on.exit` for `out_key` (propagation of amatrix-aul)
**File:** R/backend-planning.R:378–400  
**Why top-ranked:** `amatrix_dispatch_op` is the central dispatch hub called by virtually every exported operation. The `out_key` allocated at line 378 in the resident path has no `on.exit` guard. If `.amatrix_cleanup_temp_resident` at line 395 or the `backend$resident_drop` at line 400 both fail (e.g. backend itself errored), the key leaks permanently. This is the highest-traffic code path in the package.

### 2. `wrappers.R:1415–1452` — Three-key double-drop in weighted crossprod (propagation of amatrix-8kj)
**File:** R/wrappers.R:1415–1452  
**Why top-ranked:** Structurally identical to the two confirmed double-drops (amatrix-8kj). Three keys (`x_scaled_key`, `y_scaled_key`, `out_key`) allocated; cleanup drops all three unconditionally at lines 1446–1448 even if the success path already consumed/dropped some. Direct code sibling of the filed bugs.

### 3. `chol-factor.R:236–244` — Deferred factor_obj becomes dangling (propagation of amatrix-cth)
**File:** R/chol-factor.R:236–244  
**Why top-ranked:** `.amatrix_amchol_wrap_resident_result` creates a `new_adgeMatrix_deferred` and stores it in `amChol@factor_obj`. When the GPU buffer is freed (GC, overwrite, or explicit `amatrix_gc()`), the `factor_obj` becomes a dangling deferred object whose `@x = NaN`. `chol_solve` then falls back to `factor@factor` (the host matrix), which is `matrix(numeric(0), 0, 0)` for GPU-only factors — producing silently wrong solve results.

### 4. `zzz.R:17` — `.onLoad` reregisters cpu with no calibration/health invalidation (propagation of amatrix-ubq + amatrix-2nh)
**File:** R/zzz.R:17  
**Why top-ranked:** Every `library(amatrix)` or `devtools::load_all()` in a session that already has calibration data will leave stale thresholds governing the cpu backend. This affects ALL users doing interactive development and CI runs with `--no-save` session caching. The fix is a two-line addition to `amatrix_register_backend`.

### 5. `svd-factor.R:250–262` — Subspace compile operator swallows typed errors (propagation of amatrix-uu2 + amatrix-833)
**File:** R/svd-factor.R:250–262  
**Why top-ranked:** `.amatrix_subspace_compile_operator` wraps `amatrix_compile_product` in `tryCatch(... error=function(e) NULL)`. The two `stop()` calls at svd-factor.R:379 and 411 (amatrix-833's confirmed sites) are called inside `.amatrix_subspace_svd` which is called from within `amatrix_compile_product`'s scope. Any GPU OOM, backend unavailability, or numerical failure inside subspace SVD is silently converted to NULL → CPU fallback with no logged event and no condition class. Users see wrong/slow results with no diagnostic.

---

## Cross-Reference: Not Double-Counting

The following open issues already cover adjacent territory; new candidates above do NOT duplicate them:

- `amatrix-6m9` covers the full unclassed-stop audit — new candidates here are about PROPAGATION of the pattern to specific ops, not the general class.
- `amatrix-cng` covers `expect_error(class=)` test coverage — not duplicated.
- `amatrix-4rt` covers irlba GPU upload cleanup — the irlba `on.exit` gap noted above (irlba.R:503) is distinct: it's about the `out_key` NOT being guarded by `on.exit`, not about the upload path.
- `amatrix-5ni` covers `am_sweep_inplace` specifically — not duplicated by the general missing `on.exit` candidates above.

