# Round 2 Bug Hunt — Anti-Pattern Grep Sweep

Generated: 2026-04-14
Scope: `R/*.R` (all source files in amatrix package)

---

## Grep 1 — Unclassed `stop()`

**Pattern:** `stop("` in R/*.R (no single-quoted hits found)
**Total hits:** ~90 sites

All `stop()` calls use bare string messages with no `class=` argument and no `rlang::abort()`. This means callers cannot `tryCatch(expr, amatrix_error_class = function(e) ...)` by class — they must match by fragile regex on `conditionMessage()`.

| File:Line | Message | Severity |
|-----------|---------|----------|
| `R/residency.R:400` | `"deferred adgeMatrix lost its GPU resident data"` | **HIGH** — critical path; callers need to distinguish from generic errors |
| `R/residency.R:436` | `"resident backend returned an unsupported dense materialization type"` | **HIGH** — same |
| `R/resident-handle.R:109` | `"resident_handle is no longer active"` (also missing `call. = FALSE`) | **HIGH** — stack trace leaks to user; no class for catch |
| `R/resident-handle.R:170,199` | backend sweep failure stops | **HIGH** — in-place op failures uncatchable by class |
| `R/models-lm.R:39` | `"x must be a dense matrix-like object"` (missing `call. = FALSE`) | **MED** |
| `R/models-lm.R:612,616,620` | weight validation stops | **MED** |
| `R/qr.R:730` | `"singular matrix 'a' in solve"` | **MED** — callers cannot distinguish singularity from other solve failures |
| `R/chol-factor.R:127,435,499,516,563,589,701,788,820` | chol/solve validation | **MED** |
| `R/irlba.R:76` | `"Package 'irlba' is required"` — no class | **LOW** |
| `R/constructors.R:67,88` | constructor guards | **LOW** |
| `R/sinkhorn.R:49,108,126,158,161,164,167,170,173,265` | all sinkhorn validation | **LOW** |
| `R/svd-factor.R:379,411,566,570,580,583,682,734` | SVD path errors | **LOW** |
| `R/wrappers.R:1098,1271,1333,1399,1404,1482,1494,1556` | wrapper arg checks | **LOW** |
| `R/backend-registry.R:122,144,147,150,157,160,163` | registry validation | **LOW** |

**Count: ~90 unclassed stops. ~10 HIGH, ~20 MED, ~60 LOW.**

---

## Grep 2 — `tryCatch` swallowing to `NULL`

**Pattern:** `error = function(...) NULL`
**Total hits:** ~80 sites

### Safe / intentional (majority, ~55 sites)
Pattern: `.amatrix_get_backend(name)` probe — if backend not loaded, returns NULL and falls back. This is deliberate. Files: residency.R, chol-factor.R, svd-factor.R, qr.R, policy.R, resident-handle.R, backend-*.R, product-plan.R, etc.

### Actual bugs — silent swallow of real errors

| File:Line | Context | Hypothesis | Severity |
|-----------|---------|------------|----------|
| `R/chol-factor.R:175–191` | `solve_triangular_resident(...)` on error drops `out_key` (line 185), but lines 189–190 unconditionally try to drop `out_key` again — **double-drop** if error path fires | **HIGH** |
| `R/models-lm.R:579–589` | `out_key` pre-allocated at 579; tryCatch at 580 swallows broadcast_ewise error; failure falls through to `try(backend$resident_drop(out_key))` at 589 but only if `!is.null(result)` branch not taken — on error `out_key` leaks if drop at 589 is skipped | **HIGH** |
| `R/wrappers.R:1283–1299` | `scaled_key` + `out_key` pre-allocated; inside tryCatch `scaled_key` is dropped at line 1288 on success; on error fallback at 1297 tries to drop `scaled_key` again — **double-drop race** | **HIGH** |
| `R/wrappers.R:1345–1359` | Same pattern in `tcrossprod_weighted` — identical double-drop race | **HIGH** |
| `R/irlba.R:303,367,400` | `.Call("amatrix_block_reorth_bridge", ...)` C-level errors swallowed; C code may partially mutate internal Lanczos basis state before throwing | **MED** |
| `R/svd-factor.R:252–261` | `amatrix_compile_product(...)` compile failure silently becomes CPU fallback; OOM or precision error invisible | **MED** |
| `R/qr.R:453` | `backend$resident_materialize(key)` error → NULL; caller receives NULL Q/R factor silently | **MED** |
| `R/product-plan.R:101,135` | `resident_materialize` error → NULL; matmul plan silently falls back | **MED** |
| `R/models-lm.R:584` | GPU broadcast_ewise error → NULL without logging; WLS silently falls back to CPU, performance cliff invisible | **LOW** |
| `R/backend-calibration.R:326–332` | warmup loop — best-effort; intentional | OK |
| `R/backend-warmup.R:73,120` | warmup ops — best-effort; intentional | OK |
| `R/memory-stats.R:47–52` | best-effort memory query — intentional | OK |

---

## Grep 3 & 4 — `host_cache_valid` flag management

**Pattern:** `host_cache_valid <- TRUE` / `host_cache_valid <- FALSE` / `host_cache_valid` reads
**Total hits:** SET TRUE: 1 site; SET FALSE: 0 sites; READ: 2 sites

| File:Line | Direction | Context | Hypothesis | Severity |
|-----------|-----------|---------|------------|----------|
| `R/residency.R:97` | SET TRUE | `.amatrix_bind_resident()` — after uploading host data to GPU | Correct: host and device in sync at bind time | OK |
| `R/residency.R:406` | READ | `amatrix_materialize_dense()` cache-hit fast path | Correct: returns cached host data if flag is set | OK |

**Critical gap: `host_cache_valid` is NEVER set to FALSE anywhere in the codebase.**

Consequences:
1. `resident-handle.R:155–168` — `broadcast_ewise_resident_inplace_key` mutates device buffer in place. The `adgeMatrix` that spawned this handle has `host_cache_valid = TRUE` (set at bind time). Next call to `amatrix_materialize_dense()` hits the cache-hit branch at line 406 and returns **stale host data**.
2. `resident-handle.R:178–205` — `broadcast_ewise_resident_key` allocates a new resident key and updates `h$resident_key`, but the originating `adgeMatrix`'s `finalizer_env$host_cache_valid` is never cleared.
3. Any backend op that writes to a new key and rebinds the result correctly: the intermediate handle that was the *source* of the operation still has `host_cache_valid = TRUE` even though the binding for that object ID has been replaced.

**Severity: HIGH — stale host cache reads after any in-place resident mutation.**

---

## Grep 5 — Finalizer attach order bugs

**Pattern:** `reg.finalizer`
**Total hits:** 3 sites

| File:Line | Assessment | Severity |
|-----------|------------|----------|
| `R/residency.R:9` | `e$object_id <- object_id` at line 8 THEN `reg.finalizer(e, ...)` at line 9 — correct order | OK |
| `R/qr.R:19` | `state$q`, `state$q_key`, `state$backend_ops` set at lines 15–17 THEN `reg.finalizer` at 19 — correct | OK |
| `R/resident-handle.R:95` | All `h$*` fields populated at lines 86–93 THEN `reg.finalizer(h, ...)` at 95 — correct | OK |

No finalizer ordering bugs found.

---

## Grep 6 — Assign-to-self / `<<-` / S4 slot mutation in-place

**Pattern:** `<<-` and `x@...` mutations inside methods
**Total hits:** `<<-` at 1 site; `x@slot` reads pervasive but no write-back mutations found

| File:Line | Code | Hypothesis | Severity |
|-----------|------|------------|----------|
| `R/backend-explain.R:121` | `s <<- c(s, paste0(...))` | Local closure accumulator — standard idiom, not S4 mutation | OK |
| `R/residency.R:97` | `x@finalizer_env$host_cache_valid <- TRUE` | Mutates through reference env — intentional design, safe | OK |

No S4 slot direct mutation bugs found. The `x@finalizer_env` pattern is reference-env based and is the intended sharing mechanism.

---

## Grep 7 — `as.matrix()` without dimnames preservation

**Pattern:** `as.matrix(` in R/*.R
**Total hits:** ~160 sites

Most are safe. The core materialization path (`amatrix_materialize_dense` → `new("dgeMatrix", ..., Dimnames = x@Dimnames, ...)`) correctly preserves dimnames.

**Actual dimnames-drop bugs:**

| File:Line | Code | Hypothesis | Severity |
|-----------|------|------------|----------|
| `R/residency.R:396,403` | `fenv$host_x <- ... as.matrix(mat)` then `_dense_base(fenv$host_x)` | If GPU backend returns plain `matrix` with no dimnames, `fenv$host_x` lacks names. `.amatrix_dense_base(fenv$host_x)` creates `dgeMatrix` from raw data; `x@Dimnames` is not consulted in this deferred path — **names permanently lost** | **HIGH** |
| `R/qr.R:726` | `b <- as.matrix(b)` inside `am_solve` | Drops rownames from RHS before backsolve; solution matrix has no row names | **MED** |
| `R/wrappers.R:177` | `t(as.matrix(.amatrix_host_arg(value)))` | Transpose without explicitly swapping dimnames; rewrap afterward may not restore them | **MED** |
| `R/svd-factor.R:628–629` | `u <- as.matrix(svd_result$u); v <- as.matrix(svd_result$v)` | If u/v are adgeMatrix, the `as.matrix` dispatch should preserve names, but only if `amatrix_materialize_host` is invoked; direct `as.matrix` on a plain-matrix backend result loses them | **MED** |
| `R/chol-factor.R:39` | `mat <- as.matrix(factor@factor_obj)` | Factor obj dimnames silently dropped | **LOW** |
| `R/backend-cpu.R:17,58,66,78,92` | `as.matrix(x)` in backend dispatch | CPU backend discards adgeMatrix dimnames before dispatch | **LOW** |

---

## Grep 8 — `identical()` on float results

**Pattern:** `identical(.*[0-9]\.[0-9]` or `identical(.*as\.numeric`
**Total hits in R/:** 0
**Total hits in tests/:** 0

No floating-point `identical()` misuse found.

---

## Grep 9 — Missing `on.exit` cleanup

**Pattern:** `set_resident_context|acquire_lock|tempfile(` in R/*.R
**Total hits:** 0 (`set_resident_context` and `acquire_lock` not used; no `tempfile()` calls)

**`on.exit` sites found:** 14 total

| File:Line | What it guards | Assessment |
|-----------|----------------|------------|
| `R/chol-factor.R:303,338` | `rhs_arg` and `z_key` resident cleanup | Correct — uses `add = TRUE` |
| `R/sinkhorn.R:273` | resident temp cleanup | Correct |
| `R/models-lm.R:159` | deferred GPU option restore | Correct |
| `R/product-plan.R:118` | `amatrix.defer_host` option restore | Correct |
| `R/svd-factor.R:312–313` | product plan release | Correct |
| `R/policy.R:283` | policy restore | Correct |
| `R/irlba.R:155,631,634,699,702,803,804` | resident handle / Lanczos operator cleanup | Correct |

**Gap found:** `R/wrappers.R` has ~40 sites that allocate `out_key`/`scaled_key` via `_next_resident_key` inside a tryCatch block but do NOT use `on.exit` for cleanup — they rely on manual try-drop after the tryCatch. This is fragile: if a second error occurs between the tryCatch return and the manual drop (e.g., `.amatrix_rewrap_like` throws), the key leaks permanently.

Representative sites: wrappers.R:361, 402, 428, 463, 506, 534, 577, 627, 671, 701, 1156, 1283–1299, 1345–1359, 1415–1422, 1803, 2023, 2091–2098.

**Severity: MED** — GPU memory leak under error conditions.

---

## Grep 10 — Return type drift in S4 methods

**Pattern:** `setMethod` then last expression — does it return `matrix` instead of `adgeMatrix`?
**Total hits scanned:** ~120 setMethod calls

| File:Line | Method | Return type issue | Severity |
|-----------|--------|-------------------|----------|
| `R/methods-dense.R:67–69` | `%*%,numeric,adgeMatrix` | Returns result of `am_crossprod(y, matrix(x, ncol=1L))` — wraps numeric `x` as plain matrix, dispatches to am_crossprod which should return adgeMatrix. **OK** | OK |
| `R/methods-dense.R:147` | `coerce,aTransposeView,dgeMatrix` | `new_adgeMatrix(t(as.matrix(...)))` — correctly returns adgeMatrix | OK |
| `R/methods-dense.R:264–280` | `[,adgeMatrix` all variants | Delegates to `am_subset` which can return plain matrix or vector with `drop=TRUE` | **MED** — `drop=TRUE` path returns base `matrix`/`numeric`, breaking dispatch for callers expecting aMatrix |
| `R/methods-coercion.R:22` | `as.matrix,adgeMatrix` | Returns `as.matrix(dgeMatrix)` — returns plain `matrix`, **correct for coercion** | OK |
| `R/wrappers.R:978` | `am_transpose` for `adgCMatrix` | `t(as.matrix(amatrix_materialize_host(x)))` — returns plain `matrix`, then `.amatrix_rewrap_like` wraps it back. But if rewrap fails, returns raw matrix | **LOW** |
| `R/backend-cpu.R:58,66,78` | CPU backend matmul/crossprod | Returns `as.matrix(x %*% y)` — plain matrix, not adgeMatrix. This is correct since backend results are re-wrapped by the caller | OK |

---

## Top 10 Confirmed Bugs

| Rank | File:Line(s) | Pattern | Hypothesis | Severity |
|------|-------------|---------|------------|----------|
| 1 | `R/residency.R:396,403` | Dimnames loss in deferred path | GPU backend returns plain `matrix`; deferred materialization builds `dgeMatrix` from raw `fenv$host_x` without consulting `x@Dimnames`; **names permanently lost** | **HIGH** |
| 2 | `R/residency.R` (no FALSE site) | `host_cache_valid` never invalidated | In-place resident ops never set `host_cache_valid <- FALSE`; `amatrix_materialize_dense` returns stale host data after GPU mutation | **HIGH** |
| 3 | `R/chol-factor.R:175–191` | Double resident-drop on error | `solve_triangular_resident` error handler drops `out_key`; unconditional drop at lines 189–190 then fires again — double-drop corrupts backend key registry | **HIGH** |
| 4 | `R/wrappers.R:1283–1299` | Double-drop race in `crossprod_weighted` | `scaled_key` dropped inside tryCatch on success path; fallback at line 1297 tries to drop it again unconditionally | **HIGH** |
| 5 | `R/wrappers.R:1345–1359` | Double-drop race in `tcrossprod_weighted` | Same pattern as #4 | **HIGH** |
| 6 | `R/models-lm.R:579–589` | Resident key leak on broadcast error | `out_key` pre-allocated before tryCatch; error path can skip cleanup if `result` is NULL but fallback path also skips the drop | **HIGH** |
| 7 | `R/residency.R:400`, `R/resident-handle.R:109,170,199` | Unclassed critical errors | Unclassed `stop()` in core residency API; callers cannot catch GPU-specific failures by class | **HIGH** |
| 8 | `R/wrappers.R` ~40 sites | Missing `on.exit` for resident keys | Pre-allocated keys in tryCatch blocks with manual post-hoc cleanup; secondary error between tryCatch and manual drop leaks GPU memory permanently | **MED** |
| 9 | `R/qr.R:726` + `R/wrappers.R:177` | Dimnames lost through as.matrix | RHS dimnames dropped before backsolve; transpose rewrap does not restore swapped names | **MED** |
| 10 | `R/irlba.R:303,367,400` | C-level Lanczos errors swallowed | `.Call("amatrix_block_reorth_bridge", ...)` partial C-state mutation before error swallowed to NULL; Lanczos basis may be partially corrupted | **MED** |

---

## Summary Statistics

| Grep | Total Hits | High | Med | Low | Confirmed Bugs |
|------|-----------|------|-----|-----|----------------|
| 1. Unclassed stop() | ~90 | 10 | 20 | 60 | 7 unique callsites (high), ~30 med/low |
| 2. tryCatch → NULL | ~80 | 4 | 5 | 2 | 9 true bugs, rest intentional |
| 3/4. host_cache_valid | 3 | 1 | 0 | 0 | 1 (never-invalidated flag) |
| 5. Finalizer order | 3 | 0 | 0 | 0 | 0 |
| 6. <<- / slot mutation | 2 | 0 | 0 | 0 | 0 |
| 7. as.matrix dimnames | ~160 | 1 | 4 | 3 | 5 |
| 8. identical() float | 0 | 0 | 0 | 0 | 0 |
| 9. Missing on.exit | 14 ok / ~40 gaps | 0 | 1 | 0 | 1 (wrappers.R key leak pattern) |
| 10. S4 return drift | ~120 | 0 | 1 | 1 | 1 ([,drop=TRUE] returns base matrix) |
| **TOTAL** | | **~16** | **~30** | **~66** | **~24 unique bug sites** |
