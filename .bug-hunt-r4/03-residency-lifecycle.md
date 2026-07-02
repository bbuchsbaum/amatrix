# 03 Residency Lifecycle — Round 4 Bug Hunt

Generated: 2026-04-14  
Hunter: 03 (retry), round 4  
Package commit: eaf8c43

---

## (a) Drift Check

`packageVersion("amatrix")` → `0.1.0`, matches HEAD eaf8c43. OK.

---

## (b) Resident-Table Inspection API Used

- **Registry environment**: `amatrix:::.amatrix_state$residency`  
  Keys are `"obj:<object_id>"` strings. Count via:  
  `length(ls(envir = amatrix:::.amatrix_state$residency, all.names = TRUE))`
- **Per-object query**: `amatrix_residency_info(x)` (exported) — returns backend, resident_key, pinned_backend, live columns.
- **Internal helpers** (accessed via `:::`): `.amatrix_resident_entry(x)`, `.amatrix_bind_resident()`, `.amatrix_drop_resident_binding()`.
- **Note**: `am_crossprod`, `am_tcrossprod`, `am_chol`, `am_svd`, `am_solve` are present in the namespace but **not exported** in the installed NAMESPACE (stale vs HEAD). Required `:::` to reach them. Marked "lower reachability."

---

## (c) Per-Op Delta Table (CPU backend)

CPU backend is **not residency-capable** (`resident_store` absent). All ops produce zero deltas — no registry entries are ever created on this backend. MLX is registered but `available()` returns FALSE on this machine, so MLX live probes were not possible.

| op | happy-delta | error-delta | verdict | notes |
|----|-------------|-------------|---------|-------|
| `%*%` (matmul) | 0 | 0 | CLEAN | CPU, no residency |
| `am_crossprod` | 0 | 0 | CLEAN | CPU; `:::` required |
| `am_tcrossprod` | 0 | 0 | CLEAN | CPU; `:::` required |
| `am_chol` | 0 | 0 | CLEAN | CPU; `:::` required |
| `am_qr` | 0 | 0 | CLEAN | CPU |
| `am_svd` | 0 | 0 | CLEAN | CPU; `:::` required |
| `am_solve` | 0 | 0 | CLEAN | CPU; `:::` required |
| `am_sweep` | 0 | 0 | CLEAN | CPU |

**Limitation**: MLX unavailable at runtime — GPU-path deltas are untestable on this machine.  
**[STAT:n]** n=8 ops probed, n=0 leaks confirmed on CPU.

---

## (d) Confirmed Leaks with Repros

No registry-level (double-entry) leaks were confirmed. See amatrix-3ka below for a GPU-buffer-level leak in the internal path.

---

## (e) Refutation Verdicts

### amatrix-5ni — `am_sweep_inplace` leaks `new_key` on backend error

**Code path** (`resident-handle.R:258–276`):

```
old_key <- h$resident_key
new_key <- .amatrix_next_resident_key(h$backend_name)
tryCatch(
  backend$broadcast_ewise_resident(..., new_key, defer = TRUE),
  error = function(e) {
    err <<- e
    .rh_drop_key(backend, new_key)   # <-- cleanup IS present
  }
)
if (!is.null(err)) stop(err)
```

`.rh_drop_key` (lines 120-126) calls `backend$resident_drop(new_key)` if `backend$resident_has(new_key)` is TRUE. The drop is guarded on `resident_has`, so if the backend allocated `new_key` before erroring, it gets dropped; if the error fired before any allocation, there is nothing to drop.

**Verdict: REFUTED at the code level.** The error handler calls `.rh_drop_key(backend, new_key)` immediately in the `tryCatch` error branch, before re-raising. No leak path exists unless `broadcast_ewise_resident` both allocates the key AND fails without the allocation appearing in `resident_has` — a backend-internal inconsistency, not an amatrix bug.  
Live verification impossible (MLX unavailable), so this is a **code-level refutation only**.

---

### amatrix-3ka — `.amatrix_bind_resident` leaks prior key on rebind

**Code path** (`residency.R:90–116`):

```r
.amatrix_bind_resident <- function(x, backend, resident_key, sparse = FALSE) {
  object_key <- .amatrix_object_key(x)
  # ...
  assign(object_key, list(backend=backend, resident_key=resident_key, ...), 
         envir = .amatrix_state$residency)
  # No drop of old GPU buffer before overwrite
}
```

The internal function **only calls `assign()`** — it overwrites the registry entry with the new key but does **not** call `backend$resident_drop(old_resident_key)` on the displaced entry.

**Reproduction** (confirmed in this session via direct state manipulation):

```r
# Inject a prior binding for key1
assign("obj:<id>", list(backend="mlx", resident_key="key1", ...), 
       envir = .amatrix_state$residency)
# Rebind with key2
.amatrix_bind_resident(A, "mlx", "key2")
# Registry entry now shows key2 — but key1 is NEVER dropped
# dropped_keys remains empty
```

**Registry count stays at 1** (no double-entry), but the **old GPU buffer at `key1` is silently abandoned**.

**Verdict: CONFIRMED — GPU-buffer leak in `.amatrix_bind_resident` (internal).** The public `amatrix_bind_resident` (bind-resident.R:76-86) correctly drops the old key before rebinding, so the **public API is safe**. The leak only manifests when internal callers invoke `.amatrix_bind_resident` directly on an object that already has a different resident_key — e.g., `wrappers.R:170`, `wrappers.R:347`, `wrappers.R:1042`, etc. Any internal op that produces a fresh GPU result and binds it to an object that was previously bound to a *different* key will silently leak the old GPU buffer.

---

## (f) Lint/CI Proposal

Add a testthat test using a mock residency backend that instruments `resident_drop` calls; assert that every internal call to `.amatrix_bind_resident` on an already-bound object first results in one `resident_drop` call for the displaced key (or add an explicit drop-before-overwrite guard inside `.amatrix_bind_resident` itself).

```r
# In .amatrix_bind_resident, after computing object_key:
existing <- get0(object_key, envir = .amatrix_state$residency, inherits = FALSE)
if (!is.null(existing) && !identical(existing$resident_key, resident_key)) {
  # drop old GPU buffer before overwriting
  old_be <- tryCatch(.amatrix_get_backend(existing$backend), error = function(e) NULL)
  if (!is.null(old_be) && is.function(old_be$resident_drop))
    try(old_be$resident_drop(existing$resident_key), silent = TRUE)
}
```

