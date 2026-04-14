# Round 2 Bug Hunt: Resident-Handle State Coherence

## Scope

Traces every mutation path in the amatrix R package and verifies that
derived state (host_cache_valid, factorization caches, residency bindings,
finalizers, backend_state) is invalidated or preserved correctly.

---

## Invariant I1: After any in-place numerical mutation, all cached factorizations on the object MUST be cleared.

### Mutations that touch numerical content

| Operation | Entry point | Backed by |
|-----------|-------------|-----------|
| `[<-` (element assign) | `am_subassign` (wrappers.R:986) | materialises host, modifies, calls `.amatrix_rewrap_value` |
| `cbind2`/`rbind2` | `.amatrix_bind2` (wrappers.R:72) | materialises host, calls `.amatrix_rewrap_value` |
| `am_set_dimnames` | wrappers.R:2366 | materialises host, calls `.amatrix_rewrap_like` |
| Arithmetic Ops (`+`, `-`, `*`, `/`) | `am_arith` (wrappers.R) | result rewrapped via `.amatrix_rewrap_value` |

### Cache key construction

- `chol_factor` (chol-factor.R:130): key = `paste0("chol:", X@object_id)`
- `am_svd` (svd-factor.R:472): key = `"svd:" + X@object_id + ...`
- `lm_fit` / QR path (models-lm.R:186): key includes `X@object_id`

### The bug: `.amatrix_rewrap_like` always allocates a fresh `object_id`

`am_subassign` calls `.amatrix_rewrap_value` → `.amatrix_rewrap_like` → `new_adgeMatrix`
→ `.amatrix_new_dense` (constructors.R:41): `object_id <- .amatrix_next_object_id()`.

Every mutation that goes through `.amatrix_rewrap_like` produces a **new object** with a
**new object_id**. Because the cache is keyed on `object_id`, the old cached
factorization is never retrieved for the new object. This means the cache does NOT serve
stale data after a mutation — the stale entry simply leaks and is never freed until LRU
eviction.

**I1 verdict: technically no stale-read, but there IS a cache-leak violation.**
Old cache entries accumulate in `.amatrix_state$model_cache` for object_ids that no longer
correspond to any live R object. No code path calls `rm(list = "chol:<old_id>")` or any
equivalent invalidation. With default `cache_max_size = Inf`, these entries never evict.

**Violation sites:**
- R/chol-factor.R:154 — `chol:X@object_id` written; old key from prior object never removed
- R/svd-factor.R:654 — same pattern for SVD
- R/models-lm.R:300 — same pattern for QR/XtX

---

## Invariant I2: host_cache_valid coherence after device→host and host mutation

### Where host_cache_valid is SET to TRUE

`R/residency.R:97`:
```r
if (inherits(x, "adgeMatrix") && !isTRUE(sparse) && !isTRUE(x@finalizer_env$host_deferred)) {
  x@finalizer_env$host_cache_valid <- TRUE
}
```
This is called from `.amatrix_bind_resident`. It sets `host_cache_valid = TRUE` on the
object at bind time, meaning "the host @x slot is valid".

### Where host_cache_valid is NEVER reset to FALSE

Grep confirms there is no site in any R file that sets `host_cache_valid <- FALSE` or
clears it. The flag can only be set to TRUE (residency.R:97) or read (residency.R:406).

### The bug: host mutation after GPU push never clears host_cache_valid

Scenario:
1. `A <- adgeMatrix(m)` — object created, no residency, host_cache_valid = NULL (falsy)
2. `amatrix_push_gpu(A)` or an operation that calls `.amatrix_bind_resident` — sets
   `host_cache_valid = TRUE` in `A@finalizer_env`
3. At this point `am_subassign` performs a HOST mutation: `am_subassign` materialises the
   host copy, modifies it, and calls `.amatrix_rewrap_value` → `new_adgeMatrix`. The
   result is a **new** object so `host_cache_valid` on the new object starts as NULL (no
   harm).

However the symmetric danger is: a GPU-resident op writes a new value to device memory
while the `adgeMatrix` R object still has `host_cache_valid = TRUE` from the bind. If any
code path updates the device-side buffer WITHOUT creating a new object and WITHOUT clearing
`host_cache_valid`, subsequent `amatrix_materialize_dense` at residency.R:406 will return
the stale host `@x` instead of downloading the fresh device value.

**Confirmed path:** `as_adgeMatrix.resident_handle` (resident-handle.R:393–417):
- Creates a new `adgeMatrix` (fresh object_id, fresh `host_cache_valid = NULL`)
- Calls `.amatrix_bind_resident` → sets `host_cache_valid = TRUE` on the new object
- This is correct IF `defer_host = FALSE` because host data was just downloaded.
- When `defer_host = TRUE` (line 395–401): a deferred object is created. Its
  `host_deferred = TRUE`. Then `.amatrix_bind_resident` is called (line 412). In
  `.amatrix_bind_resident` (residency.R:96–97) the guard is:
  `!isTRUE(x@finalizer_env$host_deferred)` — so `host_cache_valid` is NOT set to TRUE for
  deferred objects. That part is correct.

**Actual I2 violation — residency.R:97 unconditional set-on-bind:**

`.amatrix_bind_resident` is called in scenarios beyond construction:
- After a backend operation that stores a new result key (`wrappers.R` resident paths)
- In those calls the object's `@x` slot still holds the OLD host data.
- Setting `host_cache_valid = TRUE` at bind time is premature: the host `@x` may be
  stale relative to the freshly stored resident key.

Concretely, any resident-path operation that:
1. Materialises the input, dispatches to GPU, stores result under a new resident key, and
2. calls `.amatrix_bind_resident` on an object whose `@x` still holds the pre-op data

will leave `host_cache_valid = TRUE` on the output object even though `@x` is stale.
The next `amatrix_materialize_dense` call at residency.R:406 returns the wrong matrix.

**Violation site:** R/residency.R:97

---

## Invariant I3: After backend switch, old resident handle MUST be released exactly once.

### Release mechanism

`.amatrix_release_resident` (residency.R:117) calls `backend$resident_drop` and then
`.amatrix_drop_resident_binding`. This is correct IF called.

`.amatrix_bind_resident` (residency.R:84) **overwrites** the residency registry entry
without first calling `_release_resident`. If an object already has a resident key and
`.amatrix_bind_resident` is called again (e.g., after a second GPU push or during
`as_adgeMatrix.resident_handle`), the old resident key is **silently replaced** in the
registry env without dropping the device buffer.

**The finalizer** (residency.R:9–52) uses the key stored at finalizer time — which is the
LAST key written. The first (overwritten) key is never freed.

**Violation site:** R/residency.R:84–100 — `.amatrix_bind_resident` has no guard to
release an existing entry before overwriting.

**Reproducer sketch:**
```r
A <- adgeMatrix(matrix(1:4, 2, 2), preferred_backend = "mlx")
# push to GPU — stores key K1
.amatrix_bind_resident(A, "mlx", "mlx:1")
# push again (e.g., after an op that re-uploads)
.amatrix_bind_resident(A, "mlx", "mlx:2")
# K1 is never dropped — device memory leak
```

---

## Invariant I4: After clone/copy, finalizers must be set up for BOTH copies independently.

### `.amatrix_new_dense` (constructors.R:25–55)

Every call to `new_adgeMatrix` / `.amatrix_new_dense` generates a fresh `object_id` and
fresh `finalizer_env`. This is correct.

### `as_adgeMatrix.resident_handle` (resident-handle.R:393)

When `defer_host = FALSE`:
- `new_adgeMatrix(mat, ...)` — fresh object_id, fresh finalizer_env. OK.
- `.amatrix_bind_resident(obj, h$backend_name, h$resident_key)` — binds under new id.
- `h$active <- FALSE; h$resident_key <- NULL` — handle's finalizer will not drop the key. OK.

**I4 verdict: no violation found here.** The transfer-ownership pattern is implemented
correctly for `resident_handle` → `adgeMatrix`.

---

## Invariant I5: After dimnames change, no cached dimension-dependent state may persist.

### `am_set_dimnames` path (wrappers.R:2366–2369)

```r
am_set_dimnames <- function(x, value) {
  host_x <- amatrix_materialize_host(x)
  dimnames(host_x) <- value
  .amatrix_rewrap_like(x, host_x)
}
```

`.amatrix_rewrap_like` → `new_adgeMatrix` → fresh `object_id`. The old cache entries
(chol, SVD, QR) are keyed on the OLD `object_id`. After a dimnames change:
- The new object has a new `object_id`, so it will not hit the old cached factorization.
- The old cache entry leaks (same I1 leak).

**I5 verdict: no stale-read for factorizations.** Dimnames are structurally not part of
the cache key for chol/svd/QR, so no correctness hazard. The leak issue from I1 applies.

---

## Consolidated Violation Summary

| # | Invariant | Severity | File:line | Description |
|---|-----------|----------|-----------|-------------|
| V1 | I2 | HIGH | R/residency.R:97 | `.amatrix_bind_resident` sets `host_cache_valid = TRUE` unconditionally for any non-deferred dense bind, even when the object's `@x` is stale (pre-op host data). Causes `amatrix_materialize_dense` to return wrong matrix via residency.R:406 fast-path. |
| V2 | I3 | MEDIUM | R/residency.R:84–100 | `.amatrix_bind_resident` overwrites registry entry without releasing old resident key first. Device buffer K1 is leaked when a second bind replaces K1 with K2. |
| V3 | I1/I5 | LOW | R/chol-factor.R:154, R/svd-factor.R:654, R/models-lm.R:300 | No cache invalidation function exists. Old entries accumulate under dead `object_id` keys. Default `cache_max_size = Inf` means they never evict. Not a stale-read bug (new objects get new keys) but a memory leak for long-running sessions. |

---

## Reproducer Sketches

### Sketch 1 — V1: stale host_cache_valid after GPU op (I2 violation)

```r
library(amatrix)
m <- matrix(1:4, 2, 2) * 1.0
A <- adgeMatrix(m, preferred_backend = "mlx")

# Suppose a GPU matmul stores result under key K and calls .amatrix_bind_resident.
# Internally this happens in any resident-path wrapper that returns a new adgeMatrix
# whose @x was filled with the pre-op host data but host_cache_valid is set TRUE.
# Trigger: any resident path that goes through .amatrix_bind_resident on an object
# whose @x is NaN sentinel (deferred) — but host_deferred flag was NOT set.

# Manual repro of the hazard:
# 1. Build an object with old @x
old_x <- c(1, 0, 0, 1)  # identity
A2 <- new("adgeMatrix", x = old_x, Dim = 2L:2L, Dimnames = list(NULL,NULL),
          factors = list(), preferred_backend = "cpu", policy = "eager",
          precision = "strict",
          object_id = "test:1", src_id = "",
          finalizer_env = amatrix:::.amatrix_make_finalizer_env("test:1"))
# 2. Bind a resident key that corresponds to a DIFFERENT matrix (2x identity scaled by 2)
.amatrix_bind_resident(A2, "cpu", "cpu:999")
# 3. host_cache_valid is now TRUE; @x still holds old identity matrix
# 4. amatrix_materialize_dense returns identity, not the 2x version
result <- amatrix_materialize_dense(A2)
# Expected: 2*identity; Actual: identity — WRONG
```

### Sketch 2 — V2: double-bind leaks device buffer (I3 violation)

```r
library(amatrix)
# Requires a backend with residency (mlx or arrayfire)
m <- matrix(runif(4), 2, 2)
A <- adgeMatrix(m, preferred_backend = "mlx")

# First resident push — K1 allocated
key1 <- .amatrix_next_resident_key("mlx")
bk <- .amatrix_get_backend("mlx")
bk$resident_store(key1, m)
A <- .amatrix_bind_resident(A, "mlx", key1)

# Second resident push (simulate result of another op)
key2 <- .amatrix_next_resident_key("mlx")
bk$resident_store(key2, m * 2)
A <- .amatrix_bind_resident(A, "mlx", key2)

# Now key1 is orphaned: still alive in device memory, never dropped
stopifnot(isTRUE(bk$resident_has(key1)))  # should be FALSE after proper release
```

### Sketch 3 — V3: factorization cache leak after mutation (I1/I5)

```r
library(amatrix)
m <- crossprod(matrix(rnorm(16), 4)) + diag(4)
A <- adgeMatrix(m)
fac1 <- chol_factor(A)  # cached under "chol:<A@object_id>"

# Mutate A — produces new object with new object_id
A[1, 1] <- A[1, 1] + 0.5
# Old cache entry "chol:<old_id>" is now orphaned in .amatrix_state$model_cache
# It will never be freed unless LRU eviction fires or session restarts.
# In a LOO-CV loop with n=1000, this creates 1000 leaked chol entries.
n_cached <- length(ls(amatrix:::.amatrix_state$model_cache))
stopifnot(n_cached <= 1)  # fails if leak is present
```

---

## Severity Assessment

- **V1 (HIGH)**: Can silently return wrong numerical results from `amatrix_materialize_dense`
  in any code path where `.amatrix_bind_resident` is called on an object whose `@x` is
  not current. This is the highest-risk correctness bug.

- **V2 (MEDIUM)**: Device memory leak. In iterative algorithms with many GPU ops, this
  accumulates stale resident keys that are never freed. Does not corrupt values but can
  exhaust device memory.

- **V3 (LOW)**: Host memory leak only. Values are correct (new object_id = cache miss =
  recompute). However in long sessions or LOO-CV loops the leak is O(n_mutations).
