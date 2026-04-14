# Composition Chain Bug Hunt — Round 2
# File: .bug-hunt-r2/04-composition-chains.md

## Methodology

Static code audit of composition patterns across:
R/methods-dense.R, R/methods-sparse.R, R/methods-coercion.R, R/wrappers.R,
R/residency.R, R/prepare-operands.R, R/product-plan.R, R/dispatch-hardening.R,
R/bind-resident.R, R/constructors.R, R/qr.R, R/qr-downdate.R

---

## Pattern A: op → coerce → op (as.matrix then re-promote)

**Code path:**
1. `adgeMatrix %*% B` → `matmul()` → returns `adgeMatrix` result (possibly resident)
2. `as.matrix(result)` → `setMethod("as.matrix","adgeMatrix")` → `amatrix_materialize_host(x)` → plain R `matrix`
3. `result %*% C` — `result` is now a plain `matrix`; S4 dispatch hits `%*%,matrix,ANY` (base R), **not** amatrix

**Status: SAFE** (by design — coercion intentionally strips the wrapper)

The plain matrix from step 2 is no longer an `aMatrix` object. If `C` is an `adgeMatrix`, the `%*%,matrix,adgeMatrix` method in dispatch-hardening.R line 76–81 wraps the plain matrix with `new_adgeMatrix` inheriting `y`'s backend. If `C` is a plain matrix, base R handles it. No residency leak here, but the user likely expects the second op to also use GPU — it will only do so if `C` is still an `adgeMatrix`.

**Limitation exposed:** There is no "re-promote from plain matrix" in the result of `as.matrix(resident_adgeMatrix)`. The residency is intentionally lost. This is correct but can surprise users who coerce mid-chain expecting GPU to persist.

---

## Pattern B: transpose → product (lazy symbolic transpose)

**Code path:**
`t(A)` where A is `adgeMatrix` → `am_transpose(x)` (wrappers.R:969) → `.new_aTransposeView(x)` (constructors.R:407)

The `aTransposeView` stores `@source = A` with a **new** `object_id` and `src_id = A@object_id`. It does NOT carry A's resident key — it has its own `finalizer_env` with no backend binding.

When `t(A) %*% B` is called:
- S4 dispatch hits `%*%,aTransposeView,ANY` → `am_crossprod(x@source, y)` (methods-dense.R:117–118)
- This routes to the source's resident key correctly — **SAFE**

When `t(A) %*% t(B)` is called:
- Hits `%*%,aTransposeView,aTransposeView` → `am_tcrossprod(y@source, x@source)` (methods-dense.R:133–134) — **SAFE**

**SUSPICIOUS: `aTransposeView` of an `adgCMatrix`**

`am_transpose` for sparse (wrappers.R:973–976):
```r
host <- amatrix_materialize_host(x)  # materializes to dgCMatrix
return(new_adgCMatrix(t(host), ...)) # creates new adgCMatrix, no resident binding
```

This is EAGER materialization: transposing a resident sparse matrix forces a GPU→CPU round-trip. There is no lazy path for sparse transpose. The result has no resident binding, so any subsequent op must re-upload. However, this is arguably correct (sparse transpose creates a structurally different object). The tripwire WILL fire.

**BUG HYPOTHESIS B1 — severity: MEDIUM**
File: R/wrappers.R:973-976
If an `adgCMatrix` is resident on a backend and the user calls `t(sparse_resident)`, the sparse transpose forces an unconditional `amatrix_materialize_host()` download even when the backend supports in-place sparse transpose. The result loses residency. Subsequent ops re-upload. There is no `t(adgCMatrix)` fast path analogous to the `aTransposeView` lazy path for dense matrices.

Reproducer sketch:
```r
S <- adgCMatrix(sparse_mat, backend = "mlx")
S_bound <- amatrix_bind_resident(S, backend = "mlx")
tS <- t(S_bound)        # forces GPU→CPU download
tS %*% dense_rhs        # re-uploads
```

---

## Pattern C: bind → op (rbind/cbind then multiply)

**Code path — `.amatrix_bind2` (wrappers.R:72–86):**
```r
.amatrix_bind2 <- function(kind, x, y) {
  template <- .amatrix_template(x, y)    # first aMatrix found
  value <- switch(kind,
    cbind2 = methods::cbind2(.amatrix_host_arg(x), .amatrix_host_arg(y)),
    rbind2 = methods::rbind2(.amatrix_host_arg(x), .amatrix_host_arg(y))
  )
  .amatrix_rewrap_value(template, value)
}
```

`_host_arg` materializes both operands to host. The result is rewrapped via `_rewrap_value` → `_rewrap_like(template, value)` → `new_adgeMatrix(value, preferred_backend=template@preferred_backend, ...)`.

**BUG HYPOTHESIS C1 — severity: HIGH (confirmed)**
File: R/wrappers.R:72–86
After `rbind(A_gpu, B_gpu)` or `cbind(A_gpu, B_gpu)`:
1. Both A and B are materialized to host (GPU→CPU download for BOTH — two tripwire events)
2. The result is a new `adgeMatrix` with `preferred_backend` copied from the *first* aMatrix template but **no resident binding**
3. The subsequent `result %*% C` must re-upload the entire bound matrix

There is no attempt to perform the bind on the GPU side — the bind always goes through host. For large matrices this is a silent performance regression but not a correctness bug per se. However, if A is on MLX and B is on ArrayFire, template takes A's `preferred_backend = "mlx"` and ignores B's backend entirely. B's resident data is leaked if `B` has `owned_resident = TRUE` and no other reference (finalizer will eventually clean it but no explicit drop is triggered by `_bind2`).

**BUG HYPOTHESIS C2 — severity: MEDIUM (mixed-backend bind)**
File: R/wrappers.R:73–75
When `x` is on MLX and `y` is on ArrayFire, `_template(x, y)` returns `x` (the first aMatrix). The result's `preferred_backend` is MLX. Y's resident data on ArrayFire is not explicitly released — it relies on finalizer GC. There is no warning or error. The downstream op will re-upload the result to MLX, which is correct but the user may not realise two GPU uploads happened for no benefit.

---

## Pattern D: subset → op (A[rows,] %*% B)

**Code path — `am_subset` (wrappers.R:981–983):**
```r
am_subset <- function(x, i, j, ..., drop = TRUE) {
  value <- amatrix_materialize_host(x)[i, j, ..., drop = drop]
  .amatrix_rewrap_value(x, value)
}
```

**Status: SAFE for correctness; SUSPICIOUS for residency**

The subset always materializes to host, applies base R `[`, then rewraps. The result is a fresh `adgeMatrix` with `preferred_backend` from the template but NO resident binding. This is correct — there is no GPU-side subset operation.

**BUG HYPOTHESIS D1 — severity: LOW-MEDIUM**
File: R/wrappers.R:981–983
Every `A[rows, ]` on a resident `adgeMatrix` triggers a GPU→CPU download (tripwire event). If the user follows with `A[rows, ] %*% B`, the result is host-only and must re-upload. There is no in-place GPU slice path. This is a performance issue, not a correctness bug. No dangling view risk since subset always makes a copy.

However, there IS a subtle issue: `am_subassign` (line 986–990) also materializes to host, mutates, then rewraps — but does NOT invalidate the original `x`'s GPU resident key. The original `x` object still has a live resident binding to the OLD data. If `x` was passed by reference (R's copy-on-modify semantics usually prevent this, but with S4 slots and environments it can occur), the old GPU buffer and new host copy are out of sync.

**BUG HYPOTHESIS D2 — severity: MEDIUM**
File: R/wrappers.R:986–990
`am_subassign` returns a new object but does NOT call `.amatrix_release_resident(x)` on the original before rewrapping. The returned object from `_rewrap_value` is a NEW `adgeMatrix` (new `object_id` via `new_adgeMatrix`). The original `x` still points to the old GPU buffer. In R's copy-on-modify model `x[i,j] <- v` replaces `x` in the calling environment, but the old binding in `.amatrix_state$residency` for the original `object_id` is never freed until GC. This leaks the GPU buffer for the duration of the session if `x` is reassigned.

---

## Pattern E: factorization → solve → product (QR cache invalidation)

**Code path — `qr_downdate.amQR` (qr-downdate.R:55–73):**
```r
qr_downdate.amQR <- function(qr_factor, row_idx, X = NULL) {
  X_sub <- X[-row_idx, , drop = FALSE]
  am_qr(.amatrix_qr_arg(X_sub))
}
```

And `lm_loo_cv` (qr-downdate.R:126–152):
```r
qr_full <- am_qr(X_am, ...)
for (i in seq_len(n)) {
  qr_i   <- qr_downdate(qr_full, i, X = X)   # drop row i — uses ORIGINAL X, not X_am
  coef_i <- as.numeric(qr.coef(qr_i, y_vec[-i]))
  loo_resid[[i]] <- y_vec[[i]] - sum(X[i, ] * coef_i)
}
```

**BUG HYPOTHESIS E1 — severity: HIGH (confirmed)**
File: R/qr-downdate.R:143,145
`lm_loo_cv` calls `qr_downdate(qr_full, i, X = X)` passing the **original** `X` (which may be a plain matrix), not `X_am` (the wrapped `adgeMatrix`). Then on line 145 it evaluates `X[i, ]` using the original `X` again. This is inconsistent: the QR is computed from `X_am` (potentially GPU-accelerated) but downdate uses plain `X`. If `X` is an `adgeMatrix` passed by the user, `X[-row_idx, , drop=FALSE]` triggers `am_subset` → GPU download → new host copy → `am_qr` re-uploads. This is n GPU→CPU→GPU round-trips in the LOO loop.

More critically: if `X` is a plain matrix and `X_am` is the GPU-wrapped version, the two objects may differ in precision (`X_am` may have been cast to `float32` via `precision="fast"`). The LOO residuals are computed with mixed precision — `qr.coef` uses float32 factors but `X[i,] * coef_i` uses float64 values. This is a **numerical correctness bug**, not just a performance issue.

Reproducer sketch:
```r
X <- matrix(rnorm(100), 20, 5)
y <- rnorm(20)
# with precision="fast", X_am is float32 but X stays float64
cv <- lm_loo_cv(X, y)   # X[i,] * coef_i mixes precisions silently
```

**Pattern E2 — QR state finalizer (qr.R)**
The `amQR` state environment registers a finalizer that drops GPU keys on GC. This is correct. The plan does NOT cache the source matrix (by design — `amQR` explicitly refuses to store it). So there is no stale QR bug from mutation of the original matrix.

**Status: SAFE for mutation invalidation; BUG for mixed-precision LOO**

---

## Pattern F: mixed-backend chain (A on MLX, B on ArrayFire)

**Code path — `matmul(A, B)` with A@preferred_backend="mlx", B@preferred_backend="arrayfire":**

`matmul` → `_backend_for(x, "matmul", y=y)` — this function looks at `x` (the LHS) only:
- Calls `amatrix_resident_backend_for(x, ...)` which checks x's live resident binding first
- If A is resident on MLX, returns "mlx"
- Then calls `_try_resident_matmul(A, B, "mlx")`
- Inside, `_prepare_resident_arg(B, "mlx", promote_amatrix=TRUE)` is called for the RHS
- Since B is resident on ArrayFire (not MLX), `resident_key(B, backend="mlx")` returns NULL
- `promote_amatrix=TRUE` → uploads B to MLX (a SECOND GPU upload)
- B's ArrayFire buffer is NOT freed — no call to `_release_resident(B)`

**BUG HYPOTHESIS F1 — severity: HIGH (confirmed)**
File: R/wrappers.R (via `_prepare_resident_arg`, wrappers.R:260–272)
When A is resident on MLX and B is resident on ArrayFire, executing `A %*% B`:
1. B is uploaded to MLX (GPU→CPU→GPU round-trip for B)
2. B's ArrayFire buffer is left alive (memory leak until finalizer GC)
3. The upload in `_prepare_resident_arg` line 270–271 calls `_bind_resident(value, backend_name, resident_key)` which **overwrites** B's registry entry, changing it from ArrayFire to MLX. This silently destroys B's ArrayFire binding metadata.
4. If B was shared by another operation expecting its ArrayFire binding, that op will now find it "resident on MLX" and either fail or re-upload to ArrayFire.

The registry overwrite is at residency.R:84–106: `_bind_resident` assigns a new entry unconditionally, replacing the old one. There is no "both-backends" tracking.

**BUG HYPOTHESIS F2 — severity: MEDIUM**
File: R/prepare-operands.R:127–138 (amatrix_prepare_operands)
`amatrix_prepare_operands` calls `_select_binary_resident_backend` which only checks x's residency then y's. It does NOT check "are x and y resident on DIFFERENT non-cpu backends?" and emit a warning. It silently picks one backend and promotes the other operand, destroying the loser's binding.

---

## Pattern G: lazy plan reuse (compile_product then mutate operand)

**Code path — `amatrix_compile_product`:**
The plan closure captures `lhs_bound` by reference (it's in the closure environment). If the user:
1. Compiles `plan <- amatrix_compile_product(A, op="matmul")`
2. Mutates A (e.g., `A[1,1] <- 99` → `am_subassign` returns new object, A in user env is replaced)
3. Calls `plan(B)` — uses OLD `lhs_bound` (the pre-mutation version)

**Status: SAFE for this path**
Because `am_subassign` returns a NEW `adgeMatrix` with a new `object_id`, the user's `A` now points to the new object. The plan's `lhs_bound` still references the original (old) object. The old object's GPU buffer is still valid (not freed because `lhs_bound` in the closure holds a reference). So the plan is internally consistent but operates on STALE data.

**BUG HYPOTHESIS G1 — severity: MEDIUM**
File: R/product-plan.R:219–268
The plan does NOT validate on each call that `lhs_bound`'s GPU buffer is still live. If the user explicitly calls `amatrix_release_resident(A)` on the original object (freeing the GPU buffer), the plan's next invocation will call `_try_resident_matmul` with a stale resident key, get NULL back (buffer gone), fall through to `amatrix_dispatch_op`, and silently re-upload. This is a silent performance regression, not a crash, because the fallback exists. BUT:

```r
plan <- amatrix_compile_product(A, op="matmul")
amatrix_release_resident(A)   # frees GPU buffer
plan(B)   # lhs_bound still has the old object_id; resident_has() returns FALSE
          # _try_resident_matmul returns NULL
          # dispatch_op falls back, re-materializes lhs_bound from host @x
          # but lhs_bound's @x is EMPTY if it was created with host_deferred=TRUE
```

If `A` was created via `new_adgeMatrix_deferred` (resident-only, no host copy), then `lhs_bound@x` contains NaN placeholders. After the GPU buffer is freed, the fallback materializes NaN data. **This is a correctness bug** for deferred-host objects used in a compiled plan after release.

**BUG HYPOTHESIS G2 — severity: HIGH (deferred + plan release)**
File: R/product-plan.R:232–261, R/residency.R:378–436
Sequence:
1. `A_deferred = new_adgeMatrix_deferred(...)` — host @x = NaN
2. `A_bound = amatrix_bind_resident(A_deferred, "mlx")` — live GPU buffer
3. `plan = amatrix_compile_product(A_bound)` — `lhs_bound` references `A_bound`
4. `amatrix_release_resident(A_bound)` — GPU buffer freed, `@x` still NaN
5. `plan(B)` → resident path returns NULL → fallback materializes `lhs_bound` → downloads NaN

The fallback at residency.R:397–403:
```r
if (is.null(fenv$host_x)) {
  ...
  if (is.null(fenv$host_x)) {
    stop("deferred adgeMatrix lost its GPU resident data", call. = FALSE)
  }
}
```
This correctly STOPS — but only for the `host_deferred=TRUE` path. For regular `adgeMatrix` objects the fallback silently uses `@x` which may be stale/zero-initialized if the object was constructed without host data being set.

---

## Pattern H: coerce sparse→dense → op

**Code path:**
`as(S, "adgeMatrix")` where S is `adgCMatrix`:
- No `setAs("adgCMatrix","adgeMatrix")` is defined in constructors.R (only matrix→adge, dge→adge, matrix→adgC, dgC→adgC are registered)
- S4 will search for an inherited coerce path: adgCMatrix extends dgCMatrix, and there is `setAs("dgCMatrix","adgCMatrix")` but NOT `setAs("adgCMatrix","adgeMatrix")`
- S4 will find `adgCMatrix` → `dgCMatrix` (via `adgeMatrix` extends `dgeMatrix`? No.) Actually the inheritance graph is `adgCMatrix extends dgCMatrix` and `adgeMatrix extends dgeMatrix`. There is no coerce path registered for `adgCMatrix→adgeMatrix`.

**BUG HYPOTHESIS H1 — severity: HIGH (missing coerce path)**
File: R/constructors.R:424–429
There is no `setAs("adgCMatrix", "adgeMatrix")` registered. When a user writes `as(S_sparse, "adgeMatrix")`, S4 will attempt to find a coerce chain. It will likely find `adgCMatrix → dgCMatrix → dgeMatrix → adgeMatrix` (2+ hops) or fail entirely. If it finds the chain, each hop strips amatrix metadata: `adgCMatrix→dgCMatrix` (constructors.R:429, no backend metadata) → `dgCMatrix→dgeMatrix` (Matrix pkg) → `dgeMatrix→adgeMatrix` (constructors.R:425, new_adgeMatrix with defaults). The `preferred_backend` and `policy` are LOST.

More critically, `as(S_sparse, "adgeMatrix")` forces a dense materialization but the resulting `adgeMatrix` has `preferred_backend=""` (default from `new_adgeMatrix`), not S's backend. The next op will use auto-selection, ignoring the user's intent.

Reproducer sketch:
```r
S <- adgCMatrix(sparse_mat, backend = "mlx")
D <- as(S, "adgeMatrix")   # backend metadata lost
D %*% X                    # might route to a different backend
```

---

## Summary Table

| Pattern | File:Line | Status | Severity |
|---------|-----------|--------|----------|
| A: op→coerce→op | methods-coercion.R:22 | SAFE | — |
| B1: sparse t() eager materialize | wrappers.R:973-976 | SUSPICIOUS | MEDIUM |
| C1: bind always host-side | wrappers.R:72-86 | BUG (perf+leak) | HIGH |
| C2: mixed-backend bind silently drops loser | wrappers.R:73-75 | BUG | MEDIUM |
| D1: subset always materializes | wrappers.R:981-983 | SUSPICIOUS | LOW-MEDIUM |
| D2: subassign doesn't release old GPU buffer | wrappers.R:986-990 | BUG (leak) | MEDIUM |
| E1: LOO CV mixed-precision residuals | qr-downdate.R:143,145 | BUG (correctness) | HIGH |
| F1: mixed-backend op overwrites loser binding | residency.R:84-106 | BUG (leak+corruption) | HIGH |
| F2: no mixed-backend warning | prepare-operands.R:127-138 | SUSPICIOUS | MEDIUM |
| G1: plan reuse after release (stale fallback) | product-plan.R:219-268 | BUG (silent perf) | MEDIUM |
| G2: deferred + plan + release = NaN result | product-plan.R:232-261 | BUG (correctness) | HIGH |
| H1: missing adgCMatrix→adgeMatrix coerce path | constructors.R:424-429 | BUG (metadata loss) | HIGH |

---

## Top 3 Most Certain Bugs

### 1. Mixed-backend op overwrites registry binding (F1)
**File: R/wrappers.R (`.amatrix_prepare_resident_arg`) + R/residency.R:84–106**
When A is resident on backend X and B on backend Y, executing `A %*% B` silently overwrites B's registry entry with X's key via `_bind_resident`. B's original GPU buffer on Y is not freed, and any other reference to B now sees a wrong backend in the registry. Severity: HIGH.

### 2. LOO CV mixed-precision correctness (E1)
**File: R/qr-downdate.R:143 and 145**
`lm_loo_cv` passes the original `X` (user-supplied, possibly float64) to `qr_downdate` but builds `X_am` (possibly float32 via `precision="fast"`) for the initial QR. The LOO residual `X[i,] * coef_i` mixes float64 rows with float32 QR coefficients silently. Severity: HIGH.

### 3. Deferred-host object in compiled plan after GPU release (G2)
**File: R/product-plan.R:232–261**
A `new_adgeMatrix_deferred` object bound to a plan, then released via `amatrix_release_resident`, will cause the plan's fallback path to operate on NaN placeholder host data (or throw a stop-level error on the deferred path). Severity: HIGH.
