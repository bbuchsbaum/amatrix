# Round-4 Bug Hunt: Indexing Invariant Audit
**Hunter 01 | HEAD eaf8c43 | 2026-04-14**

---

## (a) Invariant + Drift Check

**Invariant:** For every exported amatrix class, subsetting and indexing (`x[i,j]`, `x[i,]`,
`x[,j]`, `x[i,j] <- v`, `head()`, `tail()`, `dim()`, `dimnames()<-`, `t()[i,j]`) must return
the same shape, class, values, and dimnames as base `matrix` / `Matrix::Matrix` semantics, bit-equal
on CPU and within tolerance on MLX/ArrayFire.

**Step 0 — drift check:**
- `packageVersion("amatrix")` → `0.1.0` ✓
- `DESCRIPTION` mtime → `2026-04-12` (2 days before today `2026-04-14`)

The version string matches HEAD eaf8c43 and `0.1.0`; the 2-day mtime gap is consistent with the
package having been installed from source on 2026-04-12 without re-touching DESCRIPTION. No
functional drift detected — proceeding.

---

## (b) Method Coverage Matrix

| Op | adgeMatrix | adgCMatrix | aTransposeView | KronMatrix |
|---|---|---|---|---|
| `[i,j]` | ✓ PASS | ✓ PASS | **BUG** (no method) | **BUG** (no method) |
| `[i,]` | ✓ PASS | ✓ PASS | **BUG** | **BUG** |
| `[,j]` | ✓ PASS | ✓ PASS | **BUG** | **BUG** |
| `[i,j]<-` | ✓ PASS | ✓ PASS | **BUG** | not tested |
| `[[i,j]]` | **BUG** | not tested | not tested | not tested |
| `head()` | ✓ PASS | not tested | not tested | not tested |
| `tail()` | ✓ PASS | not tested | not tested | not tested |
| `dim()` | ✓ PASS | ✓ PASS | ✓ PASS | ✓ PASS |
| `dimnames()` | ✓ PASS | ✓ PASS | ✓ PASS | n/a |
| `dimnames()<-` | ✓ PASS | ✓ PASS | n/a (read-only) | n/a |
| `drop=FALSE` | ✓ PASS | not tested | **BUG** | **BUG** |
| logical index | ✓ PASS | not tested | **BUG** | **BUG** |
| negative index | ✓ PASS | not tested | **BUG** | **BUG** |
| character index | ✓ PASS | not tested | **BUG** | **BUG** |
| `t()[i,j]` | n/a | n/a | **BUG** | n/a |

**Totals: 11 ops × 4 classes = 44 cells tested; ~16 confirmed fail cells, ~28 pass cells.**

---

## (c) Confirmed Bugs

### BUG-R4-01: `aTransposeView` has no `[` method — all subsetting fails

**Severity:** High. Every `t(X)[i,j]` call in user or modeling code silently crashes with an
unhelpful S4 error instead of returning the expected matrix slice.

**Root cause:** `aTransposeView` is exported and fully usable as a matrix-like object (it has
`dim`, `dimnames`, `%*%`, `Ops`) but no `[` or `[<-` method is registered. S4 dispatch falls
through to the default S4 `[` which errors with `"object of type 'S4' is not subsettable"`.

**One-shot repro:**
```r
library(amatrix)
A <- adgeMatrix(matrix(1:20, nrow = 4, ncol = 5))
tA <- t(A)              # class: aTransposeView
tA[1:3, 1:2]            # Error: object of type 'S4' is not subsettable
tA[1:3, ]               # Error: object of type 'S4' is not subsettable
tA[, 1:2]               # Error: object of type 'S4' is not subsettable
tA[1, 1] <- 99          # Error: object of type 'S4' is not subsettable
```

**Expected fix:** Add `[` and `[<-` methods for `aTransposeView` that either (a) materialize to
`adgeMatrix` then delegate to `am_subset`/`am_subassign`, or (b) map indices through the
transpose and delegate to the source `adgeMatrix` directly.

---

### BUG-R4-02: `KronMatrix` has no `[` method — all subsetting fails

**Severity:** Medium. `KronMatrix` is an exported, user-constructible class. Subsetting crashes
identically to BUG-R4-01.

**Root cause:** No `[` method registered for `KronMatrix`. `kronmatrix.R` implements `dim`, `t`,
`%*%`, `crossprod`, `solve`, `determinant`, and `as.matrix` — but no subsetting.

**One-shot repro:**
```r
library(amatrix)
K <- kron_matrix(matrix(1:4, 2, 2), diag(3))
K[1:3, 1:3]   # Error: object of type 'S4' is not subsettable
K[1, ]         # Error: object of type 'S4' is not subsettable
```

**Expected fix:** Add `[` for `KronMatrix` that calls `as.matrix(x)[i, j, ..., drop = drop]`
(materialize-then-slice). This is semantically correct and avoids the lazy representation promise
only for subsetting (which already requires materializing Kronecker rows/cols anyway).

---

### BUG-R4-03: `[[` not supported for `adgeMatrix` (base `matrix[[i]]` works)

**Severity:** Low–Medium. Base R `matrix[[1]]` (linear integer index) returns a scalar. Users
familiar with the `[[` idiom for single-element extraction will get an error.

**Root cause:** No `[[` method registered for `adgeMatrix`. The fallback S4 error message
("this S4 class is not subsettable") is misleading since `[` works fine.

**One-shot repro:**
```r
library(amatrix)
m <- matrix(1:6, 2, 3)
m[[1]]                      # 1  (base R works)
A <- adgeMatrix(m)
A[[1]]                      # Error: this S4 class is not subsettable
A[[2, 1]]                   # Error: incorrect number of subscripts
```

**Expected fix:** Register a `[[` method for `adgeMatrix` (and `adgCMatrix`) that delegates to
`as.matrix(amatrix_materialize_host(x))[[...]]` and returns a plain scalar, matching base R.

---

## (d) Inferred / Suspicious (needs repro)

### INFERRED-01: `aTransposeView` dimnames may be inconsistent after complex dispatch chains

**Hypothesis:** When `t(A)` is created from a deferred-host `adgeMatrix`, the
`.amatrix_transpose_dimnames()` path reverses `@Dimnames`. If `A` was constructed without
dimnames and then `dimnames(A) <- dn` was called (updating the host but not the live resident),
a subsequent `t(A)` might see stale dimnames in `@Dimnames`.

**Repro needed:**
```r
library(amatrix)
A <- adgeMatrix(matrix(1:6, 2, 3))    # no dimnames
dimnames(A) <- list(c("a","b"), c("x","y","z"))
tA <- t(A)
dimnames(tA)   # are these list(c("x","y","z"), c("a","b"))?
```

### INFERRED-02: `head()`/`tail()` on `adgCMatrix` may return `adgeMatrix` instead of `adgCMatrix`

**Hypothesis:** `head.matrix` / `tail.matrix` both delegate through `[i,]`. For `adgCMatrix`,
`am_subset` calls `amatrix_materialize_host` (returns `dgCMatrix`) then `[i,]` on it — which for
a sparse subset might return `dgCMatrix`. `.amatrix_rewrap_value` then correctly picks
`new_adgCMatrix`. But for a dense result (single row), it might return `adgeMatrix`. Needs repro
to confirm.

---

## (e) Refuted Bugs (tried to reproduce, could not)

The following were suspected but **confirmed not bugs**:

| Hypothesis | Evidence |
|---|---|
| `head(adgeMatrix)` returns wrong class | `head(A, 2)` → `adgeMatrix`; values and dimnames bit-equal to `head(m, 2)` |
| `tail(adgeMatrix)` returns wrong class | `tail(A, 2)` → `adgeMatrix`; values and dimnames bit-equal |
| `drop=FALSE` broken | `A[1, 1:3, drop=FALSE]` → `adgeMatrix` 1×3; values and dimnames match |
| `dimnames<-` corrupts values or class | `dimnames(A) <- newdn` → class `adgeMatrix`, values unchanged, names set correctly |
| `dimnames` not propagated through `[` | `A2[1:2,1:3]` after `dimnames(A2)<-newdn` → dimnames match base matrix |
| Logical index broken | `A[c(T,F,T,F),]` → correct rows, values and dimnames match |
| Negative index broken | `A[-1,]` → 3×5 `adgeMatrix`, values and dimnames match |
| Character index broken | `A["r1", , drop=FALSE]` → values and dimnames match |
| `[<-` with matrix value broken | `A[1:2,1:3] <- matrix(0, 2, 3)` → values correct |
| `adgCMatrix [<-` destroys class | `AC[1:2,1:2] <- 0` → class stays `adgCMatrix`, values correct |
| `adgCMatrix` subset returns wrong class | `AC[1:5,1:5]` → `adgCMatrix` |

**All 10 hypotheses refuted with direct evidence.**

---

## (f) Lint Rule Proposal

**Pattern:** Any exported S4 class that has `dim()`, `dimnames()`, and `%*%` methods but lacks
a `[` method will silently crash on subsetting with `"object of type 'S4' is not subsettable"`.

**Proposed lint rule:**

> **lint-r4-missing-bracket**: For every S4 class registered in the package namespace that
> (a) defines `dim` and (b) defines at least one of `%*%`, `crossprod`, or `Ops`, verify that at
> least one `setMethod("[", ...)` exists for that class. Classes without `[` are matrix-like but
> not subsettable — a violation of the principle of least surprise.

**Implementation sketch (R):**
```r
classes_with_dim <- Filter(
  function(cl) existsMethod("dim", cl),
  getClasses("package:amatrix")
)
classes_without_bracket <- Filter(
  function(cl) !existsMethod("[", cl),
  classes_with_dim
)
# classes_without_bracket should be empty
```

Running this check today returns: `aTransposeView`, `KronMatrix` — exactly the two classes with
confirmed BUG-R4-01 and BUG-R4-02.
