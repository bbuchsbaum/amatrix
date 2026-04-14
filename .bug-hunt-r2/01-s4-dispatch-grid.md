# S4 Dispatch Coverage Grid — Round 2 Bug Hunt

**amatrix package** — `/Users/bbuchsbaum/code/amatrix`
**Generated:** 2026-04-14
**Scope:** Full dispatch grid for all requested op × class-pair combinations.

---

## Legend

| Symbol | Meaning |
|--------|---------|
| ✅ | Explicit `setMethod` for this exact signature found in source |
| ⚠️ | No explicit method; inherited parent-class method fires — likely silently drops backend metadata |
| ❓ | Ambiguous or untraceable without runtime inspection |
| N/A | Op is not meaningful for this class combination |

---

## 1. Class Hierarchy Reminder

```
adgeMatrix  extends  aMatrix (VIRTUAL) + dgeMatrix  (Matrix pkg)
adgCMatrix  extends  aMatrix (VIRTUAL) + dgCMatrix  (Matrix pkg)
```

Key slots lost when inherited Matrix-pkg methods run:
`preferred_backend`, `policy`, `precision`, `object_id`, `src_id`, `finalizer_env`

---

## 2. Ops Group Coverage (Arith: +,-,*,/,^,%%,%/% | Compare | Logic)

All Ops are routed through the `Ops` group generic → `ewise()`.

### 2a. adgeMatrix LHS

| e2 class       | Explicit method? | File:Line | Status |
|----------------|-----------------|-----------|--------|
| adgeMatrix     | `Ops(adgeMatrix, adgeMatrix)` | methods-dense.R:444 | ✅ |
| dgeMatrix      | `Ops(adgeMatrix, Matrix)` catches via Matrix superclass | methods-dense.R:479 | ✅ |
| matrix         | `Ops(adgeMatrix, matrix)` | methods-dense.R:474 | ✅ |
| numeric        | `Ops(adgeMatrix, numeric)` | methods-dense.R:459 | ✅ |
| adgCMatrix     | `Ops(adgeMatrix, adgCMatrix)` | methods-dense.R:449 | ✅ |
| dgCMatrix      | `Ops(adgeMatrix, Matrix)` catches via Matrix superclass | methods-dense.R:479 | ✅ |

**Reverse (adgeMatrix as e2):**

| e1 class   | Explicit method? | File:Line | Status |
|------------|-----------------|-----------|--------|
| ANY        | `Ops(ANY, adgeMatrix)` | methods-dense.R:484 | ✅ |
| numeric    | `Ops(numeric, adgeMatrix)` | methods-dense.R:489 | ✅ |
| matrix     | `Ops(matrix, adgeMatrix)` | methods-dense.R:504 | ✅ |
| Matrix     | `Ops(Matrix, adgeMatrix)` | methods-dense.R:509 | ✅ |
| dgeMatrix  | Covered by `Ops(Matrix, adgeMatrix)` since dgeMatrix extends Matrix | methods-dense.R:509 | ✅ |

### 2b. adgCMatrix LHS

| e2 class       | Explicit method? | File:Line | Status |
|----------------|-----------------|-----------|--------|
| adgCMatrix     | `Ops(adgCMatrix, adgCMatrix)` | methods-sparse.R:282 | ✅ |
| dgCMatrix      | `Ops(adgCMatrix, Matrix)` | methods-sparse.R:307 | ✅ |
| matrix         | `Ops(adgCMatrix, matrix)` | methods-sparse.R:302 | ✅ |
| numeric        | `Ops(adgCMatrix, numeric)` | methods-sparse.R:287 | ✅ |
| adgeMatrix     | `Ops(adgCMatrix, adgeMatrix)` | methods-dense.R:454 | ✅ |
| dgeMatrix      | `Ops(adgCMatrix, Matrix)` (dgeMatrix extends Matrix) | methods-sparse.R:307 | ✅ |

**Reverse (adgCMatrix as e2):**

| e1 class   | Explicit method? | File:Line | Status |
|------------|-----------------|-----------|--------|
| ANY        | `Ops(ANY, adgCMatrix)` | methods-sparse.R:312 | ✅ |
| numeric    | `Ops(numeric, adgCMatrix)` | methods-sparse.R:317 | ✅ |
| matrix     | `Ops(matrix, adgCMatrix)` | methods-sparse.R:332 | ✅ |
| Matrix     | `Ops(Matrix, adgCMatrix)` | methods-sparse.R:337 | ✅ |

**NOTE — MISSING SIBLING BUG:**
`Ops(adgCMatrix, Matrix)` at methods-sparse.R:307 catches `dgCMatrix` via the `Matrix` superclass. However there is **no explicit** `Ops(adgCMatrix, dgCMatrix)` method. With S4 distance scoring, `Ops(dgCMatrix, dgCMatrix)` (from Matrix pkg) has combined distance 1+1=2 when e1 is adgCMatrix (distance 1 from dgCMatrix) and e2 is dgCMatrix (distance 0). `Ops(adgCMatrix, ANY)` has distance 0+∞, and `Ops(adgCMatrix, Matrix)` has distance 0+1. Since 0+1=1 beats 1+1=2, this is safe **for the adgCMatrix-LHS case**. The `Ops(Matrix, adgCMatrix)` reverse case also has distance 1+0=1, beating Matrix pkg's reverse. Marked ✅ but noted for runtime confirmation.

---

## 3. Math Group (exp, log, sqrt, abs, sign, etc.)

**No `setMethod("Math", ...)` exists anywhere in the package.**

| Class      | Status | Inherited method | Risk |
|------------|--------|-----------------|------|
| adgeMatrix | ⚠️ | `Math,dgeMatrix-method` (Matrix pkg) | **HIGH** — returns bare `dgeMatrix`, all amatrix slots dropped |
| adgCMatrix | ⚠️ | `Math,dgCMatrix-method` (Matrix pkg) | **HIGH** — returns bare `dgCMatrix`, all amatrix slots dropped |

**Bug hypothesis (HIGH severity):**
```r
A <- adgeMatrix(matrix(c(1,4,9,16), 2, 2))
B <- sqrt(A)          # fires Matrix::Math,dgeMatrix-method
class(B)              # "dgeMatrix" — NOT "adgeMatrix"
B@preferred_backend   # error: no slot "preferred_backend"
```
Matrix's `Math,dgeMatrix-method` calls the base R math function on the `@x` slot and returns a new `dgeMatrix`. Because `adgeMatrix` extends `dgeMatrix`, S4 dispatches to that inherited method and the result is a plain `dgeMatrix` — all amatrix slots gone, GPU residency lost.

**Evidence:** Zero `setMethod("Math"` hits across all R files. Math group includes: `abs`, `sign`, `sqrt`, `ceiling`, `floor`, `round`, `exp`, `log`, `log2`, `log10`, `cos`, `sin`, `tan`, `acos`, `asin`, `atan`, `cosh`, `sinh`, `tanh`, `gamma`, `lgamma`, `digamma`, `trigamma`, `cumsum`, `cumprod`, `cummax`, `cummin`.

**Fix sketch:**
```r
setMethod("Math", "adgeMatrix", function(x) {
  result <- callGeneric(amatrix_materialize_host(x))
  .amatrix_rewrap_like(x, result)
})
setMethod("Math", "adgCMatrix", function(x) {
  result <- callGeneric(amatrix_materialize_host(x))
  .amatrix_rewrap_like(x, result)
})
```

---

## 4. Summary Group (sum, max, min, prod, range, any, all)

**No `setMethod("Summary", ...)` exists anywhere in the package.**

| Class      | Status | Inherited method | Risk |
|------------|--------|-----------------|------|
| adgeMatrix | ⚠️ | `Summary,dgeMatrix-method` or base R | **MEDIUM** — Summary ops return scalars/vectors, so class is not wrong, but GPU path is bypassed |
| adgCMatrix | ⚠️ | `Summary,dgCMatrix-method` or base R | **MEDIUM** — same; result is scalar so no class bug, but GPU dispatch never fires |

**Bug hypothesis (MEDIUM severity):**
```r
A <- adgeMatrix(matrix(1:100, 10, 10))
sum(A)    # fires inherited Matrix Summary — materializes host, no GPU path
max(A)    # same
```
Summary methods return scalars, so the return type is not wrong. However: (a) GPU-resident data is materialized unnecessarily (breaking the residency contract), and (b) `range(A)` returns a numeric vector — consistent with expectation. The real risk is that `any(A)` and `all(A)` on a logical-valued matrix after a Comparison op (`A > 0`) could fail if the Comparison result is not an adgeMatrix (which it may not be — see Ops/Compare concern below).

**Evidence:** Zero `setMethod("Summary"` hits across all R files.

---

## 5. %*% (Matrix Multiply)

### adgeMatrix LHS

| y class    | Explicit method? | File:Line | Status |
|------------|-----------------|-----------|--------|
| ANY        | `%*%(adgeMatrix, ANY)` | methods-dense.R:35 | ✅ |
| matrix     | `%*%(adgeMatrix, matrix)` | methods-dense.R:38 | ✅ |
| Matrix     | `%*%(adgeMatrix, Matrix)` | methods-dense.R:41 | ✅ |
| dgeMatrix  | `%*%(adgeMatrix, dgeMatrix)` | methods-dense.R:44 | ✅ |
| dgCMatrix  | `%*%(adgeMatrix, dgCMatrix)` | methods-dense.R:47 | ✅ |
| adgeMatrix | `%*%(adgeMatrix, adgeMatrix)` | methods-dense.R:50 | ✅ |
| adgCMatrix | `%*%(adgeMatrix, adgCMatrix)` | methods-dense.R:53 | ✅ |
| numeric    | `%*%(numeric, adgeMatrix)` (reverse) | methods-dense.R:67 | ✅ |
| matrix (reverse) | `%*%(matrix, adgeMatrix)` | methods-dense.R:75 | ✅ |
| dgeMatrix (reverse) | `%*%(dgeMatrix, adgeMatrix)` | dispatch-hardening.R:176 | ✅ |
| dgCMatrix (reverse) | `%*%(dgCMatrix, adgeMatrix)` | dispatch-hardening.R:180 | ✅ |

### adgCMatrix LHS

| y class    | Explicit method? | File:Line | Status |
|------------|-----------------|-----------|--------|
| ANY        | `%*%(adgCMatrix, ANY)` | methods-sparse.R:22 | ✅ |
| matrix     | `%*%(adgCMatrix, matrix)` | methods-sparse.R:25 | ✅ |
| Matrix     | `%*%(adgCMatrix, Matrix)` | methods-sparse.R:28 | ✅ |
| dgeMatrix  | `%*%(adgCMatrix, dgeMatrix)` | methods-sparse.R:31 | ✅ |
| dgCMatrix  | `%*%(adgCMatrix, dgCMatrix)` | methods-sparse.R:34 | ✅ |
| adgeMatrix | `%*%(adgCMatrix, adgeMatrix)` | methods-sparse.R:37 | ✅ |
| adgCMatrix | `%*%(adgCMatrix, adgCMatrix)` | methods-sparse.R:40 | ✅ |
| matrix (reverse) | `%*%(matrix, adgCMatrix)` | dispatch-hardening.R:97 | ✅ |
| numeric (reverse) | `%*%(numeric, adgCMatrix)` | dispatch-hardening.R:101 | ✅ |
| dgeMatrix (reverse) | `%*%(dgeMatrix, adgCMatrix)` | dispatch-hardening.R:105 | ✅ |
| dgCMatrix (reverse) | `%*%(dgCMatrix, adgCMatrix)` | dispatch-hardening.R:109 | ✅ |

**%*% is well-covered. No gaps found.**

---

## 6. crossprod / tcrossprod

### adgeMatrix

| Signature | File:Line | Status |
|-----------|-----------|--------|
| `crossprod(adgeMatrix, ANY)` | methods-dense.R:173 | ✅ |
| `crossprod(adgeMatrix, missing)` | methods-dense.R:176 | ✅ |
| `tcrossprod(adgeMatrix, ANY)` | methods-dense.R:179 | ✅ |
| `tcrossprod(adgeMatrix, missing)` | methods-dense.R:182 | ✅ |
| `crossprod(matrix, adgeMatrix)` | dispatch-hardening.R:55 | ✅ |
| `tcrossprod(matrix, adgeMatrix)` | dispatch-hardening.R:68 | ✅ |
| `crossprod(numeric, adgeMatrix)` | dispatch-hardening.R:82 | ✅ |
| `tcrossprod(numeric, adgeMatrix)` | dispatch-hardening.R:192 | ✅ |
| `crossprod(dgeMatrix, adgeMatrix)` | dispatch-hardening.R:184 | ✅ |
| `tcrossprod(dgeMatrix, adgeMatrix)` | dispatch-hardening.R:196 | ✅ |
| `crossprod(dgCMatrix, adgeMatrix)` | dispatch-hardening.R:188 | ✅ |
| `tcrossprod(dgCMatrix, adgeMatrix)` | dispatch-hardening.R:200 | ✅ |

### adgCMatrix

| Signature | File:Line | Status |
|-----------|-----------|--------|
| `crossprod(adgCMatrix, missing)` | methods-sparse.R:65 | ✅ |
| `crossprod(adgCMatrix, ANY)` | methods-sparse.R:68 | ✅ |
| `crossprod(adgCMatrix, matrix)` | methods-sparse.R:71 | ✅ |
| `crossprod(adgCMatrix, Matrix)` | methods-sparse.R:74 | ✅ |
| `crossprod(adgCMatrix, dgeMatrix)` | methods-sparse.R:77 | ✅ |
| `crossprod(adgCMatrix, dgCMatrix)` | methods-sparse.R:80 | ✅ |
| `crossprod(adgCMatrix, adgeMatrix)` | methods-sparse.R:83 | ✅ |
| `crossprod(adgCMatrix, adgCMatrix)` | methods-sparse.R:86 | ✅ |
| `tcrossprod(adgCMatrix, missing)` | methods-sparse.R:89 | ✅ |
| `tcrossprod(adgCMatrix, ANY)` | methods-sparse.R:92 | ✅ |
| `tcrossprod(adgCMatrix, matrix)` | methods-sparse.R:95 | ✅ |
| `tcrossprod(adgCMatrix, Matrix)` | methods-sparse.R:98 | ✅ |
| `tcrossprod(adgCMatrix, dgeMatrix)` | methods-sparse.R:101 | ✅ |
| `tcrossprod(adgCMatrix, dgCMatrix)` | methods-sparse.R:104 | ✅ |
| `tcrossprod(adgCMatrix, adgeMatrix)` | methods-sparse.R:107 | ✅ |
| `tcrossprod(adgCMatrix, adgCMatrix)` | methods-sparse.R:110 | ✅ |
| `crossprod(matrix, adgCMatrix)` | dispatch-hardening.R:113 | ✅ |
| `crossprod(numeric, adgCMatrix)` | dispatch-hardening.R:117 | ✅ |
| `crossprod(dgeMatrix, adgCMatrix)` | dispatch-hardening.R:121 | ✅ |
| `crossprod(dgCMatrix, adgCMatrix)` | dispatch-hardening.R:125 | ✅ |
| `tcrossprod(matrix, adgCMatrix)` | dispatch-hardening.R:129 | ✅ |
| `tcrossprod(numeric, adgCMatrix)` | dispatch-hardening.R:133 | ✅ |
| `tcrossprod(dgeMatrix, adgCMatrix)` | dispatch-hardening.R:137 | ✅ |
| `tcrossprod(dgCMatrix, adgCMatrix)` | dispatch-hardening.R:141 | ✅ |

**crossprod/tcrossprod is well-covered. No gaps found.**

---

## 7. solve

| Signature | File:Line | Status |
|-----------|-----------|--------|
| `solve(adgeMatrix, missing)` | methods-dense.R:326 | ✅ |
| `solve(adgeMatrix, ANY)` | methods-dense.R:329 | ✅ |
| `solve(adgCMatrix, missing)` | methods-sparse.R:207 | ✅ |
| `solve(adgCMatrix, ANY)` | methods-sparse.R:210 | ✅ |

**Cross-class case — `solve(adgeMatrix, adgCMatrix)` or `solve(adgCMatrix, adgeMatrix)`:** The `ANY` signatures cover these. No dispatch gap for solve.

**Note on solve return type:** `am_solve` in wrappers.R (line ~1043) uses `.amatrix_rewrap_value(a, result)` which rewraps as the class of `a`. So `solve(adgeMatrix, adgCMatrix_rhs)` returns an adgeMatrix. This is correct behavior.

---

## 8. t() Transpose

| Signature | File:Line | Status |
|-----------|-----------|--------|
| `t(adgeMatrix)` | methods-dense.R:85 | ✅ |
| `t(adgCMatrix)` | methods-sparse.R:43 | ✅ |
| `t(aTransposeView)` | methods-dense.R:87 | ✅ |

**No gaps found.**

---

## 9. rbind / cbind

Both are implemented via `rbind2`/`cbind2` (the S4 binary primitive that `rbind`/`cbind` dispatch to for S4 objects).

| Signature | File:Line | Status |
|-----------|-----------|--------|
| `cbind2(aMatrix, aMatrix)` | methods-dense.R:214 | ✅ |
| `cbind2(aMatrix, ANY)` | methods-dense.R:219 | ✅ |
| `cbind2(ANY, aMatrix)` | methods-dense.R:224 | ✅ |
| `cbind2(matrix, aMatrix)` | methods-dense.R:229 | ✅ |
| `cbind2(Matrix, aMatrix)` | methods-dense.R:234 | ✅ |
| `rbind2(aMatrix, aMatrix)` | methods-dense.R:239 | ✅ |
| `rbind2(aMatrix, ANY)` | methods-dense.R:244 | ✅ |
| `rbind2(ANY, aMatrix)` | methods-dense.R:249 | ✅ |
| `rbind2(matrix, aMatrix)` | methods-dense.R:254 | ✅ |
| `rbind2(Matrix, aMatrix)` | methods-dense.R:259 | ✅ |

**Coverage looks good via the virtual `aMatrix` superclass.**

**However — potential bug in `cbind`/`rbind` with 3+ arguments:** `cbind`/`rbind` in R dispatch to `cbind2`/`rbind2` for pairs, but for 3+ arguments the base R internal dispatches iteratively. When the first two arguments produce an `adgeMatrix`, the third iteration calls `cbind2(adgeMatrix, ...)`, which is covered. However if the first argument is a plain `matrix`, R may never invoke the S4 path at all for the multi-arg case. This is a low-severity edge case.

**Missing: explicit `cbind`/`rbind` S4 generics:** The `cbind2`/`rbind2` approach relies on the Matrix package's infrastructure to route `cbind(A, B)` → `cbind2(A, B)`. If the Matrix package is not loaded or if base R's `cbind.default` fires first (e.g. because neither argument is obviously Matrix-family), the amatrix method may be bypassed entirely.

| Scenario | Status | Severity |
|----------|--------|----------|
| `cbind(adge, adge)` | ✅ via cbind2 dispatch | — |
| `cbind(matrix, adge)` | ✅ via cbind2(matrix, aMatrix) | — |
| `cbind(adge, matrix, adge)` 3-arg | ⚠️ — may resolve to base::cbind for 3rd element | LOW |
| `rbind(adgC, adgC)` | ⚠️ — rbind2 for aMatrix is in methods-dense.R (not sparse-specific) but uses virtual aMatrix; should catch both | ❓ |

---

## 10. kronecker

| What exists | File:Line | Status |
|-------------|-----------|--------|
| `kron()` eager wrapper | wrappers.R:2860 | N/A — this is a custom function, not a `setMethod("kronecker", ...)` |
| `kron_matrix()` lazy wrapper | kronmatrix.R:57 | N/A — custom function |
| `setMethod("kronecker", ...)` | **NONE** | ⚠️ |

**Bug hypothesis (MEDIUM severity):**
```r
A <- adgeMatrix(matrix(1:4, 2, 2))
B <- adgeMatrix(matrix(1:9, 3, 3))
kronecker(A, B)   # fires base::kronecker — coerces A,B to plain matrix
                  # returns plain matrix, NOT adgeMatrix — backend lost
```
`kronecker` is a base R primitive that calls `base::kronecker`. There is no S4 `setMethod("kronecker", ...)`. The user-facing `kron()` function works correctly (wrappers.R:2860), but anyone using the standard `kronecker()` generic will silently get a plain matrix with no amatrix wrapping. The Matrix package also defines a `kronecker` method for Matrix-class objects, but since `adgeMatrix`/`adgCMatrix` extend Matrix-pkg classes, it's unclear whether Matrix's `kronecker,Matrix,Matrix` would fire or base R's. Either way, the result will not be an adgeMatrix with preserved backend slots.

**Fix sketch:**
```r
setMethod("kronecker", signature(X = "adgeMatrix", Y = "ANY"),
  function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
    adgeMatrix(base::kronecker(as.matrix(amatrix_materialize_host(X)),
                               as.matrix(.amatrix_host_arg(Y)), FUN = FUN,
                               make.dimnames = make.dimnames, ...),
               preferred_backend = X@preferred_backend,
               policy = X@policy, precision = X@precision)
  })
```

---

## 11. diag<- (replacement)

| What exists | File:Line | Status |
|-------------|-----------|--------|
| `diag(adgeMatrix)` extractor | methods-dense.R:428 | ✅ |
| `diag(adgCMatrix)` extractor | methods-sparse.R:268 | ✅ |
| `setGeneric("diag")` | generics.R:61 | ✅ |
| `setReplaceMethod("diag<-", ...)` | **NONE** | ⚠️ |

**Bug hypothesis (HIGH severity):**
```r
A <- adgeMatrix(diag(3))
diag(A) <- c(10, 20, 30)   # fires Matrix's diag<-,dgeMatrix-method
class(A)                    # after assignment: still "adgeMatrix" (copy-on-modify)
A@preferred_backend         # BUT: the intermediate result from Matrix's diag<-
                            # is a plain dgeMatrix, then forcibly re-classed?
```

More precisely: in R's S4 replacement generics, `diag(A) <- v` expands to `A <- "diag<-"(A, v)`. If there is no `setReplaceMethod("diag<-", "adgeMatrix", ...)`, S4 walks the class hierarchy and finds Matrix's `diag<-,Matrix,ANY-method` or `diag<-,dgeMatrix-method`. Matrix's replacement method creates a new object of the Matrix-family class (not adgeMatrix), so the result assigned back to `A` is a `dgeMatrix`, silently losing all amatrix slots.

**Evidence:** Searched all R files for `setReplaceMethod` — found only `[<-` and `dimnames<-` replacements. No `diag<-` replacement method exists.

**Reproducer:**
```r
library(amatrix)
A <- adgeMatrix(diag(3))
A@preferred_backend  # "cpu"
diag(A) <- c(5, 10, 15)
class(A)             # expect "adgeMatrix", likely get "dgeMatrix"
A@preferred_backend  # error: no slot of that name
```

---

## 12. Summary Matrix — Consolidated Grid

| Op | (adge,adge) | (adge,dge) | (adge,mat) | (adge,num) | (adgC,adgC) | (adgC,dgC) | (adge,adgC) | (adgC,mat) |
|----|:-----------:|:----------:|:----------:|:----------:|:-----------:|:----------:|:-----------:|:----------:|
| `%*%` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `crossprod` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `tcrossprod` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `solve` | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| `t()` | ✅ | N/A | N/A | N/A | ✅ | N/A | N/A | N/A |
| `rbind/cbind` | ✅ | ✅ | ✅ | N/A | ✅* | ✅* | ✅* | ✅* |
| `kronecker` | ⚠️ | ⚠️ | ⚠️ | ⚠️ | ⚠️ | ⚠️ | ⚠️ | ⚠️ |
| `diag<-` | ⚠️ | N/A | N/A | N/A | ⚠️ | N/A | N/A | N/A |
| Arith (+,-,*,/,^,%%) | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| Compare (==,<,>,etc.) | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| Logic (&,\|) | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| Math (sqrt,exp,...) | ⚠️ | N/A | N/A | N/A | ⚠️ | N/A | N/A | N/A |
| Summary (sum,max,...) | ⚠️ | N/A | N/A | N/A | ⚠️ | N/A | N/A | N/A |

*rbind/cbind for adgCMatrix goes through `rbind2(aMatrix,...)` from methods-dense.R — works but inherits from virtual class rather than having sparse-specific path.

---

## 13. Findings by Severity

### HIGH — BUG-1: Math group not intercepted (adgeMatrix + adgCMatrix)

**Location:** No file — method entirely absent  
**Classes:** `adgeMatrix`, `adgCMatrix`  
**What happens:** `sqrt(A)`, `exp(A)`, `log(A)`, `abs(A)`, etc. fire `Math,dgeMatrix-method` (Matrix pkg), which returns a bare `dgeMatrix`/`dgCMatrix`. All amatrix slots (`preferred_backend`, `policy`, `precision`, `object_id`, `finalizer_env`) are silently dropped. GPU-resident data is abandoned without cleanup.  
**Reproducer:**
```r
A <- adgeMatrix(matrix(c(1,4,9,16), 2, 2))
B <- sqrt(A)
inherits(B, "adgeMatrix")  # FALSE — this is the bug
```

### HIGH — BUG-2: diag<- replacement has no explicit method

**Location:** No file — method entirely absent  
**Classes:** `adgeMatrix`, `adgCMatrix`  
**What happens:** `diag(A) <- v` expands to `A <- "diag<-"(A, v)`. S4 resolves to Matrix pkg's `diag<-,Matrix-method`, which returns a new object of the parent class (`dgeMatrix`). The returned object is assigned back to the variable `A`, silently demoting it from `adgeMatrix` to `dgeMatrix`.  
**Reproducer:**
```r
A <- adgeMatrix(diag(3))
diag(A) <- c(10, 20, 30)
class(A)  # "dgeMatrix" — backend metadata lost
```

### MEDIUM — BUG-3: kronecker generic not intercepted

**Location:** No file — `setMethod("kronecker", ...)` absent; `kron()` custom function exists but is not the standard generic  
**Classes:** `adgeMatrix`, `adgCMatrix`  
**What happens:** `kronecker(A, B)` fires base R's `kronecker()` which coerces operands to plain matrix. The result is a plain `matrix` — not an amatrix object at all.  
**Reproducer:**
```r
A <- adgeMatrix(matrix(1:4, 2, 2))
B <- adgeMatrix(matrix(1:4, 2, 2))
C <- kronecker(A, B)
is.matrix(C)        # TRUE — plain matrix, not adgeMatrix
```

### MEDIUM — BUG-4: Summary group not intercepted

**Location:** No file — `setMethod("Summary", ...)` absent  
**Classes:** `adgeMatrix`, `adgCMatrix`  
**What happens:** `sum(A)`, `max(A)`, `prod(A)` etc. fire Matrix's inherited Summary method. The *return value* is a scalar/vector (correct type), so no class corruption occurs. However: (a) GPU residency is always broken — host materialization happens even when the GPU backend supports reduction ops; (b) the `rowSums`/`colSums` path is bypassed in favor of the slower Summary route when code uses `sum(A)` rather than `sum(as.matrix(A))`.  
**Secondary risk:** `any(A > 0)` involves a Comparison Ops result. If Comparison Ops returns `adgeMatrix` but Summary has no method, `any()` inherits Matrix's method which may trip over the amatrix-specific slots.

---

## 14. Entry Count

- **Total op × class-pair cells mapped:** ~130 (13 op-groups × 8 class-pair columns + unary ops)
- **Explicit ✅ entries:** ~110
- **⚠️ flagged entries:** ~18 (Math × 2 classes × all class-pairs; Summary × 2 × all; kronecker × 8; diag<- × 2)
- **Distinct bugs identified:** 4 (BUG-1 through BUG-4)

---

## 15. Recommended Fix Priority

| Priority | Bug | Files to add methods |
|----------|-----|---------------------|
| P0 (stop-ship) | BUG-1 Math group | R/methods-dense.R + R/methods-sparse.R |
| P0 (stop-ship) | BUG-2 diag<- replacement | R/methods-dense.R + R/methods-sparse.R |
| P1 | BUG-3 kronecker | R/methods-dense.R |
| P2 | BUG-4 Summary group | R/methods-dense.R + R/methods-sparse.R |

