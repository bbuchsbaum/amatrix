# Round 3 Bug Hunt — Hunter 06 — Fresh Invariant: NA / NaN / Inf Propagation

**amatrix package** — `/Users/bbuchsbaum/code/amatrix`
**Generated:** 2026-04-14
**Scope:** IEEE missing-value semantics across every exported op and the
resident/deferred paths, via a lens that was **not** audited in round 2.

---

## (a) Invariant Chosen + Justification

**Invariant (N):** *For every exported arithmetic/linear-algebra/reduction op,
the output of the backend path must match `base R` / `Matrix` IEEE semantics
for `NA_real_`, `NaN`, and `±Inf` inputs, bit-compatible on the CPU path and
within numerical tolerance on any GPU path. In particular:*

1. **Contamination:** if an input element is `NA` / `NaN`, every output
   element that reads it must be `NA` / `NaN`.
2. **Distinction:** `NA_real_` and `NaN` are not inter-convertible; the CPU
   path already blurs them (R's bitwise test), but amatrix must not *create*
   spurious NAs where base R produces a value.
3. **Reductions respect `na.rm`:** when an exported reducer takes `na.rm`,
   the resident GPU kernel must implement the same semantics, or the wrapper
   must force CPU fallback.
4. **Sentinel separation:** amatrix uses `rep(NaN, n)` (constructors.R:199)
   as the deferred-host sentinel. A *user-supplied* NaN matrix must not be
   mis-identified as "deferred + data lost" and silently re-materialized.

### Why this is high-yield and fresh

- **Not touched in round 2.** Hunters 01–06 covered dispatch, anti-patterns,
  residency coherence, composition, error taxonomy, and cache invalidation.
  None probed what happens when an input element is not a normal float.
- **IEEE / R semantics is the one invariant R users absolutely assume.** A
  GPU kernel that replaces `NaN*0 = NaN` with `NaN*0 = 0` (common in fast
  fused-multiply-add code), or a reduction that returns `0` instead of `NaN`,
  silently corrupts downstream statistics — exactly where users trust R most.
- **Massive blast radius.** 121 `stop()` sites plus 0 `anyNA` guards on the
  compute path (only sinkhorn and distance/kernel validate) means any NaN
  leakage reaches `lm_fit`, `covariance`, `correlation`, `am_rowsums`,
  `rowmeans`, `am_qr`, `am_chol`, `am_svd`, and every Ops path. Conformance
  tests **only** include a single 2×2 NA/NaN test for `Ops` — no reductions,
  no factorizations, no deferred, no residency.
- **Double-sentinel hazard.** Because the deferred-host path uses `NaN` as
  sentinel for "data only lives on GPU", *user-provided* NaN data can
  interact with residency bookkeeping in confusing ways.

---

## (b) Code Paths Audited

### B.1 — Deferred-host NaN sentinel

- `R/constructors.R:197–210` — `new_adgeMatrix_deferred` initialises the host
  `@x` with `rep(NaN, n)`. The flag that distinguishes "sentinel NaN" from
  "real NaN" is `finalizer_env$host_deferred`. The `@x` slot itself is NOT
  distinguishable.
- `R/residency.R:422–483` — `amatrix_materialize_dense`:
  - If `fenv$host_deferred` is TRUE, download from GPU or stop.
  - Otherwise, if `cache_state$host_cache_valid` is TRUE, return `@x`.
  - Otherwise, download from backend.

### B.2 — Exported reducers that take `na.rm`

- `R/wrappers.R:941` `rowsums(x, na.rm, dims)` →
  `R/wrappers.R:645` `.amatrix_try_resident_rowSums(..., na.rm, dims, ...)` →
  `backend$rowSums_resident(lhs$key, na.rm, dims)`.
- `R/wrappers.R:956` `colsums(...)` → same pattern.
- `R/wrappers.R:1679` `rowmeans(x, na.rm)` →
  - For `adgeMatrix`: `rowsums(x, na.rm = na.rm) / ncol(x)`.
  - For `adgCMatrix`: `Matrix::rowMeans(..., na.rm = na.rm)`.
- `R/wrappers.R:1692` `colmeans(x, na.rm)` → symmetric.
- `R/methods-sparse.R:133–142` — `rowSums,adgCMatrix-method` etc. pass
  `na.rm` straight through.
- `R/wrappers.R:1910` `rowsum.adgeMatrix(x, group, reorder, na.rm, ...)`
  **silently drops `na.rm`** (argument declared but never referenced in the
  body) and goes straight to `segment_sum` which has no NA story at all.
- `R/wrappers.R:1940` `rowsum.adgCMatrix` honours `na.rm` via a full
  densification only when `na.rm=TRUE` (inefficient but correct); for
  `na.rm=FALSE`, it goes into a C call `am_sparse_segment_sum_c` with no
  evidence of NaN propagation guarantees.

### B.3 — Exported aggregators with no `na.rm` path

- `R/models-lm.R:638` `colSums(x_host * as.double(weights)) / sum(weights)`
  — used inside `lm_fit` centering and weighted paths; `x_host` is produced
  from `amatrix_materialize_host` which never filters NaN.
- `R/models-lm.R:840` `cor_host[!is.finite(cor_host)] <- NA_real_` —
  `correlation()` post-hoc converts Inf/NaN to NA_real_. This is a *creation*
  of NA, not a preservation. It masks a real numerical issue (zero variance)
  but also masks a real NaN input silently.
- `R/qr.R:703` `matrix(NA_real_, nrow = p, ncol = ncol(y))` — the explicit
  QR back-substitution sentinel for singular columns. Combined with
  `qr.coef` in `lm_loo_cv` (qr-downdate.R:145), this NA can leak into
  `y_vec[[i]] - sum(X[i, ] * coef_i)`, producing `NA` for *that* LOO fold
  silently.
- `R/wrappers.R:2401` — `.am_metric_checked_matrix` rejects NA/NaN/Inf up
  front for `am_dist` / `am_kernel`. This is the *only* exported function
  that actively rejects the invariant's bad inputs. All other functions
  pass them through silently.

### B.4 — Backend-level NaN semantics

- `R/backend-cpu.R:1–138` — CPU backend wraps `base::` primitives, which
  preserve NA/NaN per IEEE. `rowSums/colSums` *do* accept `na.rm`. **CPU
  path is compliant.**
- No GPU backend source is in-tree. Instead, every GPU backend is called via
  opaque `backend$rowSums_resident(key, na.rm, dims)`. There is no contract
  test verifying NaN propagation on a resident buffer, and no call site
  asserts that the backend honors `na.rm`.

### B.5 — Conformance test coverage for the invariant

- `tests/testthat/test-conformance.R:1339–1349` — the **only** NA/NaN test
  in the entire conformance suite. It checks `x + y`, `x * y`, `x > y`,
  `x == y` for a 2×2 matrix, with no backend specified (i.e. CPU only).
  Result: NA semantics on non-CPU backends is not tested at all.
- `tests/testthat/test-cross-backend-conformance.R` — zero NA/NaN/Inf
  references (grep-confirmed).

### B.6 — Known base-R idioms that would break

These are idioms that "just work" in R and will silently mis-compute on
amatrix GPU paths:

| Idiom | Expected | Hazard |
|---|---|---|
| `rowSums(A)` where `A` has NA rows | whole-row NaN | Silent wrong answer if `rowSums_resident` does a naive reduction without NaN sanitize |
| `rowmeans(A, na.rm=TRUE)` for adge | `rowSums(A, na.rm=T) / ncol(A)` is **wrong** — ncol doesn't vary by row after dropping NAs | Wrong denominator per row |
| `mean(A)` where A has NaN | NaN | Dispatches to inherited Matrix Summary (round 2 BUG-4), materializes host, *correct* by accident |
| `A %*% B` with NaN in A | NaN rows | CPU correct; GPU depends on kernel — FMA paths vary |
| `chol(A)` with NaN in A | error "not positive definite" | GPU may silently produce NaN factor; amatrix cache stores NaN |
| `am_scale(A)` where a column is all-NaN | all NaN output | Division by NaN sd may yield 0/NaN mix |
| `lm_fit(X, y)` where `y` has NAs | `complete.cases` or NA coefs | `models-lm.R:611` only checks *weights* for NA; y is unchecked |
| deferred `A <- new_adgeMatrix_deferred(...)`, then `as.matrix(A)` after `amatrix_release_resident(A)` — if the user inserted a user-NaN row into A via `am_ewise_inplace` before release, release drops the GPU key, `host_x` is NULL, `@x` was NaN sentinel. `amatrix_materialize_dense` at residency.R:446 throws `"deferred adgeMatrix lost its GPU resident data"`. | user's real NaN row is reported as "data lost" | Diagnostic collides with legitimate data |

---

## (c) Violations Found (with confidence labels)

### V1 — `rowmeans(adgeMatrix, na.rm=TRUE)` wrong denominator — **CONFIRMED (read-only reasoning)**

**File:** `R/wrappers.R:1683–1684`
```r
if (inherits(x, "adgeMatrix")) {
  return(rowsums(x, na.rm = na.rm) / ncol(x))
}
```

`base::rowMeans(x, na.rm=TRUE)` divides row `i`'s NA-stripped sum by
`sum(!is.na(x[i, ]))`, **not** by `ncol(x)`. The amatrix dense path divides
by a constant `ncol(x)`, producing wrong means whenever any row contains at
least one NA. The adgCMatrix branch (line 1681) correctly delegates to
`Matrix::rowMeans`. `colmeans` (line 1697) has the symmetric bug with
`nrow(x)`.

**Severity:** HIGH — silent numerical wrong answer. No error, no warning.
Breaks any preprocessing pipeline that uses `rowmeans(A, na.rm=TRUE)` for
missing-data imputation or centering.

**Probe sketch (read-only):**
```r
m <- matrix(c(1,2,NA,4), 1, 4)
base::rowMeans(m, na.rm=TRUE)            # 2.333...  (= 7/3)
amatrix::rowmeans(adgeMatrix(m), na.rm=TRUE)  # 1.75   (= 7/4) — WRONG
```

---

### V2 — `rowsum.adgeMatrix` silently ignores `na.rm` — **CONFIRMED (read-only reasoning)**

**File:** `R/wrappers.R:1910–1933`
```r
rowsum.adgeMatrix <- function(x, group, reorder = TRUE, na.rm = FALSE, ...) {
  ...
  result <- segment_sum(x, labels, K)
  ...
}
```

`na.rm` is declared in the signature but **never referenced** in the body.
`segment_sum` is a GPU kernel with no documented NaN story. Users who write
`rowsum(A, g, na.rm=TRUE)` get a silent NaN-contaminated result exactly in
the cases where the base R behavior promises cleanup.

**Severity:** HIGH — `rowsum` is commonly used in SR-side preprocessing
and k-means cluster aggregation. Silent miscompute.

---

### V3 — Deferred-host NaN sentinel collides with user NaN data — **CONFIRMED (read-only reasoning, strong)**

**File:** `R/constructors.R:199`, `R/residency.R:422–450`

The deferred path marks "data lives only on device" by setting
`finalizer_env$host_deferred = TRUE` and filling `@x` with `NaN`. If the
user's matrix legitimately contains NaN, there is no way to distinguish
"this NaN came from device" from "this NaN is the sentinel" after a
`amatrix_release_resident` or a failed `resident_has`:

1. `new_adgeMatrix_deferred(dim=c(n,p))` → `@x = NaN*n*p`, `host_deferred=T`.
2. Backend op stores real matrix (with legitimate NaN row) under
   `resident_key = K`.
3. Crash, backend restart, or explicit `amatrix_release_resident(A)`.
4. `amatrix_materialize_dense(A)` hits residency.R:427: host_deferred=TRUE,
   fenv$host_x is NULL, tries to re-download. `resident_has(K)` returns
   FALSE (buffer freed). Flow falls to residency.R:446:
   `stop("deferred adgeMatrix lost its GPU resident data")`.

The stop is thrown *even when* `@x` still holds the sentinel NaN vector
that would be "technically retrievable". The diagnostic conflates
"device data was released" with "amatrix cannot recover your data", so a
user pipeline that treats NaN-as-data has no recourse.

More subtly, if the deferred path were ever to *skip* the stop and return
`@x`, the user would silently receive an all-NaN matrix with correct
dimnames and no error.

**Severity:** MEDIUM — correctness hazard limited to the deferred code path,
but that path is the default under `amatrix.defer_host=TRUE` benchmark mode.

---

### V4 — `rowSums`/`colSums` resident path forwards `na.rm` to backends without capability assertion — **INFERRED (high)**

**File:** `R/wrappers.R:645–663`
```r
.amatrix_try_resident_rowSums <- function(x, na.rm, dims, backend_name) {
  backend <- .amatrix_get_backend(backend_name)
  if (!.amatrix_backend_supports_resident_op(backend, "rowSums", x = x)) return(NULL)
  lhs <- .amatrix_prepare_resident_arg(x, backend_name, promote_amatrix = FALSE)
  if (is.null(lhs)) return(NULL)
  result <- backend$rowSums_resident(lhs$key, na.rm, dims)
  ...
}
```

`.amatrix_backend_supports_resident_op(backend, "rowSums", ...)` gates on
the *op name*, not on a capability like `"rowSums_resident_na_rm"`. Any
backend registered as supporting `rowSums` gets `na.rm` passed down
unconditionally. If a backend's kernel ignores the flag (the typical fused
reduction implementation), `na.rm=TRUE` is silently a no-op — user gets
`NaN` propagating into every downstream calculation.

The CPU fallback at line 950 correctly uses `Matrix::rowSums(..., na.rm)`,
so the bug is specifically "resident path is stricter than CPU fallback"
and hence **fails only on fast/GPU paths**.

**Severity:** HIGH — affects every GPU-resident pipeline that uses `na.rm`,
but symptom is backend-specific. The package has no conformance test that
would catch it.

**Confidence:** Inferred (cannot actually exercise a GPU backend here), but
the absence of a capability guard is visible in source.

---

### V5 — `correlation()` converts `!is.finite` to `NA_real_`, masking real NaN inputs — **CONFIRMED**

**File:** `R/models-lm.R:838–841`
```r
cor_host <- cov_host / scale_mat
cor_host[!is.finite(cor_host)] <- NA_real_
diag(cor_host) <- 1
```

This rewrite-to-NA is *after* `cov_host / scale_mat`. If the input `X`
already contained NaN, the cov would be NaN, the division would be NaN,
and this line would convert to `NA_real_`. The user sees NA (missing
data) where the actual story is "you gave me NaN inputs". Conversely, a
genuinely degenerate (zero-variance) column produces Inf, which is
*also* rewritten to NA — merging two distinct failure modes into one
silent output.

**Severity:** LOW-MEDIUM — correctness is arguably preserved, but the
NaN/Inf/NA conflation prevents callers from detecting upstream problems.

---

### V6 — `lm_fit` validates `weights` for NA but not `y` or `X` — **CONFIRMED**

**File:** `R/models-lm.R:611`
```r
if (!is.numeric(weights) || anyNA(weights) || any(weights < 0)) {
  stop("weights must be a numeric vector of non-missing non-negative values")
}
```

`X` and `y` are never checked for NA. If they contain NaN, the normal
equations path returns silently NaN coefficients; the QR path may silently
return `NA_real_` via `qr.coef`'s unpivot-fill (qr.R:507,
`.amatrix_explicit_qr_unpivot(..., fill = NA_real_)`). The function does
not refuse the input or report which rows / which columns caused the NA.

Downstream, `qr-downdate.R:145` computes `y_vec[[i]] - sum(X[i, ] * coef_i)`;
if `coef_i` has an `NA_real_` entry from the unpivot-fill, the *entire* LOO
residual for that fold becomes `NA`, and `mean(loo_resid^2)` at line 150
becomes `NA` **unless** the user passed `na.rm=TRUE` — which they cannot,
because `qr-downdate.R:150` calls `mean()` without an `na.rm` argument.

**Severity:** HIGH — an LOO CV that *should* report a poisoned fold instead
reports a single `NA` scalar MSE with no diagnostic.

---

### V7 — Sentinel NaN violates `host_cache_valid` semantics — **INFERRED**

**File:** `R/residency.R:452–454`
```r
if (isTRUE(fenv$cache_state$host_cache_valid)) {
  return(new("dgeMatrix", x = x@x, Dim = x@Dim, Dimnames = x@Dimnames, factors = x@factors))
}
```

Cross-reference with round 2 hunter 03 finding V1: `host_cache_valid` is
set TRUE at bind time regardless of whether `@x` was the deferred sentinel
or real data. If a deferred object is later bound as resident (same
`object_id`), the cache_state env is shared, and `host_cache_valid=TRUE` is
set — but `@x` is still `rep(NaN, n)`. The next `amatrix_materialize_dense`
enters the fast-path at line 452 and returns the NaN sentinel as if it were
the matrix.

Round 2 V1 (hunter 03) noted that stale `@x` is returned, but did not
recognize that in the deferred path, "stale `@x`" is specifically
`rep(NaN, n)`, not "the previous real data". The user sees an all-NaN
matrix with correct dim and dimnames, no warning.

**Severity:** HIGH in combination with hunter 03 V1 — compounding bug.

---

### V8 — Resident `rowSums_resident(key, na.rm, dims)` contract is not enforced on registration — **INFERRED**

**File:** `R/backend-registry.R:140–163`

Registration checks `capabilities`, `features`, `precision_modes` exist and
return character vectors, but does not check that `rowSums_resident` (or
any resident reducer) *accepts* an `na.rm` parameter, let alone honours it.
This is the registration-time companion of V4: even if a backend author
added `rowSums_resident`, there is no contract check that asserts
`rowSums_resident(key, na.rm=TRUE, dims)` differs from
`rowSums_resident(key, na.rm=FALSE, dims)` on a NaN-containing test matrix.

**Severity:** MEDIUM — enables the V4 bug class to land without review.

---

## (d) Recommended Fix Sketch

### Fix 1 — Correct `rowmeans` / `colmeans` denominators (V1)

```r
rowmeans <- function(x, na.rm = FALSE) {
  if (!isTRUE(na.rm)) {
    if (inherits(x, "adgeMatrix")) return(rowsums(x, na.rm=FALSE) / ncol(x))
    ...
  }
  # na.rm=TRUE path: must count non-NA entries per row
  x_host <- as.matrix(amatrix_materialize_host(x))
  base::rowMeans(x_host, na.rm = TRUE)
}
```

The cleanest implementation is: when `na.rm=TRUE`, always materialize to
host. Performance regresses on GPU, but correctness is restored.

### Fix 2 — Honour `na.rm` in `rowsum.adgeMatrix` (V2)

```r
rowsum.adgeMatrix <- function(x, group, reorder = TRUE, na.rm = FALSE, ...) {
  if (isTRUE(na.rm)) {
    return(base::rowsum(as.matrix(amatrix_materialize_host(x)),
                        group, reorder = reorder, na.rm = TRUE, ...))
  }
  ... existing segment_sum path ...
}
```

Same shape as the sparse method at line 1943.

### Fix 3 — Disambiguate deferred sentinel from user NaN (V3, V7)

The deferred path should NOT use `NaN` as the sentinel, or the sentinel
should live in a separate slot from `@x`:

```r
new_adgeMatrix_deferred <- function(dim, ...) {
  n <- as.integer(dim[1L]) * as.integer(dim[2L])
  object_id <- .amatrix_next_object_id()
  fenv <- .amatrix_make_finalizer_env(object_id)
  fenv$host_deferred <- TRUE
  fenv$host_x <- NULL
  new("adgeMatrix",
      x = double(0),        # NOT rep(NaN, n): empty vector is the sentinel
      Dim = as.integer(dim),
      ...)
}
```

Every consumer that reads `x@x` of a deferred object must go through
`amatrix_materialize_dense` (which already checks `host_deferred`), so the
empty-vector sentinel is strictly better: any accidental read fails loudly.

### Fix 4 — Gate `na.rm` on a capability bit (V4, V8)

At registration, probe the backend:

```r
.amatrix_probe_rowsums_nam_rm <- function(backend) {
  test <- matrix(c(1, NA, NA, 4), 2)
  key  <- backend$resident_store(test)
  result_rm <- backend$rowSums_resident(key, na.rm=TRUE,  dims=1L)
  result_nr <- backend$rowSums_resident(key, na.rm=FALSE, dims=1L)
  backend$resident_drop(key)
  !isTRUE(all.equal(result_rm, result_nr)) && all(is.finite(result_rm))
}
```

If the probe fails, strip `"rowSums"` from the capability set or force CPU
fallback for any `na.rm=TRUE` call. The wrapper at `wrappers.R:945` should
check: `if (isTRUE(na.rm) && !.backend_honours_na_rm(backend, "rowSums"))
   return(Matrix::rowSums(amatrix_materialize_host(x), na.rm=TRUE, dims))`.

### Fix 5 — NA validation at `lm_fit` (V6)

```r
if (anyNA(X_host) || anyNA(y_host)) {
  amatrix_abort("X and y must not contain NA/NaN",
                class = "amatrix_error_invalid_value")
}
```

And in `lm_loo_cv` (qr-downdate.R:150), change `mean(loo_resid^2)` to
`mean(loo_resid^2, na.rm=TRUE)` plus a diagnostic count of NA folds.

### Fix 6 — Add NA/NaN/Inf conformance tests across backends

Extend `tests/testthat/test-cross-backend-conformance.R` with a parametric
matrix that contains one NA row, one NaN row, and one Inf row, then
verifies, *per backend*:

- `as.matrix(A + 1)` preserves NA/NaN/Inf element-wise
- `rowsums(A, na.rm=TRUE)` matches `Matrix::rowSums` bit-for-bit on CPU and
  within tolerance on GPU
- `rowmeans(A, na.rm=TRUE)` matches `base::rowMeans` with correct
  per-row denominator
- `lm_fit(X, y)` refuses NA X or NA y with a typed error
- `correlation(X)` distinguishes "zero variance column" (Inf→NA) from
  "user provided NaN" (NaN passes through or errors)

This alone would have caught V1, V2, V5, V6.

---

## Severity Summary

| # | Invariant clause | File:Line | Severity | Confidence |
|---|---|---|---|---|
| V1 | Reductions respect `na.rm` | R/wrappers.R:1683–1700 | **HIGH** | confirmed |
| V2 | `na.rm` honoured in reducers | R/wrappers.R:1910–1933 | **HIGH** | confirmed |
| V3 | Sentinel separation | R/constructors.R:199 + residency.R:422–450 | MEDIUM | confirmed |
| V4 | `na.rm` forwarded without capability guard | R/wrappers.R:645–663 | **HIGH** | inferred |
| V5 | `!is.finite → NA_real_` masks input NaN | R/models-lm.R:840 | LOW-MED | confirmed |
| V6 | `lm_fit` validates only weights | R/models-lm.R:611 + qr.R:703 + qr-downdate.R:145,150 | **HIGH** | confirmed |
| V7 | Deferred + host_cache_valid = sentinel leak | R/residency.R:452 + constructors.R:199 | HIGH | inferred (compound) |
| V8 | Registration contract doesn't assert NaN semantics | R/backend-registry.R:140–163 | MEDIUM | inferred |

---

## Cross-References to Earlier Bug Hunts

- **Compounds with round 2 hunter 03 V1** (host_cache_valid stale): when the
  stale `@x` is the deferred NaN sentinel, the symptom changes from "last
  known real data" to "all NaN output" — strictly worse.
- **Compounds with round 2 hunter 01 BUG-4** (Summary group): `mean(A)` and
  `sum(A)` silently materialise host *and* do not distinguish NaN
  contamination from valid data, so `mean(A)` may accidentally work while
  `rowmeans(A, na.rm=TRUE)` is wrong — a user will trust the first and
  never suspect the second.
- **Compounds with amatrix-j5a / amatrix-cng** (error taxonomy): when V6
  fires, the user gets silent NA, not a typed `amatrix_error_invalid_value`.

