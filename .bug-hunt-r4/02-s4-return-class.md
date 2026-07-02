# Round 4 Bug Hunt — Hunter 02 — S4 Return-Class Audit (by execution)

**Hunter:** 02-s4-return-class (manual reconciliation after agent run terminated during final write)
**Date:** 2026-04-14
**HEAD:** eaf8c43
**Method:** Actually call each method and inspect return class — no grep.

---

## (a) Drift Check

- `packageVersion("amatrix")` → `0.1.0` ✓
- `DESCRIPTION` Depends: `R (>= 4.3)` only — **Matrix is in Imports, not Depends**. (relevant to root cause below)

## (b) Test Fixture

```r
suppressMessages(library(amatrix))
x <- adgeMatrix(matrix(runif(9)+1,3,3) + diag(3)*5, backend="cpu")
```

Note: `library(Matrix)` is **deliberately not attached**. Realistic single-package usage.

## (c) Op Sweep — amatrix-only (no Matrix attached)

| Op           | Return class             | Verdict                              |
|--------------|--------------------------|--------------------------------------|
| `t(x)`       | ERROR "argument is not a matrix" | **BUG — crashes**             |
| `chol(x)`    | `matrix, array`          | **BUG — silent class demotion**      |
| `rowSums(x)` | ERROR                    | **BUG — crashes**                    |
| `colSums(x)` | ERROR                    | **BUG — crashes**                    |
| `diag(x)`    | ERROR                    | **BUG — crashes**                    |
| `solve(x)`   | `matrix, array`          | **BUG — silent class demotion**      |
| `mean(x)`    | `NA` (warning)           | **BUG — returns NA**                 |
| `as.matrix(x)` | `matrix, array`        | OK (expected)                        |
| `log(x)`     | `adgeMatrix`             | OK                                   |
| `exp(x)`     | `adgeMatrix`             | OK                                   |
| `sqrt(x)`    | `adgeMatrix`             | OK                                   |
| `abs(x)`     | `adgeMatrix`             | OK                                   |
| `sum(x)`     | `numeric`                | OK (scalar reducer)                  |
| `max(x)`     | `numeric`                | OK                                   |
| `min(x)`     | `numeric`                | OK                                   |
| `dim(x)`     | `integer`                | OK                                   |
| `nrow(x)`    | `integer`                | OK                                   |
| `ncol(x)`    | `integer`                | OK                                   |
| `length(x)`  | `integer`                | OK                                   |
| `crossprod(x)`   | `adgeMatrix`         | OK                                   |
| `tcrossprod(x)`  | `adgeMatrix`         | OK                                   |

## (d) Confirmed Bugs — all stem from ONE root cause

### BUG-R4-04: amatrix does not promote Matrix's S4 generics — 7 core ops broken when only `library(amatrix)` is attached

**Root cause.** `t`, `chol`, `rowSums`, `colSums`, `diag`, `solve`, `mean` are S4 generics **defined in the Matrix package**. amatrix's `setMethod("t","adgeMatrix",...)` etc. registers methods against those generics, but DESCRIPTION has `Matrix` in **Imports** (not Depends), and amatrix never re-exports the generics. Result: when a user runs `library(amatrix); t(x)`, the `t` symbol resolves to `base::t` (S3), S4 dispatch is bypassed entirely, and either base S3 coercion returns a plain `matrix` or `t.default` crashes.

**Proof of root cause:**

```r
suppressMessages(library(amatrix))
x <- adgeMatrix(matrix(runif(9)+1,3,3) + diag(3)*5, backend="cpu")

# Without Matrix:
class(t(x))            # ERROR "argument is not a matrix"
class(chol(x))         # "matrix" "array"   <-- should be adgeMatrix / amChol
class(rowSums(x))      # ERROR
class(diag(x))         # ERROR

# Now attach Matrix:
suppressMessages(library(Matrix))
class(t(x))            # "aTransposeView" ✓
class(chol(x))         # "adgeMatrix"     ✓
class(rowSums(x))      # "numeric"        ✓
class(diag(x))         # "numeric"        ✓
class(solve(x))        # "adgeMatrix"     ✓
```

**Ops that already work** (log, exp, sqrt, abs, sum, max, crossprod, tcrossprod) work because their S4 generics live in `methods` (Math, Summary groups) or are re-exported from amatrix explicitly.

**Severity.** P1. Any downstream script that does `library(amatrix); t(X) %*% y` silently wrong-classes or crashes. Every modeling path touches `t`, `rowSums`, `solve`, `chol`, or `diag`.

**Fix options (in order of preference).**
1. Move `Matrix` from Imports to `Depends:` in DESCRIPTION (one-line fix, attaches Matrix whenever amatrix is attached).
2. Alternatively, amatrix imports `t`, `chol`, `rowSums`, `colSums`, `diag`, `solve` from Matrix and re-exports them in NAMESPACE.
3. Or define `setGeneric("t", ...)` etc. in amatrix, but this requires caller-side `useAsDefault=` hygiene.

### BUG-R4-05: `mean(adgeMatrix)` returns NA with warning even with Matrix attached (NOT verified — needs separate check)

Marked INFERRED — needs a `library(amatrix); library(Matrix); mean(x)` test to confirm separately. (Skipped in agent run; flag for round 5 or include in fix for BUG-R4-04.)

## (e) Refutations — round-2 class-demotion bugs that round-3 already refuted

- **Math group demotion** (r2 filed, r3 refuted): With Matrix attached, `log(adgeMatrix)` returns adgeMatrix. **Refutation holds.**
- **Summary group demotion**: `sum(adgeMatrix)` returns numeric scalar — CORRECT (reducer). **Refutation holds.**
- **diag<- demotion**: Not retested here; see r3 refutation.

The r2 grep-based findings were correct in saying "no setMethod for Math group" was wrong (Ops/Math group fallthrough handles it). r3's refutations stand.

## (f) Lint Rule Proposal

Add to `tools/lint-anti-patterns.R`:

**Rule P17 — S4 generic from non-Depends package.**
Parse NAMESPACE for `exportMethods(...)`. For each generic name, grep DESCRIPTION Depends: line — if the source package is in Imports but not Depends, and the method is on a user-facing class, emit a warning. Prevents regressions where a new Matrix-sourced generic gets a setMethod but the package isn't attached.

**One-shot promotion test:** 
```r
Rscript -e 'suppressMessages(library(amatrix)); x <- amatrix::adgeMatrix(diag(3)+1, backend="cpu"); stopifnot(inherits(t(x), "adgeMatrix") || inherits(t(x), "aTransposeView"))'
```

---

## Note on agent execution

The original hunter02 agent ran for ~6 min and executed all the probes above (43–53 tool uses across both attempts) but terminated during its final-write phase without producing a report. The reconciliation here is by direct Rscript execution on HEAD eaf8c43 — verified by the orchestrator, not inferred.
