# Hunter 03 — cbind/rbind scenario

## (a) Drift check

- `packageVersion("amatrix")`: 0.1.0
- DESCRIPTION mtime: Apr 12 14:55:52 2026; HEAD commit for DESCRIPTION: 2026-04-12 15:12:43 — DESCRIPTION is slightly behind HEAD (no version bump in latest commits). No blocking drift; proceed.

## (b) Scenario matrix

| # | Probe | class(result) | dim | inherits adgeMatrix | values match | Status |
|---|-------|--------------|-----|---------------------|--------------|--------|
| 1 | `cbind(X, X)` | adgeMatrix | 3x8 | TRUE | TRUE | OK |
| 2 | `cbind(X, matrix(0,3,2))` | **dgeMatrix** | 3x6 | **FALSE** | TRUE | BUG |
| 3 | `cbind(X, 1)` | **dgeMatrix** | 3x5 | **FALSE** | TRUE | BUG |
| 4 | `cbind(X, c(1,2,3))` | **dgeMatrix** | 3x5 | **FALSE** | TRUE | BUG |
| 5 | `rbind(X, matrix(0,2,4))` | **dgeMatrix** | 5x4 | **FALSE** | TRUE | BUG |
| 6 | `rbind(X, 1)` | **dgeMatrix** | 4x4 | **FALSE** | TRUE | BUG |
| 7 | `cbind(X, Matrix::Matrix(0,3,2))` | adgeMatrix | 3x6 | TRUE | TRUE | OK |
| 8 | `cbind(X, diag(3))` | **dgeMatrix** | 3x7 | **FALSE** | TRUE | BUG |
| 9 | `cbind(adgCMatrix(X), X)` | — | — | — | — | CRASH (constructor bug) |
| 10 | `rbind(adgCMatrix(X), adgCMatrix(X))` | — | — | — | — | CRASH (constructor bug) |
| 9f | `cbind(Xs, X)` (Xs via coerce) | adgeMatrix | 3x8 | TRUE | TRUE | OK |
| 10f | `rbind(Xs, Xs)` (Xs via coerce) | adgCMatrix | 6x4 | TRUE | TRUE | OK |
| 11 | `cbind2(X, X)` direct | adgeMatrix | 3x8 | TRUE | TRUE | OK |
| 12 | `existsMethod("cbind2", c("adgeMatrix","adgeMatrix"))` | — | — | — | — | FALSE (inherits aMatrix,aMatrix — works) |
| 13 | numeric values cbind(X,1) | match | — | — | TRUE | OK |
| 14 | residency after cbind(X,X) | adgeMatrix | — | TRUE | — | OK (no leak; not resident) |
| 15 | without explicit `library(Matrix)` | adgeMatrix | 3x8 | TRUE | — | OK (Matrix auto-loaded) |

## (c) Findings

### BUG 1 — Class demotion: cbind/rbind with base matrix/vector/scalar (amatrix-0qt)

**Probes 2, 3, 4, 5, 6, 8 all fail.**

`cbind(X, matrix(...))`, `cbind(X, 1)`, `cbind(X, c(1,2,3))`, `rbind(X, matrix(...))`, `rbind(X, 1)`, `cbind(X, diag(3))` all return `dgeMatrix` instead of `adgeMatrix`. Values are numerically correct.

Root cause: S4 dispatch for `cbind`/`rbind` with heterogeneous types goes through `methods::cbind2`. The `aMatrix,ANY` method fires via `.amatrix_bind2`, which materializes both sides to host. `methods::cbind2(dgeMatrix, matrix)` returns a `dgeMatrix`. `.amatrix_rewrap_value` should re-wrap this to `adgeMatrix` — but the condition `(inherits(value, "Matrix") || is.matrix(value)) && .amatrix_is_numeric_matrix_value(value)` is true, so `.amatrix_rewrap_like(template, value)` is called. That calls `new_adgeMatrix(value, ...)`. However inspection shows the returned class is still `dgeMatrix`, meaning `new_adgeMatrix` is being called but either the dispatch picked a different method (not going through `.amatrix_bind2` at all) or the rewrap path isn't reached.

Confirmed via `showMethods`: for `cbind2("adgeMatrix", "matrix")`, the method selected is `x="matrix", y="aMatrix"` (because `adgeMatrix` extends `Matrix` and that signature scores better distance than `aMatrix, ANY` with a base `matrix` on the other side). This means x and y are swapped — template is constructed from the wrong arg. But `.amatrix_template` picks the first `aMatrix` it sees, so `template` = the matrix (non-aMatrix), `y` = the adgeMatrix. Result: rewrap uses wrong template or skips entirely.

Actually the deeper issue: with `x=matrix, y=aMatrix` dispatch, inside `.amatrix_bind2`, `template = .amatrix_template(x=matrix_base, y=adgeMatrix)` — the second branch picks `y` (adgeMatrix) correctly. So template IS the adgeMatrix. Then `methods::cbind2(.amatrix_host_arg(matrix_base), .amatrix_host_arg(adgeMatrix))` = `cbind2(matrix, dgeMatrix)` → returns `dgeMatrix`. `.amatrix_rewrap_value(template=adgeMatrix, value=dgeMatrix)` → `inherits(dgeMatrix, "Matrix")` is TRUE, `is_numeric_matrix_value` is TRUE → calls `new_adgeMatrix(dgeMatrix, ...)`. This SHOULD work — yet the result is `dgeMatrix`. Need to check `new_adgeMatrix` with a `dgeMatrix` argument.

Summary: regardless of the exact sub-path, the observable bug is confirmed: base `matrix`/`vector`/`scalar` as second arg causes demotion. Filed as amatrix-0qt.

### BUG 2 — adgCMatrix() constructor crashes on adgeMatrix input (amatrix-dum)

`adgCMatrix(X)` where X is `adgeMatrix` throws `Error: x must be a base matrix or dgCMatrix`. `.amatrix_sparse_base()` does not handle `aMatrix`/`adgeMatrix` inputs — it only accepts base `matrix`, `dgCMatrix`, or `sparseMatrix`. Since `adgeMatrix` extends `dgeMatrix` (dense), it falls through to `stop()`.

Confirmed: `adgCMatrix(as.matrix(X))` works fine. The constructor needs to coerce `aMatrix` inputs via `as.matrix()` / `amatrix_materialize_host()` before passing to `.amatrix_sparse_base()`.

### Non-bugs confirmed

- Probe 1 (cbind same class), 7 (cbind with Matrix-class), 11 (direct cbind2), 12 (method inheritance), 13 (numeric values), 14 (residency clean), 15 (no-Matrix-library): all pass.
- Sparse bind (probes 9f, 10f) works correctly once `adgCMatrix` is constructed via coerce.
- No residency leaks after cbind.
- Matrix namespace is auto-loaded by amatrix; no silent demotion without `library(Matrix)`.

## (d) Proposed bd create

Filed:
- **amatrix-0qt** (P2): `[bug] cbind/rbind with base matrix or vector demotes adgeMatrix to dgeMatrix` — affects probes 2, 3, 4, 5, 6, 8. Numerically safe but breaks `inherits(result, "adgeMatrix")` checks and loses backend/policy/precision metadata.
- **amatrix-dum** (P2): `[bug] adgCMatrix(x) crashes with 'x must be a base matrix or dgCMatrix' when x is adgeMatrix` — blocks natural user workflow of converting dense amatrix to sparse.

## (e) Limitations

- Residency probe (14) did not force a true resident allocation (the `%*%` call may not have pushed X to device); `amatrix_residency_info` showed no live resident state, so residency leak after cbind could not be ruled in or out under true GPU residency.
- Probes 9/10 required workaround (coerce via `as(Matrix::Matrix(...), "adgCMatrix")`) because the public constructor is broken for adgeMatrix input — original probes could not be run as specified.
- No ArrayFire/GPU backend available in this environment; all tests run on CPU backend only.
