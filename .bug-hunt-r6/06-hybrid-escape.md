# Hunter 06 — hybrid escape / coercion surface

## (a) Drift check

- Working tree is ahead of the installed package and has substantial unrelated edits.
- I used fresh `Rscript --vanilla` probes in two modes:
  - `pkgload::load_all(".")` for current-checkout coercion behavior.
  - `library(amatrix)` / `library(Matrix)` for user-path spot checks.
- Important drift outcome:
  - installed `amatrix_0.1.0` still shows the old deferred silent-`NaN` coercion behavior;
  - current checkout loaded via `pkgload::load_all()` does **not** reproduce the silent path.
- Proceeded anyway because the point of this pass was to separate current-checkout behavior from installed-build folklore.

## (b) Scenario

- Primary scenario: a user crosses the amatrix boundary via coercion or escape hatches rather than staying inside amatrix-native verbs.
- Probed dense and sparse objects through:
  - `as.matrix()`, `methods::as(, "matrix")`, `as.array()`, `as.numeric()`, `as.vector()`
  - `Matrix::Matrix()`, `data.matrix()`
  - one base-R downstream step after escape: `as.matrix(X) %*% v`
- Follow-up scenario: the deferred-object failure family from round 5, but through alternate coercers rather than only `as.matrix()`.
- Spot-check scenario: worker processes (`mclapply`, PSOCK) to see whether escape failures were really parallel bugs or just the known generic-attachment bug.

## (c) Findings

### H1 — Ordinary hybrid escape paths mostly hold

- Dense `adgeMatrix` behaved correctly through `as.matrix`, `methods::as(, "matrix")`, `as.array`, `as.numeric`, `as.vector`, `data.matrix`, and `as.matrix(X) %*% v`.
- Sparse `adgCMatrix` also behaved correctly through `as.matrix`, `methods::as(, "matrix")`, `as.numeric`, `as.vector`, and `data.matrix`.
- Names / dimnames were preserved on the successful coercion paths.

### H2 — `Matrix::Matrix(adgCMatrix)` is a false lead, not an amatrix bug

- Initial probe suggested `Matrix::Matrix(adgCMatrix(...))` silently densifies to `dgeMatrix`.
- Control probe showed plain `dgCMatrix` takes the same path: `Matrix::Matrix(dgCMatrix)` also returns `dgeMatrix`.
- Conclusion: this is Matrix-package behavior, not a net-new amatrix regression.

### H3 — Deferred coercion bug extends beyond `as.matrix`

- Installed-build probe (`library(amatrix)`) still shows the old failure family:
  - `as.matrix(new_adgeMatrix_deferred(...))` silently returns all `NaN`.
  - `as.array(new_adgeMatrix_deferred(...))` also silently returns all `NaN`.
  - `as.numeric()` and `as.vector()` on the same object error with `deferred adgeMatrix lost its GPU resident data`.
- Current-checkout probe (`pkgload::load_all(".")`) no longer reproduces the silent path:
  - `as.matrix(...)` errors;
  - `as.array(...)` errors;
  - so the working tree appears to have fixed or hardened this already.
- Conclusion: this pass found a broader historical symptom surface for the old deferred-coercion bug, but not a live regression in the current checkout.

### H4 — Parallel-worker failures were duplicates of `amatrix-1ha`, not new parallel bugs

- `mclapply(..., rowSums(x))` and PSOCK `parLapply(..., rowSums(x))` fail under `library(amatrix)` alone with:
  - `'x' must be an array of at least two dimensions`
- Adding `library(Matrix)` inside the worker makes the same probes pass.
- Conclusion: this is more evidence for the existing generic-attachment / `Depends` bug (`amatrix-1ha`), not a distinct worker-only defect.

## (d) Proposed bd create

- No net-new issue filed from this pass.
- No reopen recommendation from this pass because the current checkout does not reproduce the silent path under `pkgload::load_all()`.
- If release notes or changelog work is being done, the installed-vs-checkout discrepancy is worth documenting so users do not mistake an old installed build for a regression.

## (e) Limitations

- CPU-only environment. No MLX / ArrayFire resident objects were available, so I could not test whether escape hatches mishandle live device-backed sparse objects differently from CPU-backed ones.
- I did not probe every downstream consumer after escape (`stats`, `MatrixModels`, `irlba`, etc.) once the core coercions themselves looked healthy.
- The worker probes were intentionally minimal and used `rowSums` only; they were sufficient to distinguish “parallel bug” from “known missing generic attachment.”
- I did not bisect which exact local change fixed the deferred silent-coercion path; only the runtime fact of “installed reproduces, current checkout errors correctly” was established.
