# Hunter 03 — REPL surface
## (a) Drift check
- `bd search "print("`, `bd search "summary("`, `bd search "head("`, and `bd search "format("` returned no matching open issue.
- Follow-up duplicate checks `bd search "cov2cor"`, `bd search "crossprod"`, and `bd search "square numeric matrix"` also returned no matching open issue.

## (b) Scenario
- Fresh-process probes via `Rscript -e` with `pkgload::load_all('.')` and no attached `Matrix`.
- First-touch REPL probes on `adgeMatrix`: `print()`, `str()`, `summary()`, `format()`, `head()`, `tail()`, `as.character()`, `dput()`.
- Follow-up composed-workflow probe from the REPL surface:
- `X <- adgeMatrix(matrix(1:9 + 0.0, 3, 3))`
- `cov2cor(crossprod(X))`
- Controls:
- same expression on base `matrix`
- same expression on `Matrix::dgeMatrix`
- workaround `cov2cor(as.matrix(crossprod(X)))`

## (c) Findings
- Core REPL methods were mostly uneventful: `print`, `str`, `summary`, `head`, `tail`, `as.character`, and `dput` all returned without error.
- New bug found in the composed REPL workflow `cov2cor(crossprod(X))`.
- For `adgeMatrix`, `crossprod(X)` returns class `adgeMatrix`, `is.numeric(crossprod(X))` is `FALSE`, and `cov2cor(crossprod(X))` errors with `'V' is not a square numeric matrix`.
- The same workflow succeeds for `Matrix::dgeMatrix`, where `crossprod(X)` returns `dpoMatrix` and `cov2cor(...)` returns a `3 x 3` correlation matrix.
- The explicit escape hatch `cov2cor(as.matrix(crossprod(X)))` succeeds for `adgeMatrix`, so the failure is in the returned class / missing downstream compatibility path, not the numeric values.

## (d) Proposed bd create
- Filed `amatrix-af1`:
- `cov2cor(crossprod(adgeMatrix)) fails because crossprod returns non-numeric adgeMatrix`

## (e) Limitations
- I did not probe every downstream consumer of `crossprod(adgeMatrix)`; other base or stats helpers that expect numeric matrices may fail similarly.
- This report captures the user-visible symptom and comparative control, not the implementation fix.
