# Hunter 02 — formula interface
## (a) Drift check
- `bd search "model.matrix"`, `bd search "lm("`, and `bd search "formula interface"` returned no matching open issue.

## (b) Scenario
- Fresh-process probes via `Rscript -e` with `pkgload::load_all('.')`.
- Tested exact round-6 user-facing paths:
- `model.matrix(~ X)` where `X` is an `adgeMatrix`
- `df <- data.frame(y = y); df$X <- X; model.matrix(y ~ X, data = df)`
- `lm(y ~ X, data = df)`
- Control path `lm.fit(as.matrix(X), y)`
- Compared behavior against base `matrix` and `Matrix::dgeMatrix`.

## (c) Findings
- No amatrix-specific bug found.
- `model.matrix(~ X)` errors for `adgeMatrix` with `invalid type (S4) for variable 'X'`, but the same error occurs for `Matrix::dgeMatrix`.
- `df$X <- X; lm(y ~ X, data = df)` is not a valid control for dense matrix predictors because base R also rejects the resulting object (`'data' must be a data.frame, not a matrix or an array`).
- `lm.fit(as.matrix(X), y)` works as expected.

## (d) Proposed bd create
- None.

## (e) Limitations
- I did not probe higher-level wrappers such as `recipes`, `parsnip`, or `MatrixModels`; this pass stayed on base `stats` formula/model-frame surfaces.
- If amatrix wants first-class formula support, that is a feature/design gap, not a confirmed regression relative to `Matrix`.
