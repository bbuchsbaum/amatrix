# Hunter 04 — dimnames preservation

## (a) Drift check

- `bd search dimnames` found no existing dimnames-specific open issue.
- Initial bare `library(amatrix)` probe re-hit known root cause `amatrix-1ha`
  (`rowSums`, `colSums`, `t`, `%*%` dispatch through base/Matrix attach gap).
- Re-ran the dimnames probe with `library(Matrix)` loaded to isolate
  dimnames-preservation from the known generic-registration issue.

## (b) Scenario

- Fresh-process probe at `tmp/bug-hunt/dimnames/probe2.R`.
- Matrix: dense `3 x 4`, deterministic `set.seed(1)`, row names `r1:r3`,
  column names `c1:c4`.
- Exercised:
  `t()`, 2d subsets, row/column extraction, `rowSums`, `colSums`,
  `rowMeans`, `colMeans`, arithmetic, `%*%`, `crossprod`, `tcrossprod`,
  `head()`, `diag(crossprod(.))`, name-based indexing, and
  `dimnames<-` / `rownames<-` / `colnames<-`.

## (c) Findings

- No new dimnames-preservation bug reproduced in this scenario.
- With `Matrix` attached, dense `adgeMatrix` preserved row/column names across
  all exercised operations and matched `Matrix`/base host behavior.
- The only failures in the first probe were duplicates of known issue
  `amatrix-1ha`, not dimnames drift.

## (d) Proposed bd create

- None. Do not file a new issue from this scenario.

## (e) Limitations

- This was dense-only. I did not probe sparse `adgCMatrix`, deferred-host
  objects, or resident-handle dimnames.
- The scenario did not include cross-process serialization or GPU-resident
  objects; it only checked ordinary dense host-backed behavior.
