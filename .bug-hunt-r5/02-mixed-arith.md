# Hunter 02 — Mixed-operand arithmetic scenario

## (a) Drift check

- Installed `packageVersion("amatrix")`: `0.1.0`
- `DESCRIPTION` mtime: 2026-04-12 14:55:52 (file modified on disk)
- HEAD commit for `DESCRIPTION`: `26f3a20` 2026-04-12 15:12:43 — i.e. the on-disk mtime predates the last git commit by ~17 minutes. The installed package may not reflect the latest source; probes used `library(amatrix)` which loads the installed binary. Discrepancy noted; probes are valid against the installed version.

## (b) Scenario matrix (what tested)

All 15 probes run with `Rscript -e '...'` against installed amatrix 0.1.0.
X = `as_adgeMatrix(matrix(1:12, 3, 4))`, ref = `as.matrix(X)`.

| # | Expression | class(result) | values match | Notes |
|---|-----------|--------------|--------------|-------|
| 1 | `X + 1` | adgeMatrix | TRUE | OK |
| 2 | `X + 1L` | adgeMatrix | TRUE | OK |
| 3 | `X * c(1,2,3)` | adgeMatrix | TRUE | OK |
| 4 | `c(1,2,3) * X` | adgeMatrix | TRUE | OK |
| 5 | `X + matrix(1:12,3,4)` | adgeMatrix | TRUE | OK |
| 6 | `matrix(1:12,3,4) + X` | adgeMatrix | TRUE | OK |
| 7 | `X + Matrix::Matrix(1:12,3,4)` | adgeMatrix | TRUE | OK |
| 7r | `dge + X` (reverse) | adgeMatrix | TRUE | OK |
| 8 | `X * X` | adgeMatrix | TRUE | OK |
| 9 | `X + adgeMatrix(matrix(1:12,3,4))` | adgeMatrix | TRUE | OK |
| 10a | `X > 5` | **lgeMatrix** | TRUE | **BUG: class demotion** |
| 10b | `X == 0` | **lgeMatrix** | TRUE | **BUG: class demotion** |
| 11a | `-X` | adgeMatrix | TRUE | OK |
| 11b | `!X` | **lgeMatrix** | TRUE | **BUG: class demotion** |
| 12 | `X %*% c(1,2,3,4)` | **numeric** | TRUE | **BUG: shape lost** |
| 12b | `X %*% matrix(c(1,2,3,4),4,1)` | adgeMatrix | TRUE | OK (col-matrix ok) |
| 13 | `crossprod(X)` (just library(amatrix)) | adgeMatrix | TRUE | OK (Matrix auto-imported) |
| 13b | `tcrossprod(X)` | adgeMatrix | TRUE | OK |
| 14 | `X / 0` | adgeMatrix | TRUE (Inf present) | OK |
| 15 | `X * .Machine$integer.max` | adgeMatrix | TRUE | OK |
| — | `X * c(1,2,3,4)` (length mismatch) | ERROR (correct) | — | Good guard |
| — | residency after ops | cpu preserved | — | OK |

## (c) Findings

### BUG 1 — Comparison/logical Ops silently demote adgeMatrix to lgeMatrix (amatrix-ol8)

`X > 5`, `X == 0`, `!X` return a bare Matrix-package `lgeMatrix`. The result is **not** an `aMatrix` subclass — `is(r, "aMatrix")` is `FALSE`. The amatrix metadata (backend, policy, precision, object_id) is lost entirely.

Root cause: `.amatrix_rewrap_value()` (wrappers.R:53–58) guards with `.amatrix_is_numeric_matrix_value()`, which returns `FALSE` for logical matrices. So the logical `lgeMatrix` from the CPU fallback is returned unwrapped.

Cascading failure: `(X > 5) + 0L` returns `dgeMatrix` (Matrix-package class), not `adgeMatrix` — a double demotion.

### BUG 2 — `X %*% plain_vector` returns shapeless numeric, not column matrix (amatrix-qic)

`X %*% c(1,2,3,4)` returns `class="numeric"` with `dim=NULL`. The matrix shape is lost. Base R `matrix(1:12,3,4) %*% c(1,2,3,4)` returns a 3x1 matrix; the amatrix path instead produces a plain vector. The amatrix wrapper is also lost.

Root cause: the matmul fallback does not coerce a vector rhs to a column matrix; the CPU operation returns a base numeric vector, which fails the numeric-matrix guard in `.amatrix_rewrap_value()`.

## (d) Proposed bd create

Both issues filed:

- `amatrix-ol8` (P2) — `[bug] Comparison ops (>, ==, !) silently demote adgeMatrix to bare lgeMatrix`
- `amatrix-qic` (P2) — `[bug] X %*% plain_vector returns numeric, not matrix/adgeMatrix`

## (e) Limitations

- Only CPU backend tested (no MLX/ArrayFire available in this environment); GPU resident paths for `ewise` were not exercised.
- Probes used installed package binary; DESCRIPTION mtime vs HEAD discrepancy noted — latest source commits may have already addressed some issues.
- Vector recycling for non-conformant lengths raises an error (correct behaviour); exact recycling semantics (column-major) not exhaustively validated for non-trivial lengths.
- `!X` logical not considered separately above — it returns `lgeMatrix` (same demotion as `>` / `==`).
