# Hunter 04 — apply/NA/integer idiom scenario

## (a) Drift check

- `packageVersion("amatrix")`: 0.1.0
- DESCRIPTION mtime: 2025-04-12 14:55 (matches HEAD commit 2695f6e)
- No drift detected.

## (b) Scenario matrix

| Probe | Description | Expected | Result |
|-------|-------------|----------|--------|
| 1a | `rowMeans(X)` no na.rm | NA propagation | ERROR: 'x' must be an array |
| 1b | `rowMeans(X, na.rm=TRUE)` | c(3, 4) | ERROR: 'x' must be an array |
| 1c | `colMeans(X, na.rm=TRUE)` | c(1, 3.5, 5) | ERROR: 'x' must be an array |
| 1d | `rowSums(X, na.rm=TRUE)` | c(9, 4) | ERROR: 'x' must be an array |
| 1e | `colSums(X, na.rm=TRUE)` | c(1, 7, 5) | ERROR: 'x' must be an array |
| 2 | NaN rowMeans (no/with na.rm) | NaN propagation / na.rm removes NaN | ERROR (same dispatch failure) |
| 3a | `apply(X, 2, mean)` | NA, 3.5, NA | CORRECT (via `[,j]` dispatch) |
| 3b | `apply(X, 1, sum)` | 9, NA | CORRECT (via `[i,]` dispatch) |
| 3c | `apply(X, 2, function(v) sum(v^2))` | NA, 25, NA | CORRECT |
| 4 | Integer matrix `adgeMatrix(1:12 matrix)` | coerce to double | CORRECT (`typeof(Xi@x) == "double"`) |
| 4b | `Xi %*% c(1,2,3,4)` | c(70, 80, 90) | CORRECT |
| 4c | `rowSums(Xi)` integer matrix | c(22, 26, 30) | ERROR (same dispatch failure) |
| 5 | All-NA row `rowMeans(na.rm=TRUE)` | NaN, 2 | ERROR (same dispatch failure) |
| 6 | `X[X > 3]` logical subscript | c(4, 5, 6) | CORRECT |
| 7 | `X[is.na(X)] <- 0` logical sub-assignment | replaces NAs with 0 | ERROR: subscript out of bounds |
| 8 | Probe 1 without explicit `library(Matrix)` | same as probe 1 | same (Matrix auto-loaded by amatrix) |

## (c) Findings

### BUG 1 (P1): rowMeans/colMeans/rowSums/colSums S4 dispatch fails completely

**Repro:**
```r
library(amatrix)
M <- matrix(c(1,NA,3,4,5,NA), 2, 3)
X <- adgeMatrix(M)
rowMeans(X)           # ERROR: 'x' must be an array of at least two dimensions
rowSums(X, na.rm=TRUE) # same error
colMeans(X)           # same error
colSums(X)            # same error
```

**Root cause:** `rowSums`, `colSums`, `rowMeans`, `colMeans` are base R primitives, not S4 generics, from the user's perspective. `Matrix` provides S4 versions under `Matrix::rowSums` etc., but these live in the Matrix namespace. `amatrix` registers S4 methods (`setMethod("rowSums", "adgeMatrix", ...)`) and calls `setGeneric("rowSums")` in `generics.R`, but at runtime `isS4(rowSums)` is `FALSE` — the base primitive is found first on the search path, bypassing S4 dispatch entirely.

The direct internal functions `amatrix:::rowmeans(X)` and `amatrix:::rowsums(X)` work correctly. Only the S4-generic-dispatched public path (`rowMeans(X)`) is broken.

Affects: ALL reductions — `rowSums`, `colSums`, `rowMeans`, `colMeans` — for both NA and non-NA inputs.

**Filed:** `amatrix-juq` (P1)

### BUG 2 (P2): `X[is.na(X)] <- 0` logical sub-assignment fails with "subscript out of bounds"

**Repro:**
```r
library(amatrix)
M <- matrix(c(1,NA,3,4,5,NA), 2, 3)
X <- adgeMatrix(M)
X[is.na(X)] <- 0   # Error in X[is.na(X)] <- 0 : subscript out of bounds
```

**Expected:** replaces NAs in place, matching `M[is.na(M)] <- 0` semantics.

`is.na(X)` itself returns a logical matrix correctly. The sub-assignment handler `am_subassign` materialises the host and assigns, but when `i` is a logical matrix (not a vector and not a row index), the two-arg `host_x[i, , ...]` form is used incorrectly — it interprets the logical matrix as a row subscript instead of a flat logical index.

**Filed:** `amatrix-e97` (P2)

### Passing probes

- `apply(X, 2, mean)` / `apply(X, 1, sum)` / `apply(X, 2, fun)` — all correct (dispatch works via `[,j]`/`[i,]` which extract numeric vectors).
- Integer matrix input: `adgeMatrix(matrix(1L, ...))` coerces to double silently as expected; `%*%` gives correct result.
- `X[X > 3]` logical extraction — correct.

## (d) Proposed bd create

Both issues filed:
- `amatrix-juq` (P1): rowMeans/colMeans/rowSums/colSums S4 dispatch fails — base primitive bypasses S4 methods
- `amatrix-e97` (P2): `X[is.na(X)] <- 0` logical sub-assignment errors with subscript out of bounds

## (e) Limitations

- Did not test `adgCMatrix` (sparse) path for these reductions — may have same or different failure mode.
- Did not test `rowMeans` with `dims > 1` parameter.
- The sub-assignment bug was only tested with a logical-matrix index; scalar/integer index sub-assignment may work.
- NaN vs NA distinction could not be probed due to the upstream dispatch failure (both fail identically).
