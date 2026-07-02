# Round 3 Bug Hunt — Method Sweep Report

**Generated:** 2026-07-01 21:33:17
**Total failures recorded:** 2

## Methods Tested

### adgeMatrix (dense)
Arith (+,-,*,/,^) mat×mat and mat×scalar and scalar×mat; Compare (==,!=,<,>,<=,>=) vs scalar; Math group: abs, sqrt, exp, log, log2, log10, ceiling, floor, sign, cos, sin, tan, cosh, sinh, tanh, cumsum, cumprod, cummax, cummin; Summary: sum, max, min, prod, range; rowSums, colSums, rowMeans, colMeans; dim, nrow, ncol; t(); [i,j,drop=FALSE] and [i,j] (no drop); [<-; dimnames<-, rownames<-, colnames<-; diag (extractor), diag<- (replacement); as.matrix, as.numeric, as.vector, as.array; %*%, crossprod, tcrossprod, solve; rbind, cbind; kronecker; norm (1/I/F/M), det, svd, qr, chol

### adgCMatrix (sparse)
Arith (+,-,*,/) mat×mat and mat×scalar; Compare (==,!=,<,>) vs scalar; Math group: abs, sign, ceiling, floor (sparsity-preserving); exp, cosh, cos, sin, tan (sparsity-breaking → dense); sqrt, log (on positive input); cumsum, cumprod; Summary: sum, max, min; rowSums, colSums, rowMeans, colMeans; t(); [i,j,drop=FALSE] and [i,j] (no drop); dimnames<-; diag (extractor), diag<- (replacement); as.matrix; %*%, crossprod, tcrossprod, solve; rbind, cbind; kronecker

---

## Failures by Symptom

### ERROR (2)

| Method | Expected | Actual | Note |
|--------|----------|--------|------|
| `det[dense/cpu]` | `no error` | `no applicable method for 'determinant' applied to an object of class "c('adgeMatrix', 'aMatrix', 'dgeMatrix', 'unpackedMatrix', 'ddenseMatrix', 'generalMatrix', 'dMatrix', 'denseMatrix', 'Matrix')"` |  |
| `det[dense/mlx]` | `no error` | `no applicable method for 'determinant' applied to an object of class "c('adgeMatrix', 'aMatrix', 'dgeMatrix', 'unpackedMatrix', 'ddenseMatrix', 'generalMatrix', 'dMatrix', 'denseMatrix', 'Matrix')"` |  |

---

## Confirmed Bugs (Round-2 Issues Executed)

Confirmed 0 failure(s) that match round-2 issue patterns.


---

## NEW Bugs (Not in Round-2 Issue List) — 2 found

- **NEW** `det[dense/cpu]` → ERROR (expected `no error`, got `no applicable method for 'determinant' applied to an object of class "c('adgeMatrix', 'aMatrix', 'dgeMatrix', 'unpackedMatrix', 'ddenseMatrix', 'generalMatrix', 'dMatrix', 'denseMatrix', 'Matrix')"`)
- **NEW** `det[dense/mlx]` → ERROR (expected `no error`, got `no applicable method for 'determinant' applied to an object of class "c('adgeMatrix', 'aMatrix', 'dgeMatrix', 'unpackedMatrix', 'ddenseMatrix', 'generalMatrix', 'dMatrix', 'denseMatrix', 'Matrix')"`)

---

## Additional Findings (Hypothesis Refutations)

The following round-2 hypotheses were WRONG — these methods work correctly:

- **Math group on adgeMatrix** (amatrix-86l): `abs`, `sqrt`, `exp`, `log`, `ceiling`, `floor`, `sign`, `cos`, `sin`, `tan`, `cosh`, `sinh`, `tanh` all return `adgeMatrix` with backend preserved. The bug hypothesis was incorrect for the dense class.
- **Math group on adgCMatrix (sparsity-preserving)**: `abs`, `sign`, `ceiling`, `floor` return `adgCMatrix`. Correct.
- **diag<- on adgeMatrix** (amatrix-j5a): `diag(A) <- v` preserves `adgeMatrix` class and `preferred_backend` slot. The replacement method works correctly.
- **diag<- on adgCMatrix** (amatrix-j5a): `diag(S) <- v` preserves `adgCMatrix` class and `preferred_backend` slot. Correct.

---

## Top 3 Failures by Impact

1. **`det[dense/cpu]`** — ERROR: expected `no error`, got `no applicable method for 'determinant' applied to an object of class "c('adgeMatrix', 'aMatrix', 'dgeMatrix', 'unpackedMatrix', 'ddenseMatrix', 'generalMatrix', 'dMatrix', 'denseMatrix', 'Matrix')"`
2. **`det[dense/mlx]`** — ERROR: expected `no error`, got `no applicable method for 'determinant' applied to an object of class "c('adgeMatrix', 'aMatrix', 'dgeMatrix', 'unpackedMatrix', 'ddenseMatrix', 'generalMatrix', 'dMatrix', 'denseMatrix', 'Matrix')"`

---

_End of Round 3 Method Sweep Report_
