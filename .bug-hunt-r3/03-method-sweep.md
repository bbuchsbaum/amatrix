# Round 3 Bug Hunt — Method Sweep Report

**Generated:** 2026-04-14 08:25:41
**Total failures recorded:** 28

## Methods Tested

### adgeMatrix (dense)
Arith (+,-,*,/,^) mat×mat and mat×scalar and scalar×mat; Compare (==,!=,<,>,<=,>=) vs scalar; Math group: abs, sqrt, exp, log, log2, log10, ceiling, floor, sign, cos, sin, tan, cosh, sinh, tanh, cumsum, cumprod, cummax, cummin; Summary: sum, max, min, prod, range; rowSums, colSums, rowMeans, colMeans; dim, nrow, ncol; t(); [i,j,drop=FALSE] and [i,j] (no drop); [<-; dimnames<-, rownames<-, colnames<-; diag (extractor), diag<- (replacement); as.matrix, as.numeric, as.vector, as.array; %*%, crossprod, tcrossprod, solve; rbind, cbind; kronecker; norm (1/I/F/M), det, svd, qr, chol

### adgCMatrix (sparse)
Arith (+,-,*,/) mat×mat and mat×scalar; Compare (==,!=,<,>) vs scalar; Math group: abs, sign, ceiling, floor (sparsity-preserving); exp, cosh, cos, sin, tan (sparsity-breaking → dense); sqrt, log (on positive input); cumsum, cumprod; Summary: sum, max, min; rowSums, colSums, rowMeans, colMeans; t(); [i,j,drop=FALSE] and [i,j] (no drop); dimnames<-; diag (extractor), diag<- (replacement); as.matrix; %*%, crossprod, tcrossprod, solve; rbind, cbind; kronecker

---

## Failures by Symptom

### CLASS_DEMOTION (28)

| Method | Expected | Actual | Note |
|--------|----------|--------|------|
| `Compare_==[dense/cpu]` | `aMatrix-derived` | `lgeMatrix` | Compare returns lgeMatrix — aMatrix wrapping lost |
| `Compare_!=[dense/cpu]` | `aMatrix-derived` | `lgeMatrix` | Compare returns lgeMatrix — aMatrix wrapping lost |
| `Compare_<[dense/cpu]` | `aMatrix-derived` | `lgeMatrix` | Compare returns lgeMatrix — aMatrix wrapping lost |
| `Compare_>[dense/cpu]` | `aMatrix-derived` | `lgeMatrix` | Compare returns lgeMatrix — aMatrix wrapping lost |
| `Compare_<=[dense/cpu]` | `aMatrix-derived` | `lgeMatrix` | Compare returns lgeMatrix — aMatrix wrapping lost |
| `Compare_>=[dense/cpu]` | `aMatrix-derived` | `lgeMatrix` | Compare returns lgeMatrix — aMatrix wrapping lost |
| `[submat_nodrop][dense/cpu]` | `adgeMatrix` | `dgeMatrix` | [i,j] without drop=FALSE falls through to dgeMatrix |
| `kronecker[dense/cpu]` | `adgeMatrix` | `dgeMatrix` | confirmed amatrix-jnd: kronecker returns dgeMatrix not adgeMatrix |
| `Compare_==[sparse/cpu]` | `aMatrix-derived` | `lgCMatrix` |  |
| `Compare_!=[sparse/cpu]` | `aMatrix-derived` | `lgeMatrix` |  |
| `Compare_<[sparse/cpu]` | `aMatrix-derived` | `lgeMatrix` |  |
| `Compare_>[sparse/cpu]` | `aMatrix-derived` | `lgCMatrix` |  |
| `[submat_nodrop][sparse/cpu]` | `adgCMatrix` | `dgCMatrix` | [i,j] without drop=FALSE falls through to dgCMatrix |
| `kronecker[sparse/cpu]` | `aMatrix-derived` | `dgCMatrix` | confirmed amatrix-jnd: kronecker returns dgCMatrix not adgCMatrix |
| `Compare_==[dense/mlx]` | `aMatrix-derived` | `lgeMatrix` | Compare returns lgeMatrix — aMatrix wrapping lost |
| `Compare_!=[dense/mlx]` | `aMatrix-derived` | `lgeMatrix` | Compare returns lgeMatrix — aMatrix wrapping lost |
| `Compare_<[dense/mlx]` | `aMatrix-derived` | `lgeMatrix` | Compare returns lgeMatrix — aMatrix wrapping lost |
| `Compare_>[dense/mlx]` | `aMatrix-derived` | `lgeMatrix` | Compare returns lgeMatrix — aMatrix wrapping lost |
| `Compare_<=[dense/mlx]` | `aMatrix-derived` | `lgeMatrix` | Compare returns lgeMatrix — aMatrix wrapping lost |
| `Compare_>=[dense/mlx]` | `aMatrix-derived` | `lgeMatrix` | Compare returns lgeMatrix — aMatrix wrapping lost |
| `[submat_nodrop][dense/mlx]` | `adgeMatrix` | `dgeMatrix` | [i,j] without drop=FALSE falls through to dgeMatrix |
| `kronecker[dense/mlx]` | `adgeMatrix` | `dgeMatrix` | confirmed amatrix-jnd: kronecker returns dgeMatrix not adgeMatrix |
| `Compare_==[sparse/mlx]` | `aMatrix-derived` | `lgCMatrix` |  |
| `Compare_!=[sparse/mlx]` | `aMatrix-derived` | `lgeMatrix` |  |
| `Compare_<[sparse/mlx]` | `aMatrix-derived` | `lgeMatrix` |  |
| `Compare_>[sparse/mlx]` | `aMatrix-derived` | `lgCMatrix` |  |
| `[submat_nodrop][sparse/mlx]` | `adgCMatrix` | `dgCMatrix` | [i,j] without drop=FALSE falls through to dgCMatrix |
| `kronecker[sparse/mlx]` | `aMatrix-derived` | `dgCMatrix` | confirmed amatrix-jnd: kronecker returns dgCMatrix not adgCMatrix |

---

## Confirmed Bugs (Round-2 Issues Executed)

Confirmed 4 failure(s) that match round-2 issue patterns.

- **CONFIRMED** `kronecker[dense/cpu]` → CLASS_DEMOTION (expected `adgeMatrix`, got `dgeMatrix`) — confirmed amatrix-jnd: kronecker returns dgeMatrix not adgeMatrix
- **CONFIRMED** `kronecker[sparse/cpu]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `dgCMatrix`) — confirmed amatrix-jnd: kronecker returns dgCMatrix not adgCMatrix
- **CONFIRMED** `kronecker[dense/mlx]` → CLASS_DEMOTION (expected `adgeMatrix`, got `dgeMatrix`) — confirmed amatrix-jnd: kronecker returns dgeMatrix not adgeMatrix
- **CONFIRMED** `kronecker[sparse/mlx]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `dgCMatrix`) — confirmed amatrix-jnd: kronecker returns dgCMatrix not adgCMatrix

---

## NEW Bugs (Not in Round-2 Issue List) — 24 found

- **NEW** `Compare_==[dense/cpu]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgeMatrix`) — Compare returns lgeMatrix — aMatrix wrapping lost
- **NEW** `Compare_!=[dense/cpu]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgeMatrix`) — Compare returns lgeMatrix — aMatrix wrapping lost
- **NEW** `Compare_<[dense/cpu]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgeMatrix`) — Compare returns lgeMatrix — aMatrix wrapping lost
- **NEW** `Compare_>[dense/cpu]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgeMatrix`) — Compare returns lgeMatrix — aMatrix wrapping lost
- **NEW** `Compare_<=[dense/cpu]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgeMatrix`) — Compare returns lgeMatrix — aMatrix wrapping lost
- **NEW** `Compare_>=[dense/cpu]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgeMatrix`) — Compare returns lgeMatrix — aMatrix wrapping lost
- **NEW** `[submat_nodrop][dense/cpu]` → CLASS_DEMOTION (expected `adgeMatrix`, got `dgeMatrix`) — [i,j] without drop=FALSE falls through to dgeMatrix
- **NEW** `Compare_==[sparse/cpu]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgCMatrix`)
- **NEW** `Compare_!=[sparse/cpu]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgeMatrix`)
- **NEW** `Compare_<[sparse/cpu]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgeMatrix`)
- **NEW** `Compare_>[sparse/cpu]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgCMatrix`)
- **NEW** `[submat_nodrop][sparse/cpu]` → CLASS_DEMOTION (expected `adgCMatrix`, got `dgCMatrix`) — [i,j] without drop=FALSE falls through to dgCMatrix
- **NEW** `Compare_==[dense/mlx]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgeMatrix`) — Compare returns lgeMatrix — aMatrix wrapping lost
- **NEW** `Compare_!=[dense/mlx]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgeMatrix`) — Compare returns lgeMatrix — aMatrix wrapping lost
- **NEW** `Compare_<[dense/mlx]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgeMatrix`) — Compare returns lgeMatrix — aMatrix wrapping lost
- **NEW** `Compare_>[dense/mlx]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgeMatrix`) — Compare returns lgeMatrix — aMatrix wrapping lost
- **NEW** `Compare_<=[dense/mlx]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgeMatrix`) — Compare returns lgeMatrix — aMatrix wrapping lost
- **NEW** `Compare_>=[dense/mlx]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgeMatrix`) — Compare returns lgeMatrix — aMatrix wrapping lost
- **NEW** `[submat_nodrop][dense/mlx]` → CLASS_DEMOTION (expected `adgeMatrix`, got `dgeMatrix`) — [i,j] without drop=FALSE falls through to dgeMatrix
- **NEW** `Compare_==[sparse/mlx]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgCMatrix`)
- **NEW** `Compare_!=[sparse/mlx]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgeMatrix`)
- **NEW** `Compare_<[sparse/mlx]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgeMatrix`)
- **NEW** `Compare_>[sparse/mlx]` → CLASS_DEMOTION (expected `aMatrix-derived`, got `lgCMatrix`)
- **NEW** `[submat_nodrop][sparse/mlx]` → CLASS_DEMOTION (expected `adgCMatrix`, got `dgCMatrix`) — [i,j] without drop=FALSE falls through to dgCMatrix

---

## Additional Findings (Hypothesis Refutations)

The following round-2 hypotheses were WRONG — these methods work correctly:

- **Math group on adgeMatrix** (amatrix-86l): `abs`, `sqrt`, `exp`, `log`, `ceiling`, `floor`, `sign`, `cos`, `sin`, `tan`, `cosh`, `sinh`, `tanh` all return `adgeMatrix` with backend preserved. The bug hypothesis was incorrect for the dense class.
- **Math group on adgCMatrix (sparsity-preserving)**: `abs`, `sign`, `ceiling`, `floor` return `adgCMatrix`. Correct.
- **diag<- on adgeMatrix** (amatrix-j5a): `diag(A) <- v` preserves `adgeMatrix` class and `preferred_backend` slot. The replacement method works correctly.
- **diag<- on adgCMatrix** (amatrix-j5a): `diag(S) <- v` preserves `adgCMatrix` class and `preferred_backend` slot. Correct.

---

## Top 3 Failures by Impact

1. **`Compare_==[dense/cpu]`** — CLASS_DEMOTION: expected `aMatrix-derived`, got `lgeMatrix`. Compare returns lgeMatrix — aMatrix wrapping lost
2. **`Compare_!=[dense/cpu]`** — CLASS_DEMOTION: expected `aMatrix-derived`, got `lgeMatrix`. Compare returns lgeMatrix — aMatrix wrapping lost
3. **`Compare_<[dense/cpu]`** — CLASS_DEMOTION: expected `aMatrix-derived`, got `lgeMatrix`. Compare returns lgeMatrix — aMatrix wrapping lost

---

_End of Round 3 Method Sweep Report_
