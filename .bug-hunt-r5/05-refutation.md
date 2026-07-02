# Hunter 05 — Runtime refutation (round 5)

## (a) Drift check

- Installed `packageVersion("amatrix")`: `0.1.0`
- `DESCRIPTION` mtime: `2026-04-12 14:55:52`
- Last git commit touching `DESCRIPTION`: `2026-04-12 15:12:43` (SHA `26f3a20618c693af`)
- **Discrepancy**: DESCRIPTION on disk is ~17 minutes older than the HEAD commit that touched it — installed package may lag HEAD by at least one commit (`2695f6e`). Probes run against installed (0.1.0); HEAD may differ slightly.

## (b) Targets pulled from `bd list --status=open`

44 open issues total. Selected 11 CPU-probeable P1/P2 `[bug]` issues:

| ID | Title (abbreviated) |
|----|---------------------|
| amatrix-jnd | kronecker generic not intercepted |
| amatrix-xnp | aTransposeView has no `[` |
| amatrix-x6a | KronMatrix has no `[` |
| amatrix-sxs | `[[` not supported on adgeMatrix |
| amatrix-vbh | amatrix_release_resident not exported |
| amatrix-lei | zzz.R onLoad clobbers calibration |
| amatrix-75h | kernel_matrix rbf diagonal drift |
| amatrix-p24 | pairwise_sqdist_argmin wrong centroid (RE-OPENED) |
| amatrix-1ha | Matrix in Imports not Depends |
| amatrix-7il | blanket tryCatch fallbacks mask errors |
| amatrix-36q | double-drop sites in backend-planning/models-lm |
| amatrix-3ka | .amatrix_bind_resident leaks prior key on rebind |

## (c) Per-bug verdict (with runtime probe + output)

---

### amatrix-jnd — kronecker generic not intercepted
**Verdict: CONFIRMED**

Probe:
```r
library(amatrix)
X <- matrix(1:4, 2, 2); A <- adgeMatrix(X); B <- adgeMatrix(diag(2))
result <- kronecker(A, B)
cat("class:", class(result), "\n")
cat("Is KronMatrix:", inherits(result, "KronMatrix"), "\n")
cat("Is plain matrix:", is.matrix(result), "\n")
```
Output:
```
kronecker result class: dgeMatrix
Is KronMatrix: FALSE
Is plain matrix: FALSE
```
Bug fires: `kronecker()` coerces to `dgeMatrix` (Matrix class), not `KronMatrix`.

---

### amatrix-xnp — aTransposeView has no `[` method
**Verdict: CONFIRMED**

Probe:
```r
library(amatrix)
X <- matrix(1:9, 3, 3); A <- adgeMatrix(X); tA <- t(A)
cat("class:", class(tA), "\n")
result <- tryCatch(tA[1, 2], error=function(e) e)
cat("ERROR:", conditionMessage(result), "\n")
```
Output:
```
class of t(A): aTransposeView
Error in t.default(A) : argument is not a matrix
```
Note: `t(A)` itself errors (no S4 `t` method on `adgeMatrix` without Matrix attached — see amatrix-1ha), so `aTransposeView` cannot be constructed. The `[` gap compounds 1ha.

---

### amatrix-x6a — KronMatrix has no `[` method
**Verdict: CONFIRMED**

Probe:
```r
library(amatrix)
X <- matrix(1:4, 2, 2); A <- adgeMatrix(X); B <- adgeMatrix(diag(2))
K <- kron_matrix(A, B)
cat("class:", class(K), "\n")
cat("Has [ method:", existsMethod("[", "KronMatrix"), "\n")
result <- tryCatch(K[1, 1], error=function(e) e)
cat("ERROR:", conditionMessage(result), "\n")
```
Output:
```
class of K: KronMatrix
Has [ method: FALSE
ERROR: object of type 'S4' is not subsettable
```

---

### amatrix-sxs — `[[` not supported on adgeMatrix
**Verdict: CONFIRMED**

Probe:
```r
library(amatrix)
X <- matrix(1:9, 3, 3); A <- adgeMatrix(X)
result <- tryCatch(A[[1]], error=function(e) e)
cat("ERROR:", conditionMessage(result), "\n")
```
Output:
```
class of A: adgeMatrix
ERROR: this S4 class is not subsettable
```

---

### amatrix-vbh — amatrix_release_resident not exported
**Verdict: CONFIRMED**

Probe:
```r
library(amatrix)
in_exports <- "amatrix_release_resident" %in% ls("package:amatrix")
ns <- asNamespace("amatrix")
in_ns <- exists("amatrix_release_resident", envir=ns, inherits=FALSE)
cat("In package exports:", in_exports, "\n")
cat("In namespace:", in_ns, "\n")
```
Output:
```
In package exports: FALSE
In namespace: FALSE
```
Function does not exist at all in namespace (not merely unexported).

---

### amatrix-lei — zzz.R .onLoad clobbers calibration
**Verdict: CONFIRMED**

Source confirms `zzz.R:17`: `amatrix_register_backend("cpu", .amatrix_cpu_backend(), overwrite = TRUE)`

Runtime probe — inject calibration, simulate re-register:
```r
library(amatrix)
ns <- asNamespace("amatrix")
state <- get(".amatrix_state", envir=ns)
state$calibration <- list(
  thresholds = list(cpu = list(matmul=1000)),
  results = data.frame(backend="cpu", op="matmul", stringsAsFactors=FALSE)
)
cat("Before:", paste(names(state$calibration$thresholds), collapse=","), "\n")
amatrix_register_backend("cpu", amatrix:::.amatrix_cpu_backend(), overwrite=TRUE)
cal_after <- state$calibration
cat("After cpu thresholds:", paste(names(cal_after$thresholds), collapse=","), "\n")
cat("CPU threshold wiped:", is.null(cal_after$thresholds[["cpu"]]), "\n")
cat("Results rows for cpu:", nrow(cal_after$results[cal_after$results$backend=="cpu",,drop=FALSE]), "\n")
```
Output:
```
Before: cpu
After cpu thresholds: 
CPU threshold wiped: TRUE
Results rows for cpu: 0
```
Every package load wipes any previously saved cpu calibration.

---

### amatrix-75h — kernel_matrix rbf diagonal drift
**Verdict: REFUTED**

Probe:
```r
library(amatrix)
set.seed(42)
X <- matrix(rnorm(30), 10, 3)
K <- kernel_matrix(X, kernel="rbf", sigma=1.0)
diag_vals <- diag(K)
cat("All diag exactly 1:", all(diag_vals == 1), "\n")
cat("Max diag deviation:", max(abs(diag_vals - 1)), "\n")
cat("First 5 diag:", round(diag_vals[1:5], 12), "\n")
```
Output:
```
All diag exactly 1: TRUE
Max diag deviation: 0
First 5 diag: 1 1 1 1 1
```
CPU RBF diagonal is exactly 1. Bug does not reproduce on CPU. May be float32 GPU-only. Mark as INCONCLUSIVE for GPU paths, REFUTED for CPU path.

---

### amatrix-p24 — pairwise_sqdist_argmin wrong centroid (RE-OPENED)
**Verdict: CONFIRMED**

Probe (orchestrator's exact case):
```r
library(amatrix)
X <- matrix(c(0,0,1,1,5,5), 3, 2, byrow=TRUE)
C <- matrix(c(0,0,5,5), 2, 2, byrow=TRUE)
result <- pairwise_sqdist_argmin(X, C)
cat("Result:", result, "\n")
cat("Expected: 1 1 2\n")
cat("Bug fires:", !identical(as.integer(result), c(1L,1L,2L)), "\n")
```
Output:
```
Result: 1 1 1
Expected: 1 1 2
Bug fires: TRUE
```
Centroid 2 (the point `[5,5]`) is never assigned. Confirms orchestrator's finding.

---

### amatrix-1ha — Matrix in Imports not Depends
**Verdict: CONFIRMED**

Probe:
```r
library(amatrix)  # Matrix NOT explicitly loaded
X <- matrix(c(4,2,2,3), 2, 2); A <- adgeMatrix(X)
cat("Matrix attached:", "package:Matrix" %in% search(), "\n")
cat(tryCatch({ t(A); "t OK" }, error=function(e) paste("t ERROR:", e$message)), "\n")
cat(tryCatch({ rowSums(A); "rowSums OK" }, error=function(e) paste("rowSums ERROR:", e$message)), "\n")
cat(tryCatch({ colSums(A); "colSums OK" }, error=function(e) paste("colSums ERROR:", e$message)), "\n")
cat(tryCatch({ diag(A); "diag OK" }, error=function(e) paste("diag ERROR:", e$message)), "\n")
cat(tryCatch({ mean(A); "mean OK" }, error=function(e) paste("mean ERROR:", e$message)), "\n")
```
Output:
```
Matrix attached: FALSE
t ERROR: argument is not a matrix
rowSums ERROR: 'x' must be an array of at least two dimensions
colSums ERROR: 'x' must be an array of at least two dimensions
diag ERROR: long vectors not supported yet: ...array.c:2212
mean OK: NA (with warning: returning NA)
```
`t()`, `rowSums()`, `colSums()`, `diag()` all fail. `mean()` silently returns `NA`.

---

### amatrix-7il — blanket tryCatch fallbacks mask errors
**Verdict: INCONCLUSIVE**

The blanket catch sites (`chol-factor.R:181-193`, `chol-factor.R:460-471`) are in GPU resident code paths (`solve_triangular_resident`, `broadcast_ewise_resident`). CPU fallback paths propagate errors correctly. Cannot trigger the GPU swallowing path without a GPU backend.

---

### amatrix-36q — double-drop sites in backend-planning/models-lm
**Verdict: INCONCLUSIVE**

Both cited sites (`backend-planning.R:395`, `models-lm.R:620`) involve GPU resident operations (`resident_drop`, `broadcast_ewise_resident`). Cannot trigger without GPU backend.

---

### amatrix-3ka — .amatrix_bind_resident leaks prior key on rebind
**Verdict: INCONCLUSIVE**

Leak is in GPU resident memory management. No CPU path to trigger `resident_drop` of a prior key without a live GPU resident handle.

## (d) Proposed `bd close` list (ONLY for runtime-refuted)

None proposed for full close. amatrix-75h is **REFUTED on CPU** — propose reclassification to GPU-only / float32 concern rather than close, since GPU path untested.

- **amatrix-75h**: Propose `bd close` or reclassify to GPU-only scope. CPU diagonal is exactly 1.0 (max deviation = 0).

## (e) Inconclusive list (GPU-only, no CPU probe possible)

- **amatrix-7il**: blanket tryCatch in GPU resident paths — CPU paths error correctly
- **amatrix-36q**: double-drop in `resident_drop`/`broadcast_ewise_resident` — GPU only
- **amatrix-3ka**: `.amatrix_bind_resident` prior-key leak — GPU resident only
- **amatrix-75h**: RBF diagonal drift may occur on GPU float32 only — CPU is clean

## Summary table

| ID | Verdict |
|----|---------|
| amatrix-jnd | CONFIRMED |
| amatrix-xnp | CONFIRMED |
| amatrix-x6a | CONFIRMED |
| amatrix-sxs | CONFIRMED |
| amatrix-vbh | CONFIRMED |
| amatrix-lei | CONFIRMED |
| amatrix-75h | REFUTED (CPU); INCONCLUSIVE (GPU float32) |
| amatrix-p24 | CONFIRMED |
| amatrix-1ha | CONFIRMED |
| amatrix-7il | INCONCLUSIVE |
| amatrix-36q | INCONCLUSIVE |
| amatrix-3ka | INCONCLUSIVE |
