# ALTREP Integration — Research Note and Deferral Rationale

## What ALTREP Is

ALTREP (Alternative Representation) is an R 3.5+ mechanism that allows a
package to register a custom representation for an R vector.  The key hook is
`DATAPTR` (and its read-only variant `DATAPTR_OR_NULL`): when R code accesses
the raw memory of a vector — as `matrix` arithmetic, `lm()`, `prcomp()`, etc.
all do — ALTREP intercepts the call and can materialise the data on demand.

In theory this means an `adgeMatrix` wrapped in an ALTREP vector could be
passed directly to *any* base-R function without pre-materialisation, and the
GPU→CPU copy would happen only if and when the function actually touches the
bytes.

## Why It Is Deferred

### 1. adgeMatrix is S4 with a 2-D dim — ALTREP is 1-D

ALTREP classes wrap a single R **vector** (1-D).  A matrix in R is a vector
with a `dim` attribute.  Emulating a 2-D matrix via ALTREP requires:

* registering a 1-D ALTREP class whose length is `nrow * ncol`,
* overriding `Attrib` and `SetAttrib` to synthesise the `dim` and `dimnames`
  attributes on the fly,
* maintaining a parallel S4 object so that our S4 generics (`%*%`, `crossprod`,
  …) still dispatch correctly.

That dual-representation bookkeeping is non-trivial and fragile across R
version updates.

### 2. DATAPTR materialises the entire matrix

The moment any C function calls `DATAPTR(x)` — not just `x[i,j]` in R — the
full `nrow × ncol × 8` bytes must be copied from GPU to CPU.  For the intended
use case of amatrix (large matrices resident on GPU for chain computations)
this is exactly the overhead we want to *avoid*.  ALTREP defers the copy but
cannot eliminate it when the callee is opaque C code.

### 3. S4 dispatch already covers the high-value paths

The S4 generics registered in `methods-dense.R` and `dispatch-hardening.R`
cover every operator that users are likely to call explicitly:
`%*%`, `crossprod`, `tcrossprod`, `solve`, `chol`, `qr`, `svd`, `eigen`,
`rowSums`, `colSums`, `Ops`, `[`, `[<-`, `t`, `diag`.

The gaps that remain — functions that call raw C internals without going
through an S4 generic — are exactly the cases where ALTREP's `DATAPTR` hook
would trigger a full materialisation anyway.  There is no net benefit over the
current explicit materialisation in the fallback path.

### 4. ALTREP API stability risk

The ALTREP C API (`R_altrep_class_t`, `R_make_altreal_class`, `ALTREP_DATA1`,
`ALTREP_DATA2`) has changed between minor R releases and is not part of the
stable API guaranteed by Writing R Extensions.  Maintaining a C-level ALTREP
shim across R 4.x releases would add ongoing maintenance burden with no clear
user-visible benefit given point 3 above.

## When to Revisit

ALTREP becomes worthwhile if amatrix ever needs to support:

* **Read-only lazy slicing** — exposing a row/column slice of a GPU-resident
  matrix to base-R code without copying the full matrix.  ALTREP's
  `DATAPTR_OR_NULL` returning `NULL` (for read-only access refusal) plus a
  custom `Extract_subset` method could make this zero-copy for sequential
  access patterns.

* **Integration with packages that accept arbitrary numeric vectors but go
  through `[` rather than raw DATAPTR** — e.g., some tidyverse or data.table
  paths.  But those packages typically have S3/S4 hooks that are cheaper to
  target directly.

Until one of those use-cases materialises with a concrete benchmark showing
ALTREP would help, the S4 dispatch + explicit materialisation fallback is the
right architecture.

## References

* Luke Tierney — *ALTREP: Alternate Representation of Basic R Objects*
  (useR! 2018 keynote)
* `src/include/R_ext/Altrep.h` in the R source tree
* `?ALTREP` in package **altrep** (CRAN, illustrative)
* Writing R Extensions §5.15 "Registering native routines" — note ALTREP is
  *not* in this stable API section.
