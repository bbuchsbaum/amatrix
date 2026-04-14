# ALTREP Integration — Adoption Multiplier, Not The Primary Speed Engine

## What ALTREP Can Actually Buy

ALTREP matters in `amatrix` for one reason: downstream code that defensively
calls `as.matrix(x)` or otherwise asks for plain R vector memory currently
collapses residency and forces host materialization. An ALTREP-backed dense
payload could defer that copy until bytes are actually required, which makes
near-zero-change adoption more real for Matrix-aware third-party code.

That makes ALTREP an adoption multiplier. It is not the main benchmark lever.

## What ALTREP Does Not Solve

### 1. It does not replace kernel or factor work

`amatrix` wins against CPU or `torch` on this roadmap through shared-`X`
many-`Y` QR, Cholesky plus batched solve, cached factorizations, and explicit
high-value kernels. ALTREP does not create those wins by itself.

### 2. It does not make opaque C consumers GPU-native

The moment an external C path calls `DATAPTR(x)`, the full dense payload still
has to exist on host. ALTREP can defer the copy, but it cannot eliminate it for
opaque LAPACK- or C-level consumers.

### 3. It does not solve common algebraic syntax by itself

The right fix for `t(A) %*% B` and `A %*% t(B)` is a Matrix-compatible
transpose view inside `amatrix` dispatch. The right fix for full `alpha`,
`beta`, and accumulator control is an explicit `am_gemm()` surface. ALTREP is
orthogonal to both.

## Current Stepping Stone

The current code has a useful but incomplete shortcut:

* `t()` on a resident `adgeMatrix` returns a new `adgeMatrix` carrying `src_id`
  so `%*%` can route to `crossprod_resident()` or `tcrossprod_resident()`
  without re-uploading the source matrix.
* That closes a real performance hole for GPU-resident transpose products.
* It is still not a true structural view. The transposed result keeps material
  host data in `@x`, so it still pays `O(nm)` host transpose work and memory.

That shortcut is a valid stepping stone. The cleaner next step is a dedicated
`aTransposeView`-style class that behaves like a matrix without pretending to be
fully materialized dense host storage.

## Why ALTREP Is Still Deferred

### 1. The surface contract should stabilize first

The current roadmap is to keep public materialized values boring, allow
specialized S4 structural or factor intermediates where they buy real speed,
and make the transpose and GEMM boundary explicit before adding another
representation layer.

### 2. adgeMatrix is S4 with a 2-D dim, while ALTREP is fundamentally 1-D

ALTREP classes wrap a vector. A matrix in R is a vector plus `dim` and
`dimnames`. Bridging that into `amatrix` means:

* registering a 1-D ALTREP payload of length `nrow * ncol`,
* synthesizing `dim` and `dimnames` correctly,
* keeping S4 dispatch coherent for `%*%`, `crossprod`, `qr`, and friends,
* and making serialization and fallback deterministic.

That is possible, but it is a large surface to stabilize.

### 3. The first transparent win should be a proper transpose view, not ALTREP

The highest-value syntax gap is still transpose-heavy algebra. The current
`src_id` shortcut is a good interim fix, but a dedicated transpose-view class
would remove the remaining host transpose cost and make later ALTREP work
cleaner.

### 4. ALTREP still carries maintenance and API-risk cost

The ALTREP C API is not the most stable part of the R extension surface.
Maintaining a dense ALTREP shim across R releases is reasonable only once the
expected user benefit is concrete enough to justify it.

## Relationship To The Main Roadmap

The intended sequence is:

1. Keep the current Matrix-compatible S4 contract explicit.
2. Replace the current `src_id` transpose shortcut with a proper structural
   transpose view.
3. Expose explicit kernels such as `am_gemm()` for full BLAS-style control.
4. Revisit ALTREP once the flagship workflows and structural-view layer are
   stable.

That sequencing keeps ALTREP in the right role: it broadens transparent
adoption, but it is not the primary speed story.

## When To Revisit

ALTREP becomes worth serious implementation effort if:

* downstream packages that should be easy adoption targets keep defeating
  residency by calling `as.matrix()` defensively,
* a benchmark shows that deferred host realization would materially improve a
  real adoption path,
* and the transpose-view plus explicit-kernel surface is already stable.

## References

* Luke Tierney — *ALTREP: Alternate Representation of Basic R Objects*
  (useR! 2018 keynote)
* `src/include/R_ext/Altrep.h` in the R source tree
* `?ALTREP` in package **altrep** (CRAN, illustrative)
