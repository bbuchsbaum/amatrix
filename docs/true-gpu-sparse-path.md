# True GPU Sparse Kernel Path

## Current State

`amatrix` already has a backend-spanning sparse-dense product pathway, but the
Apple `mlx` implementation is not a true GPU sparse kernel today.

Local facts in this repo:

- [`backends/amatrix.mlx/configure`](../backends/amatrix.mlx/configure) links
  against `mlx-c`, not the C++ extension/tooling surface.
- [`backends/amatrix.mlx/src/init.c`](../backends/amatrix.mlx/src/init.c)
  registers only C bridge entry points.
- [`backends/amatrix.mlx/src/amatrix_mlx_matmul.c`](../backends/amatrix.mlx/src/amatrix_mlx_matmul.c)
  implements sparse `%*%` as CPU-side CSC loops, even when the backend name is
  `mlx`.
- Focused sparse benchmarks show the current MLX sparse path still loses to the
  CPU baseline, so sparse MLX routing is disabled by default.

That means the current MLX sparse story is:

- useful for validating backend seams and resident sparse reuse
- not yet a real Apple GPU sparse acceleration path

## What The Upstream Stack Allows

Official MLX docs currently expose:

- the general MLX runtime and C++ API:
  https://ml-explore.github.io/mlx/build/html/index.html
- custom extensions in MLX:
  https://ml-explore.github.io/mlx/build/html/dev/custom_extensions.html
- custom Metal kernels:
  https://ml-explore.github.io/mlx/build/html/dev/custom_metal_kernels.html

That is promising, but it matters where those hooks live:

- MLX documents custom Metal kernels and extension primitives on the MLX
  Python/C++ side.
- `amatrix.mlx` currently integrates through `mlx-c`, not through MLX C++
  extension primitives.
- So a true GPU sparse path is not a small patch to the current bridge. It
  requires either:
  1. going below MLX and calling Metal directly
  2. adding a new C++ extension layer around MLX

Apple's documented MPS matrix layer is dense-oriented, for example:

- `MPSMatrixBinaryKernel`:
  https://developer.apple.com/documentation/metalperformanceshaders/mpsmatrixbinarykernel

So there is no obvious Apple-supplied sparse matrix multiply API we can simply
bind and call as a drop-in replacement.

## The Real Contract Gap

The current sparse resident hook is:

```r
spmm_resident = function(sp_key, B, trans_lhs = FALSE) NULL
```

That returns a host matrix. It is fine for opportunistic sparse support, but it
caps the upside of a true GPU kernel because:

- dense RHS may still arrive from host each call
- results come back to host immediately
- repeated sparse-dense products cannot stay fully resident on device

For a real GPU path, the missing backend seam is an optional resident-to-resident
sparse multiply surface:

```r
spmm_resident_key = function(sp_key, y_key, out_key,
                             trans_lhs = FALSE,
                             defer = FALSE) NULL

spmm_resident_host = function(sp_key, y, trans_lhs = FALSE) NULL
```

`spmm_resident_host` is the compatibility fallback.
`spmm_resident_key` is the performant surface.

Without `spmm_resident_key`, a true GPU sparse kernel will still be throttled by
host round-trips.

## Apple-Specific Options

### Option A: Direct Metal Sparse Kernel

This is the most credible path for maximum Apple performance.

Shape:

- implement the sparse kernel in ObjC++ or C++ plus Metal
- keep the current R-side backend API
- add Metal-backed sparse resident storage and sparse-dense resident multiply

Recommended resident storage:

- convert public CSC input to CSR once at `sparse_resident_store()`
- optionally keep CSC too for transpose-heavy cases
- store float32 values for `fast` mode first

Why CSR first:

- `X %*% B` is row-oriented at execution time
- CSR maps naturally to row-parallel GPU work
- a row-tiled dense RHS kernel is easier to write and tune

Recommended first kernel:

- `SpMM`: sparse CSR `[m, k]` times dense resident RHS `[k, n]`
- target `n >= 8` many-RHS cases first
- output a dense resident matrix

Do not start with:

- `SpGEMM`
- sparse QR/Cholesky/solve
- host-returning kernels only

### Option B: MLX Extension Primitive Plus Custom Metal Kernel

This is architecturally cleaner if the goal is to stay inside MLX graphs.

Shape:

- add an MLX C++ extension primitive for sparse-dense multiply
- implement `eval_gpu` with a Metal kernel
- expose that primitive through a new native bridge layer

Advantages:

- better alignment with MLX lazy execution
- potential reuse of MLX stream/device infrastructure
- better long-term fit if more custom MLX kernels are expected

Costs:

- larger build-system change than the current `mlx-c` bridge
- requires C++ extension plumbing, not just C wrappers
- higher initial integration risk for an R package

## Recommendation

If the goal is maximum Apple sparse performance, prefer Option A first:

1. keep MLX for dense kernels
2. add a direct Metal sparse resident kernel path for `SpMM`
3. bridge that path through the existing `amatrix` backend contract

Reason:

- it targets the real bottleneck directly
- it avoids waiting on MLX sparse abstractions that do not currently exist in
  the `mlx-c` path
- it gives a clean benchmark-backed answer quickly

If Option A wins convincingly, the project can later decide whether to:

- keep it as an Apple-specific sparse fast path
- or fold the same kernel idea into a richer MLX extension layer

## Minimum Viable Implementation

Phase 1:

- add optional `spmm_resident_key`
- keep existing `spmm_resident` as host fallback
- update wrappers to use resident dense RHS when available

Phase 2:

- add Metal sparse resident store:
  - CSC input from R
  - resident CSR buffers for forward `SpMM`
  - optional CSC twin for transpose cases

Phase 3:

- implement `SpMM` kernel for:
  - `adgCMatrix %*% adgeMatrix`
  - `crossprod(adgCMatrix, dense)`
  - `tcrossprod(adgCMatrix, dense)`

Phase 4:

- benchmark crossover by:
  - `nnz`
  - density bucket
  - RHS width
  - cold versus resident

Success criterion:

- resident sparse `SpMM` must beat CPU `Matrix` on a meaningful many-RHS slice
  before it is enabled by default

## What Not To Do

- Do not re-enable MLX sparse routing by default without winning benchmarks.
- Do not start with `SpMV` if `SpMM` still loses.
- Do not expand to `SpGEMM` before a resident sparse-dense GPU path wins.
- Do not tie public sparse semantics to an experimental backend path.

## Immediate Next Step

The next worthwhile spike is:

1. extend the backend contract with optional `spmm_resident_key`
2. wire wrappers to use it when available
3. prototype a direct Metal `SpMM` kernel for resident sparse CSR times resident
   dense RHS
4. benchmark only the resident many-RHS cases first

That is the shortest path from today's CPU sparse bridge to a real GPU sparse
kernel on Apple hardware.
