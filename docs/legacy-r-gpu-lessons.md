# Lessons From The Old R GPU Ecosystem

This note records a design constraint for `amatrix`.

The older R GPU package landscape is useful mainly as a warning document:

- many packages wrapped isolated functions such as `gpuLm()`, `gpuCor()`, or `gpuDist()`
- acceleration claims were often package-level rather than operation- and workload-level
- performance depended heavily on matrix size and data-movement overhead
- package users often had little visibility into what was actually accelerated

`amatrix` should not repeat that pattern.

## Core Interpretation

The project should be infrastructure, not a zoo of GPU-flavored wrappers.

That means:

- one object and kernel substrate
- optional backend engines
- factor-first runtime objects
- a small number of benchmark-backed flagship workflows

The goal is not to accumulate `gpuFoo()` entry points for every statistical primitive. The goal is to make many algorithms ride on the same matrix, factor, residency, and backend-planning layer.

## Design Rules

### 1. One substrate, not many wrappers

Prefer:

- `adgeMatrix`
- backend-neutral `am_*` kernels
- factor objects such as `amQR`, `amChol`, and `amSVD`
- model-core helpers such as `many_lm()`

Avoid:

- one-off public wrappers whose main story is “this is the GPU version of function X”

### 2. Capability must be empirical

A backend exposing an API is not enough.

For each public operation, the docs and tests should make clear whether it is:

- accelerated and validated
- supported but CPU fallback
- experimental
- unavailable

The package should describe what is true on the current product path, not what might be possible because a vendor library has a primitive.

### 3. Workflow benchmarks matter more than kernel bragging

Kernel benchmarks still matter, but they are not the product story.

The main benchmark surfaces should be workloads users actually feel:

- shared-`X`, many-`Y` QR regression
- SPD factor plus batched solve
- randomized SVD
- covariance and correlation blocks

### 4. Keep a low-level escape hatch

Package authors should be able to adopt `amatrix` without waiting for a dedicated high-level wrapper.

That is why the factor-first and kernel-first APIs matter:

- `am_qr()`
- `am_chol_factor()`
- `am_chol_solve()`
- `am_svd_factor()`
- `am_matmul()`, `am_crossprod()`, `am_solve()`

### 5. Stay honest about workload shape

`amatrix` should not promise “GPU acceleration everywhere.”

It should promise something more useful:

- when the computation has the right shape, the package will take the fast path
- when it does not, the package will stay boring and correct on CPU

Dense, repeated, high-arithmetic-intensity workloads are the wedge. Sparse remains conservative until there is a validated sparse story.

## Product Consequence

The durable `amatrix` shape is:

- Matrix-compatible object layer
- backend registry and capability planning
- resident arrays and cached factors
- a few obvious hero workflows

If the package stays centered there, it avoids the main failure mode of the old R GPU ecosystem: demos and wrappers without a reusable runtime underneath.
