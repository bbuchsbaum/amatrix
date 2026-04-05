# MLX Compact QR Plan

## Purpose

This note fixes the next QR milestone for `amatrix.mlx`.

The current MLX QR path is useful, but it is not the end state:

- factorization calls `mlx_linalg_qr()` and gets explicit `Q` and `R`
- `Q` stays resident
- helper paths such as `qr.coef()` and `qr.qty()` can reuse that resident `Q`
- this is already good enough to win on shared-`X`, many-`Y` workloads at larger RHS counts

It is still not a true compact QR factor path. The next step is to replace the explicit-`Q` hot path with a compact factor representation for `pivot = FALSE`.

## Why Research First

This is not a wrapper task.

The installed MLX C surface exposes `mlx_linalg_qr()`, but it does not expose a LAPACK-style compact QR factor API such as `geqrf` with Householder vectors and `tau`. That means a true compact QR path requires a custom algorithm on top of MLX primitives and, where needed, custom MLX extensions and Metal kernels.

The algorithm choice matters because the workload target is not generic `qr(A)` in isolation. The target is:

- `qr.coef()`
- `qr.qty()`
- repeated many-RHS solves
- `many_lm(..., method = "qr")`

Those workloads do not need an explicit full `Q` in the hot path.

## Current Constraints

### MLX API reality

- MLX C exposes `mlx_linalg_qr()` and the rest of the dense linalg family.
- MLX C also exposes compilation and custom extension hooks.
- MLX documentation still states that `float64` is CPU-only, so accelerated QR must remain a `fast`-mode target first.

### Current repo reality

- `many_lm(..., method = "qr")` is now the flagship workflow.
- The current resident-`Q` MLX path already wins once the number of right-hand sides is large enough.
- The remaining architectural weakness is dependence on explicit `Q` for the hot helper path.

## Algorithm Choice

### Public contract

Keep:

- `qr(x)` as the public S4 entry point
- `am_qr(x, ...)` as the internal wrapper
- `amQR` as the factorization object

Add a true MLX compact representation under `amQR`, rather than exposing it as a separate public class.

### Fast-path scope

The first accelerated compact QR target is:

- dense
- `pivot = FALSE`
- `thin = TRUE`
- `mode = "fast"`

Fallback rules remain:

- `exact` / `strict`: CPU QR
- `pivot = TRUE`: CPU or hybrid fallback

### Shape-specialized plan

1. Tall-skinny matrices: TSQR / CAQR-style reduction
- Primary target for regression design matrices and the current flagship workload.
- Store reflector metadata for leaves and merges rather than explicit `Q`.

2. General dense matrices: blocked Householder QR with compact WY form
- Use blocked panel factorization and GEMM-heavy trailing updates.
- This matches where MLX is already strongest.

3. Many small matrices: batched QR later
- Important, but not the first milestone.

## Internal Representation

The first true compact MLX QR payload should carry:

- `representation = "mlx_compact_qr"`
- backend-pinned resident state
- source dimensions
- effective rank
- thin/pivot metadata
- method metadata
- compact factor buffers for reflector data and update metadata

Helper methods should operate directly on that factor:

- `qr.qty()`
- `qr.qy()`
- `qr.coef()`
- `qr.solve()`
- `qr.resid()`

Only `qr.Q()` should force explicit `Q` materialization.

## Performance Rules

1. No repeated host copies in the hot path
- Factor once
- keep factor state resident
- reuse for many RHS

2. Backend pinning once resident
- No cross-accelerator hopping for a live factor object
- CPU fallback only when the pinned backend cannot serve the requested helper

3. Compile and cache by stable signature
- cache by backend, method, shape, dtype, and pivot mode
- distinguish cold compile/setup from hot steady-state timings

4. Benchmark evaluated work only
- force evaluation
- synchronize the stream
- measure cold and hot paths separately

## Execution Stages

### Stage 1

Freeze the compact-factor contract:

- add `mlx_compact_qr` representation under `amQR`
- route helper methods through representation-specific internals
- keep current resident-`Q` path as the control

### Stage 2

Implement a first tall-skinny compact path for:

- `pivot = FALSE`
- `thin = TRUE`
- `fast`

This is the first serious attempt to beat the resident-`Q` path on the flagship workload.

### Stage 3

Add blocked Householder / compact WY for general dense matrices.

### Stage 4

Use custom MLX extensions and Metal kernels only where the stock primitive path is not enough:

- panel kernels
- reflector application
- fused helper paths where profiling justifies it

## Success Criteria

The compact QR project is successful when:

1. `am_qr_info()` reports a native compact MLX factor source rather than `bridge_compact`.
2. `many_lm(..., method = "qr")` on MLX no longer depends on explicit `Q` in the hot path.
3. Many-RHS QR workflows beat both:
   - cached base QR
   - the current resident-`Q` explicit MLX path
4. `exact` mode remains honest and falls back cleanly.

## Sources

- [MLX C linalg surface](https://ml-explore.github.io/mlx-c/build/html/ops.html)
- [MLX custom extensions](https://ml-explore.github.io/mlx/build/html/dev/extensions.html)
- [MLX data types](https://ml-explore.github.io/mlx/build/html/python/data_types.html)
- [MLX unified memory](https://ml-explore.github.io/mlx/build/html/usage/unified_memory.html)
- [Communication-avoiding QR factorizations](https://digicoll.lib.berkeley.edu/record/136322/files/EECS-2008-74.pdf)
- [WY representation for Householder products](https://ecommons.cornell.edu/items/92a11030-dca1-45d4-a0ba-732cf962b2b2)
