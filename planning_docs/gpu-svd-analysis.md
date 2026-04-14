# GPU SVD Analysis: Why Naive Wrappers Fail and What Actually Works

## Current Implementation Status (2026-04)

| Algorithm | Status | Backend |
|-----------|--------|---------|
| `am_rsvd` (GPU-native randomized SVD) | **Implemented** | MLX C bridge |
| `am_irlba` (Lanczos wrapper, GPU matvecs) | Implemented | MLX + ArrayFire (slow — see below) |
| `am_irlba_native` (Lanczos in ArrayFire C) | Implemented | ArrayFire |
| `am_block_lanczos` (GEMV→GEMM upgrade) | **Implemented prototype** | MLX via current matmul/crossprod kernels |
| `am_rsvd` on ArrayFire | Planned | ArrayFire |

`am_rsvd` via the MLX C bridge is the correct GPU decomposition path. Use it for performance.
`am_irlba` is a compatibility wrapper — correct, not fast.

---

## Summary

Wrapping irlba or svdr with GPU-dispatched `%*%` methods still does not beat native CPU LAPACK at moderate matrix sizes. That part of the analysis remains true.

What changed is the product path: `amatrix` now has a device-native MLX randomized SVD bridge, and `am_svd_factor()` can route to it automatically for fast, low-rank dense factorizations. On Apple Silicon this path has a one-time MLX warm-up cost, but once warmed it can beat CPU `svdr` on medium and large matrices. The current default policy is conservative and keeps small problems on CPU.

Historical benchmark tables below describe the wrapper approaches and why they fail. They do not describe the current C-level MLX `am_rsvd` path.

## Current Product State

- `am_irlba` is still a compatibility path. It is not the performance path on Apple Silicon.
- `am_block_lanczos()` is now implemented as an explicit GEMM-oriented prototype in `R/irlba.R`.
- The current block defaults are speed-oriented, not strict replacements for CPU `irlba`: on representative dense MLX cases they currently trade about 4.2-4.8% relative singular-value error for materially lower wall time than `am_irlba()`. After switching the right-looking step to a resident MLX transpose, routing both left and right block products through direct resident matmul helpers, reusing the per-step left reorthogonalization coefficients to build the projected small matrix without a full extra `Q_L^T A Q_R` pass, and retuning the default block heuristic so smaller ranks stay mildly oversampled while the current `k = 20` path uses `block_size = k`, the current tuned settings now beat CPU `irlba::irlba()` at both `3000x1200` and `5000x2000` on this machine.
- `am_rsvd` is implemented in `backends/amatrix.mlx/src/amatrix_mlx_matmul.c`.
- Standalone `am_rsvd()` now dispatches to a backend path only in `fast` precision; under strict precision it falls back to the CPU reference path.
- `am_svd_factor()` now supports `method = "auto" | "exact" | "rsvd"`.
- In `auto` mode, `am_svd_factor()` selects `rsvd` only when:
  - precision is `fast`
  - a backend with `rsvd` support is available
  - `min(dim(x)) >= 400`
  - `k / min(dim(x)) <= 0.20`
- On Apple Silicon, the current implementation will prefer MLX for truncated SVD even when the input object was requested with `preferred_backend = "arrayfire"` and ArrayFire itself lacks an `rsvd` kernel.

Use `tools/benchmark-svd-factor.R` for CPU/reference numbers. For MLX steady-state timing on Apple Silicon, run direct `Rscript -e '...'` commands from the shell; nested or file-entry MLX benchmark harnesses are currently unstable on this machine. `tools/print-svd-factor-calibration.R` prints a parameter grid of those direct commands for the current crossover sizes.

For the block Lanczos prototype, use:

```sh
Rscript -e 'source("tools/benchmark-block-lanczos.R", local = TRUE); print_block_lanczos_benchmark()'
```

The benchmark script is intentionally source-driven because direct file-entry MLX benchmark runs are unstable on this machine.

For a stage breakdown of the current implementation, use:

```sh
Rscript -e 'source("tools/profile-block-lanczos.R", local = TRUE); print_block_lanczos_profile(n = 3000L, p = 1200L, k = 20L, block_size = 20L, n_steps = 4L)'
```

Current source-driven block Lanczos benchmarks on this machine:

| Size | CPU `irlba::irlba` | `am_irlba` | `am_block_lanczos` | `max_rel_sv_err` |
|------|--------------------|------------|---------------------|------------------|
| 1200×600, `k=10`, default `block_size=11` | 0.008s | 0.205s | 0.015s | 0.0423 |
| 3000×1200, `k=20`, default `block_size=20` | 0.123s | 0.479s | 0.071s | 0.0479 |
| 5000×2000, `k=20`, default `block_size=20` | 0.400s | 0.713s | 0.114s | 0.0462 |

Recent spot checks used to calibrate the cutoff:

| Size | CPU `svdr` | MLX `rsvd` | `auto` method | Outcome |
|------|------------|------------|---------------|---------|
| 300×240 | 0.014s | 0.029s | exact | CPU wins |
| 400×320 | 0.021s | 0.031s | exact | CPU wins |
| 500×400 | 0.031s | 0.031s | rsvd | parity |
| 700×560 | 0.057s | 0.031s | rsvd | MLX wins |
| 1000×800 | 0.108s | 0.031s | rsvd | MLX wins |

The default cutoff therefore now starts `auto` at `min(dim) >= 400`: it keeps the known losing small cases on CPU while allowing `500x400` and above onto the faster low-rank path.

Example steady-state MLX spot benchmark:

```sh
Rscript -e 'orig <- .libPaths(); .libPaths(c("/tmp/amatrix-lib", orig)); pkgload::load_all(".", quiet = TRUE); invisible(loadNamespace("amatrix.mlx")); options(amatrix.mlx.available = TRUE); n <- 1000L; p <- 800L; k <- 20L; n_oversamples <- 10L; n_iter <- 2L; set.seed(20260405L + n + p + k); host <- matrix(rnorm(n * p), nrow = n, ncol = p); invisible(amatrix.mlx:::amatrix_mlx_rsvd(matrix(rnorm(32L * 16L), nrow = 32L, ncol = 16L), k = 5L, n_oversamples = 4L, n_iter = 1L)); t_cpu <- system.time(cpu <- irlba::svdr(host, k = k, extra = n_oversamples, it = n_iter)); t_auto <- system.time(fac_auto <- am_svd_factor(adgeMatrix(host, preferred_backend = "mlx", precision = "fast"), k = k, method = "auto", n_oversamples = n_oversamples, n_iter = n_iter)); ref <- base::svd(host, nu = k, nv = k)$d[seq_len(k)]; cat(sprintf("svdr=%.3f\nauto=%.3f\nauto_method=%s\nauto_err=%.4f\n", unname(t_cpu[["elapsed"]]), unname(t_auto[["elapsed"]]), amatrix:::.amatrix_svd_factor_plan(adgeMatrix(host, preferred_backend = "mlx", precision = "fast"), k, "auto", n_oversamples, n_iter)$method, max(abs(fac_auto@d - ref) / pmax(abs(ref), 1e-12))))'
```

## Benchmarks

Machine: Apple Silicon (M-series), MLX backend.

### am_irlba vs base irlba (k=20 singular values)

| Size       | Base (C fastpath) | Base (R path) | am_irlba (MLX) | vs C fastpath |
|------------|-------------------|---------------|----------------|---------------|
| 2000×1000  | —                 | 1.16s         | 1.73s          | —             |
| 5000×2000  | 0.93s             | 6.1s          | 9.0s           | 9.7× slower   |
| 10000×3000 | —                 | 24.6s         | 32.1s          | —             |

irlba with `fastpath=FALSE` performs ~608 matrix-vector products (GEMVs) for a 5000×2000
problem with k=20. The MLX path does 686 — 13% more — because float32 arithmetic degrades
Lanczos orthogonality, forcing additional restarts. Each GEMV is ~20M flops — low arithmetic
intensity, memory-bound.

After fixing the A re-upload bug (vector promotion in am_matmul), A is uploaded once and kept
resident. Singular values match to 1.47e-8. But the per-step sync overhead (~13ms vs ~1.4ms
CPU BLAS) and extra f32 restarts together keep wall time at ~9.7× the C fastpath — worse than
the pre-fix 8.2s because the extra restarts from f32 now dominate.

### Historical svdr Wrapper Benchmark (k=20, 100 iterations)

| Size       | Base (CPU) | MLX    | Speedup |
|------------|-----------|--------|---------|
| 2000×1000  | 6.6s      | 7.4s   | 0.89×   |
| 5000×2000  | 43.6s     | 46.2s  | 0.94×   |
| 10000×3000 | 152.8s    | 155.3s | 0.98×   |

The speedup converges toward 1.0× as matrix size grows — near-parity is the ceiling for the
wrapper approach, not a floor. As GEMM size increases, the QR materialization penalty scales
proportionally, so the round-trip overhead never becomes negligible.

svdr uses GEMMs (A %*% Q where Q has k+extra columns) which have higher arithmetic intensity. But this wrapper path still loses.

## Root Cause: The QR Round-Trip

Every iteration of both algorithms forces a CPU-GPU round-trip at the orthogonalization step:

```
GPU: x %*% Q          → A %*% Q GEMM on device
     ↓ device→host    ← materialize result (qr() forces as.matrix())
CPU: qr.Q(qr(Z))       ← QR factorization on host
     ↓ host→device    ← upload Q for next step
GPU: crossprod(x, Q)  → t(A) %*% Q GEMM on device
     ↓ device→host    ← materialize again
CPU: qr.Q(qr(Z))       ← QR on host
... × 100 iterations
```

For a 5000×2000 matrix with k+extra=30 work vectors:
- Result Z of each GEMM: 5000×30 = 1.2MB
- Q uploaded next step: 5000×30 = 1.2MB
- Per iteration transfer: ~3MB
- 100 iterations: 300MB total transfers

At Apple Silicon memory bandwidth, this is ~6ms of transfer per iteration — comparable to the GEMM compute time. The GPU never gets a chance to amortize its launch overhead.

MLX's documentation explicitly warns against evaluating inside a tight loop. Each `mlx_array_eval()` call (triggered by `as.matrix()` / host materialization) has fixed overhead and forces a CPU-GPU synchronization.

## Why IRLBA Is Especially Hard

irlba's algorithmic contract is built around single matrix-vector products (GEMV):

- One vector per Lanczos step: `A %*% v` (m×n %*% n×1)
- Arithmetic intensity: ~0.25 flop/byte — purely memory-bound
- 608 sequential GEMVs for a 5000×2000 problem with k=20
- Each step depends on the previous (no batching)
- irlba's own vignette documents this: "only uses matrix vector products"

A GPU thrives on high arithmetic intensity (≥4 flop/byte). Single-vector GEMVs at these sizes are CPU-optimal workloads.

## What Actually Works: Fully Device-Native RSVD

The correct GPU decomposition avoids per-step host materialization entirely. A randomized SVD expressed as device-native operations:

```
Algorithm: GPU-native randomized SVD (Halko et al. 2011)

1. Ω = randn(n, k+p)              — on device (mlx_random_normal)
2. Y = A %*% Ω                    — GEMM on device (n×k+p columns → high intensity)
3. Q = qr(Y).Q                    — MLX QR on device (mlx_linalg_qr)

Power iteration (q passes, optional):
4.  Z = A^T %*% Q                 — GEMM on device
    Q = qr(Z).Q                   — MLX QR on device
    Z = A %*% Q                   — GEMM on device
    Q = qr(Z).Q                   — MLX QR on device

5. B = Q^T %*% A                  — GEMM on device (k×m, small)
6. Û, s, V = svd(B)               — MLX SVD on device (k×k, tiny)
7. U = Q %*% Û                    — GEMM on device

mlx_array_eval(U, s, V)           — SINGLE eval(), one device→host copy
```

Key properties:
- Zero per-step CPU-GPU syncs (all QR and matmul on device)
- Single `mlx_array_eval()` at the boundary
- MLX lazy evaluation composes all ops into a single fused graph
- Arithmetic intensity: GEMM with k+p columns — GPU-friendly at any matrix size

MLX C API provides all required primitives: `mlx_linalg_qr`, `mlx_linalg_svd`, `mlx_matmul`, `mlx_random_normal`.

## Implementation Status: am_rsvd

`am_rsvd` is now implemented as a C-level MLX bridge. The current bridge uses device-native sketching and power-iteration matmuls, Cholesky QR for device-friendly orthogonalization, and a small eigendecomposition on `B %*% t(B)` instead of a full SVD of `B`.

Observed product profile:
- Arithmetic intensity ≫ 1 flop/byte (GEMM with k+p columns)
- CPU-GPU synchronization only at explicitly materialized small intermediates and final outputs
- One-time MLX compile/warm-up cost on the first call
- Steady-state speedups versus CPU `svdr` at medium and large sizes on this machine; use direct `Rscript -e` spot benchmarks because file-entry MLX harnesses are unstable here

## What am_irlba Is Good For

`am_irlba` (R/irlba.R) is a valid convenience wrapper: it ensures `A` stays resident on device and routes `%*%` dispatch through the GPU. For very large sparse matrices where `A` does not fit in CPU cache, and for users who already have irlba code, it avoids copying `A` off device. But it will not beat `irlba(A_cpu, fastpath=TRUE)` at any realistic size.

Use `am_irlba` for code compatibility. Use `am_svd_factor(..., method = "auto")` or `am_rsvd()` for performance.

## Priority Order for Decompositions

1. ~~**`am_rsvd`**~~ — **Done.** GPU-native randomized SVD via MLX C bridge. Beats LAPACK rsvd at m,n≥3000 with k≤100.

2. ~~**Block Lanczos**~~ — **Prototype implemented.** `am_block_lanczos()` now provides a GEMM-oriented truncated-SVD surface in `R/irlba.R`. The current tuned settings route the right-looking step through a resident MLX transpose, use direct resident matmul helpers on both sides of the block iteration, and skip the final global QR when the concatenated block bases are already orthonormal within tolerance. That moves the MLX prototype to about `12.2x`, `5.4x`, and `4.9x` faster than `am_irlba()` at `1200x600`, `3000x1200`, and `5000x2000`, respectively, with CPU wins at both `3000x1200` (`0.091s` vs `0.126s`) and `5000x2000` (`0.147s` vs `0.421s`). The prototype still carries about `4.6-4.8%` relative singular-value error, so it remains explicitly approximate rather than a strict replacement for CPU `irlba`.

3. **`am_rsvd` on ArrayFire** — Port the MLX Halko algorithm to ArrayFire. All required primitives (GEMM, QR, small SVD) are available in the ArrayFire backend. Completes cross-platform rsvd parity.

4. **`am_irlba`** — CPU-compatible reference. Correct, not fast. Keep for code compatibility.

---

## Block Lanczos Design

The fundamental problem with `am_irlba` is that it issues k sequential GEMVs per restart, each with arithmetic intensity ~0.25 flop/byte. Block Lanczos upgrades this to batched GEMMs.

### Algorithm Sketch

```
Block Lanczos with block size b (e.g., b=16):

Input: A (m×n), k desired singular values
Block size b: process b Lanczos vectors per step

1. Q_0 = randn(n, b)                     # b starting vectors (on device)
2. For j = 1, 2, ..., ceil(k/b) + extra:
   Z = A %*% Q_{j-1}                     # GEMM: m×n * n×b = m×b  (single dispatch)
   Q_j = qr(Z).Q                         # MLX QR on device: m×b
   W = A^T %*% Q_j                       # GEMM: n×m * m×b = n×b
   Q_j_right = qr(W).Q                   # MLX QR on device: n×b

3. Collect basis: [Q_1, Q_2, ..., Q_J]  # n × (J*b) tall thin matrix
4. B = Q^T %*% A %*% Q_right            # small (J*b)×(J*b) bidiagonal form
5. SVD of B                              # tiny, on device
6. U = Q %*% U_B  ;  V = Q_right %*% V_B

mlx_array_eval(U, s, V)                 # SINGLE eval, one device→host copy
```

Key improvements over `am_irlba`:
- Each step issues **one GEMM** (`A %*% Q_{j-1}` where Q has b columns) vs. b sequential GEMVs
- Arithmetic intensity: GEMM with b columns = b× higher than single-vector GEMV
- MLX lazy evaluation fuses the entire block into a single graph before `mlx_array_eval()`
- Reorthogonalization is a batched GEMM against existing basis, not k sequential dot products
- Expected speedup: 5–20x over `am_irlba` at m,n≥3000

### Current Implementation

`am_block_lanczos()` is currently implemented in `R/irlba.R` on top of the existing dense-kernel surface:

- `A %*% Q` routes through the resident backend GEMM path
- on MLX fast mode, the right-looking `t(A) %*% Q` step now uses a resident-transposed MLX operator built once per solve, and both left/right block products use direct resident `matmul` helpers instead of the higher-overhead wrapped `%*%` result path
- the projected small matrix now reuses the same left reorthogonalization coefficients already computed inside the iteration, so only the final right block needs an extra `A %*% Q` pass
- the final global QR on concatenated block bases is skipped when the accumulated basis is already orthonormal within tolerance
- per-step block reorthogonalization uses two-pass CGS against the accumulated basis
- small block QR factorizations and the projected SVD remain on CPU

This is good enough to beat the compatibility Lanczos wrapper on the measured MLX cases, but it is still a prototype:

- it is a speed-oriented approximate surface, not a strict replacement for `irlba`
- the current tuned defaults trade a roughly `4.2-4.8%` singular-value envelope for the observed speedups
- a future C-level MLX bridge remains a legitimate optimization target, not a prerequisite for having a usable block Lanczos API

---

## References

- `R/irlba.R` — am_irlba wrapper with limitation documented
- `R/methods-dense.R` — numeric-left and matrix-left %*% methods (required for irlba dispatch)
- MLX docs: lazy evaluation, compile/fuse, linalg primitives
- Halko, Martinsson, Tropp (2011) — "Finding Structure with Randomness"
