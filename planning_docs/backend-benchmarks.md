# Backend Benchmark Notes

This note records the current benchmark-driven routing policy for `amatrix` dense backends on this machine.

## Environment

- Host: Apple Silicon macOS
- Core package: `amatrix`
- Optional backends:
  - `amatrix.mlx`
  - `amatrix.arrayfire`

## Summary

- There are now two distinct performance regimes:
  - cold single-op execution on fresh host objects
  - chained steady-state execution where dense MLX residency can be reused
- Cold single-op execution is still mostly CPU-routed on this machine.
- Chained dense execution is where the current MLX backend starts to look like the actual end-state product.
- `arrayfire` remains a real native dense backend, but its current value is still primarily `matmul`, not chained execution.
- `arrayfire` now exposes resident hooks too, but the product still leans on MLX for the better-validated chained dense story on this machine.
- The current routed default path remains conservative:
  - `mlx` for dense `matmul` from `128x128` upward
  - `arrayfire` for dense `matmul` from `512x512` upward
  - CPU for `crossprod`, `tcrossprod`, `ewise`, `rowSums`, and `colSums` unless thresholds are overridden

## Current Default Thresholds

### `amatrix.mlx`

- `matmul`: `128`
- `crossprod`: `2048`
- `tcrossprod`: `2048`

### `amatrix.arrayfire`

- `matmul`: `512`
- `crossprod`: `2048`
- `tcrossprod`: `2048`
- `ewise`: `4096`
- `rowSums` / `colSums`: `4096`

These thresholds are intentionally conservative. They describe when the backend is allowed onto the default path, not whether the native implementation exists.

## Observed Cold Single-Op Benchmarks

These come from fresh-object runs via `tools/profile-backends.R`. They are meant to answer "what happens on one isolated operation?" rather than "what happens in a backend-resident chain?"

### Current Cold Routing

- `256x256`
  - `cpu matmul`: `0.0097 s`
  - `mlx matmul`: routed to `cpu`, `0.0083 s`
  - `arrayfire matmul`: routed to `cpu`, `0.0080 s`
- `512x512`
  - `cpu matmul`: `0.0517 s`
  - `mlx matmul`: routed to `cpu`, `0.0537 s`
  - `arrayfire matmul`: routed to `cpu`, `0.0550 s`
- `1024x1024`
  - `cpu matmul`: `0.4013 s`
  - `mlx matmul`: routed to `cpu`, `0.3790 s`
  - `arrayfire matmul`: routed to `cpu`, `0.4037 s`

In other words: the cold path is still effectively CPU-first with the current installed policy state and bridge costs.

### Cold `crossprod` / `tcrossprod`

- `1024x1024 cpu crossprod`: `0.6660 s`
- `1024x1024 mlx crossprod`: routed to `cpu`, `0.6510 s`
- `1024x1024 arrayfire crossprod`: routed to `cpu`, `0.6570 s`
- `1024x1024 cpu tcrossprod`: `0.2087 s`
- `1024x1024 mlx tcrossprod`: routed to `cpu`, `0.2067 s`
- `1024x1024 arrayfire tcrossprod`: routed to `cpu`, `0.2097 s`

## Observed Reduction / Elementwise Benchmarks

These are also cold-path numbers from `tools/profile-backends.R`.

- `512x512 ewise multiply`
  - CPU: about `0.0023 s`
  - ArrayFire path: still routed to `cpu`
- `1024x1024 rowSums`
  - CPU: about `0.0023 s`
  - ArrayFire path: still routed to `cpu`
- `1024x1024 colSums`
  - CPU: about `0.0023 s`
  - ArrayFire path: still routed to `cpu`
- `2048x2048 rowSums`
  - CPU: about `0.0063 s`
  - ArrayFire path: still routed to `cpu`

The current reduction and `ewise` bridge is therefore correctness-first, not performance-first.

## Observed Chained Dense Benchmarks

These come from `tools/benchmark-chained-dense.R`. Unlike the cold-path profile, they intentionally reuse persistent `amatrix` objects so backend residency can matter.

- `matmul_chain = ((x %*% y) * 2) + z`
- `cross_chain = (crossprod(x) * 2) + diag(ncol(x))`

Representative timings:

- `256x256`
  - `matmul_chain cpu`: `0.0160 s`
  - `matmul_chain mlx`: `0.0920 s`
  - `cross_chain cpu`: `0.0157 s`
  - `cross_chain mlx`: `0.0077 s`
- `512x512`
  - `matmul_chain cpu`: `0.0627 s`
  - `matmul_chain mlx`: `0.0113 s`
  - `cross_chain cpu`: `0.0887 s`
  - `cross_chain mlx`: `0.0133 s`
- `1024x1024`
  - `matmul_chain cpu`: `0.4043 s`
  - `matmul_chain mlx`: `0.0457 s`
  - `cross_chain cpu`: `0.6977 s`
  - `cross_chain mlx`: `0.0333 s`

This is the main current result for the project direction: cold isolated ops are not yet the compelling story, but chained dense execution with MLX residency is already materially better at medium and large sizes.

## Observed Shared-X Model-Core Benchmarks

These come from `tools/benchmark-model-core.R`. They are meant to answer a different question:

- what happens when one design matrix `X` is reused across many fits?

Current script setup:

- `n = 2000`
- `p = 32`
- `k = 8`
- `8` response matrices

Observed timings on this machine:

- `shared_x_lm`
  - `normal_cache_off`: about `0.0743 s`
  - `normal_cache_on`: about `0.0393 s`
  - `qr_cache_off`: about `0.0617 s`
  - `qr_cache_on`: about `0.0270 s`
  - `weighted_qr_cache_on`: about `0.0433 s`
- `shared_x_ridge`
  - `cache_off`: about `0.0777 s`
  - `cache_on`: about `0.0483 s`

This is the first direct evidence that the internal shared-`X` cache is doing real work for the product wedge the PRD cares about. It also shows that the QR-backed repeated-response path is already competitive on this workload while giving a more numerically robust least-squares route than normal equations.

## Observed SPD Cholesky Workflow Benchmarks

These come from `tools/benchmark-cholesky-runtime.R`.

On this machine, the reliable invocation is:

```bash
Rscript -e 'source("tools/benchmark-cholesky-runtime.R", local = TRUE)'
```

Direct `Rscript tools/benchmark-cholesky-runtime.R` still trips an MLX startup crash on this Apple Silicon setup, so the `-e` source path is the one to use for now.

Observed timings on this machine:

- `ridge_spd`, `768x768`, `rhs_cols = 64`
  - `factor`
    - `cpu`: about `0.0653 s`
    - `mlx`: about `0.0187 s`
  - `batched_solve`
    - `cpu`: about `0.0154 s`
    - `mlx`: about `0.0084 s`
  - `factor_plus_batched_solve`
    - `cpu`: about `0.0800 s`
    - `mlx`: about `0.0287 s`
- `kernel_spd`, `640x640`, `rhs_cols = 32`
  - `factor`
    - `cpu`: about `0.0360 s`
    - `mlx`: about `0.0108 s`
  - `batched_solve`
    - `cpu`: about `0.0062 s`
    - `mlx`: about `0.0042 s`
  - `factor_plus_batched_solve`
    - `cpu`: about `0.0424 s`
    - `mlx`: about `0.0166 s`

Quality note from the same benchmark run:

- factor reconstruction residuals were about `2e-7`
- solve errors versus the CPU reference were about `5e-7`

Interpretation:

- On representative ridge-like and kernel-like SPD workloads, MLX is already the faster end-to-end path on this machine.
- The biggest win is the factor step; the many-RHS triangular solve also helps, but less dramatically.
- ArrayFire remains explicit CPU-backed fallback/stub territory for this workflow.

## Observed Similarity Workload Benchmarks

These also come from `tools/benchmark-model-core.R`. They are meant to answer:

- what does the similarity / covariance workload cost on the current substrate?

Current script setup:

- `n = 4000`
- `p = 64`

Observed timings on this machine:

- `similarity`
  - `covariance`: about `0.0143 s`
  - `weighted_covariance`: about `0.0170 s`
  - `correlation`: about `0.0150 s`

These are CPU-path timings today, but they matter for product direction because they show that covariance/correlation can sit naturally on the same substrate and kernel layer as the regression workflows. The heavy second-moment step is already routed through `am_crossprod()`, so this workload is positioned to benefit from backend improvements without changing the public API.

## Observed QR Runtime Benchmarks

These come from `tools/benchmark-qr-runtime.R`. They compare four QR-related regimes:

- `base_r`
  - cached base QR helper calls such as `base::qr.coef(base_fac, y)`
- `mlx_resident_qr`
  - MLX QR factorization that keeps `Q` resident and returns only a resident key plus host `R`
- `mlx_native_resident`
  - MLX explicit QR payloads using backend-native helper ops against resident `Q`
- `mlx_compact`
  - the compact MLX QR path
  - on tall-skinny matrices this now means a TSQR-style blocked compact factor prototype
  - on other shapes it now keeps a lazy host compact factor and only materializes it if a compact helper actually needs it

The benchmark now also varies the number of right-hand sides (`rhs_cols`) to answer the practical question:

- when does the resident QR path start to matter for the shared-`X`, many-`Y` workload?

Current observed timings on this machine:

- `512x64`
  - `qr.factor_cold`
    - `base_r`: about `0.00140 s`
    - `mlx_resident_qr`: about `0.00265 s`
  - `qr.Q_materialize`
    - `base_r`: about `0.00235 s`
    - `mlx_native_resident`: about `0.00070 s`
  - `qr.coef`
    - `base_r`: about `0.00035 s`
    - `mlx_native_resident`: about `0.00145 s`
    - `mlx_compact`: about `0.00240 s`
  - `qr.qty`
    - `base_r`: about `0.00040 s`
    - `mlx_native_resident`: about `0.00120 s`
    - `mlx_compact`: about `0.00240 s`
  - `am_lm_fit(method = "qr")`, hot shared-`X` path
    - `base_qr_cached`: about `0.00040 s`
    - `mlx_native_resident`: about `0.00150 s`
    - `mlx_compact`: about `0.00120 s`

- `1024x128`
  - `qr.factor_cold`
    - `base_r`: about `0.01100 s`
    - `mlx_resident_qr`: about `0.01540 s`
  - `qr.Q_materialize`
    - `base_r`: about `0.01960 s`
    - `mlx_native_resident`: about `0.00080 s`
  - `qr.coef`
    - `base_r`: about `0.00140 s`
    - `mlx_native_resident`: about `0.00140 s`
    - `mlx_compact`: about `0.00360 s`
  - `qr.qty`
    - `base_r`: about `0.00160 s`
    - `mlx_native_resident`: about `0.00140 s`
    - `mlx_compact`: about `0.00380 s`
  - `am_lm_fit(method = "qr")`, hot shared-`X` path
    - `base_qr_cached`: about `0.00160 s`
    - `mlx_native_resident`: about `0.00160 s`
    - `mlx_compact`: about `0.00240 s`

Interpretation:

- The new compact MLX QR path is now a real algorithmic prototype for tall-skinny matrices, not just a relabeled bridge fallback.
- It is correct and factor-first.
- Moving the leaf-block work onto resident MLX `Q` factors materially improved the compact path at larger RHS counts.
- Even after that change, the resident-`Q` MLX path remains the speed winner on the flagship many-RHS workflow.
- That means the Stage 2 prototype is now structurally credible, but the next speed step still has to push more of the compact helper path below R and closer to backend-native execution.

### Many-RHS Scaling

The more interesting result is what happens as the number of response columns grows for a fixed `1024x128` design:

- `rhs_cols = 8`
  - `am_many_lm(method = "qr")`
    - `base_qr_cached`: about `0.0014 s`
    - `mlx_native_resident`: about `0.0016 s`
    - `mlx_compact`: about `0.0024 s`
- `rhs_cols = 32`
  - `am_many_lm(method = "qr")`
    - `base_qr_cached`: about `0.0056 s`
    - `mlx_native_resident`: about `0.0018 s`
    - `mlx_compact`: about `0.0064 s`
- `rhs_cols = 128`
  - `am_many_lm(method = "qr")`
    - `base_qr_cached`: about `0.0215 s`
    - `mlx_native_resident`: about `0.0020 s`
    - `mlx_compact`: about `0.0225 s`

That is the first genuinely interesting QR result in the repo:

- at small RHS counts, both MLX QR paths still pay fixed overhead relative to cached base QR
- once the number of right-hand sides grows, the resident helper path dominates the cached base QR workflow
- the new compact TSQR prototype improved materially after moving leaf-block helper work onto resident MLX factors, but it is still slower than the resident-`Q` MLX path
- the direct `am_many_lm()` QR path still shows the architecture paying off at the workflow level, and the compact path is now close enough that further native work is justified

This is much closer to the actual product wedge than one-off `qr.coef()` calls. It says the current architecture is beginning to pay off exactly where it should: shared design, many responses.

### Reduced Tall-Skinny TSQR Harness

`tools/benchmark-qr-tsqr.R` is the faster iteration harness for the compact QR project. It isolates:

- tall-skinny designs only
- `many_lm(..., method = "qr")`
- `qr.coef()`
- forced compact TSQR vs the current resident-`Q` path

Current observed timings on this machine:

- `1024x128`, `rhs_cols = 8`
  - `qr.coef`
    - `base_qr_cached`: about `0.0015 s`
    - `mlx_native_resident`: about `0.001625 s`
    - `mlx_compact_tsqr`: about `0.003125 s`
  - `am_many_lm(method = "qr")`
    - `base_qr_cached`: about `0.0015 s`
    - `mlx_native_resident`: about `0.001625 s`
    - `mlx_compact_tsqr`: about `0.003125 s`

- `1024x128`, `rhs_cols = 32`
  - `qr.coef`
    - `base_qr_cached`: about `0.0052 s`
    - `mlx_native_resident`: about `0.0018 s`
    - `mlx_compact_tsqr`: about `0.0028 s`
  - `am_many_lm(method = "qr")`
    - `base_qr_cached`: about `0.0052 s`
    - `mlx_native_resident`: about `0.0016 s`
    - `mlx_compact_tsqr`: about `0.0034 s`

- `1024x128`, `rhs_cols = 128`
  - `qr.coef`
    - `base_qr_cached`: about `0.020667 s`
    - `mlx_native_resident`: about `0.001667 s`
    - `mlx_compact_tsqr`: about `0.003667 s`
  - `am_many_lm(method = "qr")`
    - `base_qr_cached`: about `0.020667 s`
    - `mlx_native_resident`: about `0.0020 s`
    - `mlx_compact_tsqr`: about `0.0040 s`

- `4096x128`, `rhs_cols = 32`
  - `qr.coef`
    - `base_qr_cached`: about `0.0224 s`
    - `mlx_native_resident`: about `0.0020 s`
    - `mlx_compact_tsqr`: about `0.0048 s`
  - `am_many_lm(method = "qr")`
    - `base_qr_cached`: about `0.0222 s`
    - `mlx_native_resident`: about `0.0020 s`
    - `mlx_compact_tsqr`: about `0.0046 s`

- `4096x128`, `rhs_cols = 128`
  - `qr.coef`
    - `base_qr_cached`: about `0.085333 s`
    - `mlx_native_resident`: about `0.002667 s`
    - `mlx_compact_tsqr`: about `0.006333 s`
  - `am_many_lm(method = "qr")`
    - `base_qr_cached`: about `0.085667 s`
    - `mlx_native_resident`: about `0.002667 s`
    - `mlx_compact_tsqr`: about `0.006333 s`

Interpretation:

- The compact TSQR path remains materially better on raw `qr.coef()` than it was before the native top-stage coefficient reduction.
- After moving TSQR factor construction into a single native build bridge, making `r_stack` and top `R` resident-lazy, and pushing the full-rank TSQR coefficient solve into one native bridge, the compact path improved again while keeping the public factor surface unchanged.
- The end-to-end compact workflow still tracks the compact `qr.coef()` timings instead of collapsing back to cached base QR.
- The resident-`Q` MLX path is still the speed winner on the flagship workload.
- The remaining compact-path bottleneck is now mostly the cold factor build, not the hot many-RHS solve.
- The core test harness now sees the installed optional MLX backend instead of skipping it, so the TSQR path is exercised in the main suite rather than only in backend-local tests.

`tools/profile-many-lm-qr.R` shows the same picture directly for `4096x128`, `rhs=128`, `block_rows=512`:

- `cache`: about `0.137 s`
- `solve`: about `0.027 s`
- `assemble`: about `0.001 s`

Source-loaded spot checks after the lazy-factor cleanup in `R/qr.R` and `backends/amatrix.mlx/R/backend.R` show the current state more clearly:

- `768x256` non-TSQR `qr.factor_cold`
  - `mlx_resident_qr`: about `0.004 s`
  - `mlx_compact`: about `0.004 s`
  - both now leave the compact factor unmaterialized at `qr(x)` time
- `4096x128`, `rhs=128`, `block_rows=512` via `tools/profile-many-lm-qr.R`
  - native resident path: `cache ≈ 0.066 s`, `solve ≈ 0.005 s`, `assemble ≈ 0.003 s`
  - compact TSQR path: `cache ≈ 0.122 s`, `solve ≈ 0.017 s`, `assemble ≈ 0.014 s`

Current conclusion:

- the eager host compact-factor work is no longer the hidden cold-path cost on the MLX QR boundary
- the remaining compact cold cost on the flagship tall-skinny workflow is the TSQR build itself
- the resident-`Q` MLX path is still the premium QR implementation today
- the compact TSQR path is still structurally credible, but it is no longer the next blocking flagship issue

## Inspecting Execution Mode

Use the core introspection helpers to see which regime an object is in:

```r
x <- adgeMatrix(matrix(rnorm(16), 4, 4), preferred_backend = "mlx")
y <- (x %*% diag(4)) * 2

amatrix_residency_info(y)
amatrix_backend_plan(y, "matmul", y = diag(4))
amatrix_execution_info(y, ops = c("matmul", "crossprod"), y_map = list(matmul = diag(4)))
```

Important fields:

- `amatrix_backend_plan(... )$chosen_path`
  - `"cold"` or `"resident"`
- `amatrix_backend_matrix(... )$resident_reuse`
  - `TRUE` when the backend is being selected because the object is already resident
- `amatrix_backend_status()$residency_capable`
  - whether the backend implements the residency hooks at all

## Reproducing

Run:

```sh
Rscript /Users/bbuchsbaum/code/amatrix/tools/profile-backends.R
Rscript /Users/bbuchsbaum/code/amatrix/tools/benchmark-chained-dense.R
Rscript /Users/bbuchsbaum/code/amatrix/tools/benchmark-model-core.R
Rscript /Users/bbuchsbaum/code/amatrix/tools/benchmark-qr-runtime.R
Rscript /Users/bbuchsbaum/code/amatrix/tools/benchmark-qr-tsqr.R
Rscript /Users/bbuchsbaum/code/amatrix/tools/benchmark-svd-factor.R
```

The first script reports cold single-op routing and timings. The second reports chained steady-state behavior on persistent dense objects.
The third reports repeated shared-`X` model-core workloads with cache reuse enabled and disabled.
`benchmark-svd-factor.R` currently reports CPU/reference factor timings only; for Apple Silicon MLX steady-state SVD timing, use the direct `Rscript -e` commands documented in [gpu-svd-analysis.md](/Users/bbuchsbaum/code/amatrix/docs/gpu-svd-analysis.md). `tools/print-svd-factor-calibration.R` prints a small calibration grid of those commands.

## ArrayFire rsvd Product Status: Experimental

**Decision (2026-04-06): experimental — correct algorithm, unvalidated in CI.**

### Evidence

- The implementation (`backends/amatrix.arrayfire/R/backend.R:234`) is a standard randomized SVD with power iteration. The algorithm is algebraically correct.
- The large matmuls (sketch, power iterations, projection) use `amatrix_arrayfire_matmul_correct_bridge` and `amatrix_arrayfire_crossprod_correct_bridge` — the "correct" bridge variants that were fixed in commit c21561d (June 2026).
- The conformance test (`tests/testthat/test-cross-backend-conformance.R:296`) is guarded by `skip_if_not(amatrix_arrayfire_is_available())` and therefore always skips on this machine and in CI. The quality checks have not run against real AF hardware.
- No performance benchmark comparing AF rsvd vs MLX rsvd has been collected.
- The active threshold is `min(nrow, ncol) >= 400` (configurable via `amatrix.arrayfire.rsvd_min_dim`).

### Status

- **Not fallback-only**: the path is real — it will activate when AF hardware is present and the threshold is met.
- **Not product-ready**: no real-hardware quality or performance validation.
- **Experimental**: users with AF hardware can enable and exercise it; the recommended path for Apple Silicon is the MLX backend.

No roadmap document should claim ArrayFire rsvd as either "not implemented" or "fully validated at parity with MLX". It is a functioning but unvalidated feature gate.
