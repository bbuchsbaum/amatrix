# amatrix Benchmark Results

**Machine:** Apple Silicon (M-series), macOS 23.3  
**R version:** 4.4+ with Accelerate BLAS  
**amatrix.mlx:** available (float32 fast path)  
**Date:** 2026-04-06  

All timings are median wall-clock seconds unless noted otherwise.
Reproduce with the scripts in `tools/` — see `docs/quality-tracking.md §3` for protocol.

---

## 1. many_lm — batched OLS with QR cache

**Script:** `tools/benchmark-flagship-many-y.R`  
**Reference:** `lm.fit()` loop over k response columns  
**Model:** `am_many_lm(X, Y, method="qr", cache=TRUE)`

The cache stores the QR decomposition on first call and reuses it across
all k responses. The speedup vs the lm_loop reference grows with k.

| Scenario | n | p | k | lm_loop | amatrix CPU | amatrix MLX | CPU speedup | MLX speedup |
|----------|---|---|---|---------|-------------|-------------|-------------|-------------|
| small | 1 000 | 50 | 100 | 0.185 s | 0.007 s | 0.004 s | 26× | 53× |
| large | 5 000 | 100 | 500 | 16.337 s | 0.297 s | 0.017 s | 55× | **993×** |

**Interpretation:** The dominant cost in `lm_loop` is re-factorizing X for each
of the k responses. `am_many_lm` computes the QR once and solves k RHS in a
batched triangular solve. MLX benefits additionally from float32 GEMM throughput
at large n×k. The 993× figure is real and repeatable — it reflects both the
cache and the GPU compute advantage at n=5000.

---

## 2. Cholesky factor + batched solve

**Script:** `tools/benchmark-cholesky-runtime.R`  
**Workloads:**
- `ridge_spd`: A = X'X + λI, 768×768, 64 RHS columns  
- `kernel_spd`: A = K_rbf + δI, 640×640, 32 RHS columns  

MLX factorizes in float32; solve is GPU triangular solve. CPU uses Accelerate
LAPACK (float64). Relative errors are vs exact CPU solve.

| Workload | Size | Phase | CPU | MLX | MLX speedup | MLX rel_err |
|----------|------|-------|-----|-----|-------------|-------------|
| ridge_spd | 768×768 | factor | 71 ms | 20 ms | **3.6×** | 1.9e-7 |
| ridge_spd | 768×768 | batched_solve (64 RHS) | 16 ms | 8 ms | **2.0×** | 5.7e-7 |
| ridge_spd | 768×768 | factor + solve | 88 ms | 30 ms | **3.0×** | 5.7e-7 |
| kernel_spd | 640×640 | factor | 41 ms | 13 ms | **3.1×** | 2.0e-7 |
| kernel_spd | 640×640 | batched_solve (32 RHS) | 6 ms | 4 ms | **1.5×** | 5.0e-7 |
| kernel_spd | 640×640 | factor + solve | 46 ms | 19 ms | **2.5×** | 5.0e-7 |

**Interpretation:** Cholesky is the first operation where MLX shows consistent
speedup at these sizes. The float32 precision loss is ~1e-7 (relative to CPU
float64) — acceptable for all use cases that don't require double-precision
factorization. Use `precision = "strict"` to force CPU float64 when needed.

---

## 3. SVD factor — truncated rank-k decomposition

**Script:** `tools/benchmark-svd-factor.R`  
**Method:** rank-20 truncated SVD, n_oversamples=10, n_iter=2  
**Reference:** `base::svd` (full), `irlba::svdr` (randomized)

| Case | base::svd | irlba::svdr | am_svd_factor(exact) | svdr speedup vs svd | svdr rel_sv_err |
|------|-----------|-------------|---------------------|---------------------|-----------------|
| 300×240 | 0.036 s | 0.011 s | 0.088 s | 3.3× | 0.059 |
| 500×400 | 0.163 s | 0.028 s | 0.165 s | 5.8× | 0.079 |
| 1000×800 | 1.292 s | 0.105 s | 1.306 s | 12.3× | 0.087 |
| 2000×1600 | 10.335 s | 0.407 s | 10.316 s | 25.4× | 0.096 |

**Note:** MLX rows are omitted here — nested Rscript MLX benchmarking is
unstable on this machine (see `docs/gpu-svd-analysis.md` for MLX spot numbers).
`am_svd_factor(exact)` wraps `base::svd` and adds factor-object overhead; the
`rsvd`/`auto` methods via `irlba::svdr` deliver 3–25× faster truncated SVD at
the cost of ~10% singular-value relative error (acceptable for most downstream
uses like PCA projection and low-rank approximation).

---

## 4. Dense matrix chains — matmul, crossprod, ewise

**Script:** `tools/benchmark-chained-dense.R`  
**Workloads:**  
- `matmul_chain`: `(X %*% Y) * 2 + Z`  
- `cross_chain`: `crossprod(X) * 2 + diag(p)`

| Workload | n | CPU | MLX | MLX speedup |
|----------|---|-----|-----|-------------|
| matmul_chain | 256 | 10 ms | 36 ms | 0.28× |
| matmul_chain | 512 | 42 ms | 41 ms | 1.0× |
| matmul_chain | 1024 | 307 ms | 307 ms | 1.0× |
| cross_chain | 256 | 12 ms | 11 ms | 1.1× |
| cross_chain | 512 | 69 ms | 69 ms | 1.0× |
| cross_chain | 1024 | 551 ms | 567 ms | 0.97× |

**Interpretation:** Apple Silicon's Accelerate BLAS is extremely competitive
with MLX for square matrix operations up to n=1024. No MLX speedup is observed
here. This is expected — see `docs/backend-benchmarks.md` for the routing
thresholds: GPU dispatch is only beneficial for large matrices or for ops where
float32 throughput matters (many_lm, Cholesky). At n=256 the Metal kernel launch
overhead (~30ms) dominates, making MLX materially slower.

---

## 5. Regression baseline (cold + warm)

**Script:** `tools/benchmark-regression.R`

The regression harness covers 6 ops × 3 sizes × 2 variants (cold/warm) × 2
backends. Selected warm-path results at medium size (1024×128):

| Op | CPU warm | MLX warm | MLX speedup |
|----|----------|----------|-------------|
| matmul | 5.9 ms | 5.6 ms | 1.05× |
| crossprod | 9.9 ms | 9.6 ms | 1.03× |
| covariance | 12.2 ms | 12.3 ms | 0.99× |
| many_lm (cached) | 1.2 ms | 1.3 ms | 0.93× |
| rsvd (k=10) | 13.8 ms | 13.8 ms | 1.00× |

At large size (4096×128), warm-path speedups remain ~1.0× for all ops. This is
consistent with item 4 above — Accelerate BLAS matches MLX for non-square,
moderately sized dense products.

---

## Summary

| Op family | Headline speedup | Notes |
|-----------|-----------------|-------|
| `am_many_lm` (cached, k=500) | **993× vs lm_loop** | QR cache + MLX batched solve |
| `am_many_lm` (cached, k=100) | **53× vs lm_loop** | QR cache dominates |
| Cholesky (factor + solve) | **2.5–3.0× MLX vs CPU** | float32; set `precision="strict"` for float64 |
| Truncated SVD (irlba path) | **3–25× vs base::svd** | CPU-only; scales with matrix size |
| Dense BLAS (matmul, crossprod) | **~1.0× at n≤4096** | Accelerate BLAS = MLX at these sizes |

The clear performance story: **amatrix wins through caching and algorithm
selection**, not raw GPU GEMM throughput. `am_many_lm` with `cache=TRUE` is the
flagship case — the QR decomposition is amortized across all k responses, which
is the dominant factor. MLX provides additional leverage at large n where GPU
compute throughput exceeds what Accelerate delivers in float64.
