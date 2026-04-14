# Fast Singular Value Decomposition on GPU
**S9226 — GPU Technology Conference**
Lung-Sheng Chien & Samuel Rodriguez Bernabeu, NVIDIA

---

## Outline

1. Issues of GESVD
2. Approximate SVD (GESVDA)
3. Randomized SVD (SVDR)
4. Conclusions

---

## Part 1: Issues of GESVD

### General SVD

Standard form: A = U S V^T

- **LAPACK GESVD**: most popular routine, based on QR iteration
- **cuSOLVER** provides two routines:
  - `GESVD`: same algorithm as LAPACK
  - `GESVDJ`: two-sided Jacobi method
- **Tall skinny SVD** is the common data-analytics use case:
  - Singular vectors required
  - Only a few large singular values needed
  - Typical size: 1e6 rows × 100 columns

### Strategy of GESVD

For tall skinny (m >> n):

1. **QR factorization** to preprocess: A (m×n) = Q (m×n) · R (n×n)
2. **SVD of square R**: R = U · S · V^T  via GEBRD + BDSQR + ORGBR

**Problems on GPU:**
- GPU does not perform well on tall-skinny QR factorization
- GPU does not perform well on QR iteration for small matrix

### Performance on Square Matrix (V100 vs MKL/8-core i9-7900X)

Runtime in milliseconds. SGEMM = theoretical peak.

| n    | cuSOLVER SGESVD | cuSOLVER SGESVDJ | MKL SGESVD | SGEMM  |
|------|-----------------|------------------|------------|--------|
| 32   | 0.04            | 0.12             | 0.11       | 1      |
| 64   | 0.13            | 0.45             | 0.12       | 7      |
| 128  | 0.48            | 1.63             | 1.16       | 74     |
| 256  | 1.31            | 4.84             | 2.80       | 558    |
| 512  | 5.22            | 13.56            | 10.26      | 2,828  |
| 1024 | 18.97           | 27.73            | 8.33       | 8,586  |
| 2048 | 63.15           | 40.90            | 19.80      | 12,514 |
| 4096 | 152.07          | 58.08            | 8.03       | 13,366 |
| 8192 | 264.11          | 49.19            | 5.52       | 13,956 |

- SVD flop count is 2N^3, same as SGEMM — but SVD runtime is ~50× slower
- Jacobi (GESVDJ) faster than QR iteration (GESVD) for n < 1024

### Performance on Tall Skinny Matrix (N=32, M varies)

| M         | SGEQRF (sec) | SGESVDJ (sec) | QR ratio |
|-----------|-------------|---------------|----------|
| 1,000     | 0.00021     | 0.00128       | 0.17     |
| 10,000    | 0.00058     | 0.00147       | 0.40     |
| 100,000   | 0.00524     | 0.00654       | 0.80     |
| 1,000,000 | 0.05897     | 0.06336       | 0.93     |

QR factorization complexity is 2MN^2 flops (proportional to M).
As M grows, QR becomes the bottleneck (93% of total time at M=1e6).

### Weakness of QR Factorization (M=8192, N varies)

| N    | SGEQRF (Gflops) |
|------|-----------------|
| 32   | 32.6            |
| 64   | 66.3            |
| 128  | 116.6           |
| 256  | 159.5           |
| 512  | 338.1           |
| 1024 | 627.6           |
| 2048 | 990.8           |
| 4096 | 2487.6          |

Despite 2MN^2 complexity, runtime is proportional to N because:
- Only the **trailing matrix** uses BLAS-3 (which is negligible for tall skinny A)
- Runtime dominated by **panel factorization** (mainly BLAS-1)

---

## Part 2: Approximate SVD (GESVDA)

### Key Insight

Instead of QR-factorization path: A = U S V^T

Use the **normal equations** path:  A^T A = V S^2 V^T

Strategy: **GEMM + EIG** (replaces QR factorization with matrix multiply)

### Technical Issues

**Rounding error:**
- A^T A has rounding error proportional to ||A||^2
- Use **double-precision GEMM (DGEMM)** to control rounding error

**Performance:**
- A^T A is N×N — small matrix compared to tall skinny A
- Regular GEMM doesn't perform well at this size
- Need special GEMM to improve performance
- Use **Jacobi method (DSYEVJ)** for eigenpairs of A^T A (faster than QR on small matrix)

### GESVDA Algorithm

1. **B = A^T A** via DGEMM  (reduces rounding errors vs SGEMM)
2. **(S, V) = eig(B)** via DSYEVJ (Jacobi); adjust stopping criteria for performance
3. **U = A V S^{-1}** via DGEMM + scaling
   - Note: left singular vectors inaccurate when singular value is small

### Quality of Solution

- Right singular vectors: always accurate to 1e-6
- Singular values and left singular vectors depend on M and N

For N ≤ 100 and M ≤ 176,000 × N:
- If S_l(A) ≥ 2.65e-3 · ||A||_F, then singular values and vectors accurate to 1e-6
- Example: largest singular value S_l(A) ≥ (1/sqrt(N)) · ||A||_F is always accurate to 1e-6

### Performance of GESVDA (N=32 fixed)

| M         | SGEQRF (sec) | SGESVDJ (sec) | SGESVDA (sec) | QR ratio | Speedup |
|-----------|-------------|---------------|---------------|----------|---------|
| 1,000     | 0.00021     | 0.00128       | 0.00077       | 0.17     | 1.67x   |
| 10,000    | 0.00058     | 0.00147       | 0.00078       | 0.40     | 1.89x   |
| 100,000   | 0.00524     | 0.00654       | 0.00118       | 0.80     | 5.55x   |
| 1,000,000 | 0.05897     | 0.06336       | 0.00376       | 0.93     | **16.84x** |

Speedup comes entirely from replacing QR factorization with GEMM.

### GESVDA Breakdown (N=32, as % of total time)

| M         | DGEMM | DSYEVJ | Compute U | Residual |
|-----------|-------|--------|-----------|----------|
| 1,000     | 0.13  | 0.78   | 0.03      | 0.10     |
| 10,000    | 0.15  | 0.74   | 0.03      | 0.10     |
| 100,000   | 0.33  | 0.46   | 0.13      | 0.10     |
| 1,000,000 | 0.40  | 0.16   | 0.35      | 0.11     |

- DGEMM is only ~40% of total time
- DSYEVJ cost is O(N^3), independent of M — fraction decreases as M grows
- "Compute U" is slower than "Residual" because it requires double precision

### Performance of Batched GESVDA (N=35 fixed)

| M         | batchSize | SGESVDJ (sec) | SGESVDA (sec) | Speedup |
|-----------|-----------|---------------|---------------|---------|
| 16,384    | 1         | 0.0022        | 0.0012        | 1.77x   |
| 16,384    | 32        | 0.0562        | 0.0076        | 7.41x   |
| 65,536    | 1         | 0.0051        | 0.0015        | 3.41x   |
| 65,536    | 32        | 0.0977        | 0.0141        | 6.94x   |
| 1,048,576 | 1         | 0.0770        | 0.0062        | 12.45x  |
| 1,048,576 | 16        | 0.5329        | 0.0775        | 6.87x   |

- SGESVDJ uses OpenMP (CPU threads); GESVDA uses multi-stream GPU
- OpenMP can parallelize batchSize GEQRF in parallel (GEQRF is 40%+ of GESVDJ), so batch32 speedup is limited

### Batched GESVDA Breakdown (N=35, as % of total)

| M         | batchSize | DGEMM | DSYEVJ | Compute U | Residual |
|-----------|-----------|-------|--------|-----------|----------|
| 16,384    | 32        | 0.21  | 0.50   | 0.13      | 0.15     |
| 65,536    | 32        | 0.32  | 0.22   | 0.29      | 0.17     |
| 1,048,576 | 16        | 0.36  | 0.03   | 0.43      | 0.17     |

- In batched mode, DSYEVJ drops to 3% (custom batched Jacobi avoids kernel-launch overhead)
- "Compute U" becomes the bottleneck in large batches

### Double-Double (fp128) GEMM for Higher Accuracy

Goal: low-rank SVD accurate to 1e-14 (vs 1e-6 for standard GESVDA)

- Uses QD package: `LGEMM: C(dd) += A(d) * B(dd)` 
- Only useful when M > 100,000

| M         | DGEQRF (sec) | DGEMM / QR ratio | LGEMM / QR ratio |
|-----------|-------------|------------------|------------------|
| 1,000     | 0.00023     | 1.49             | 0.38             |
| 10,000    | 0.00077     | 7.39             | 0.96             |
| 100,000   | 0.00822     | 26.19            | 1.88             |
| 1,000,000 | 0.08447     | 14.53            | 5.34             |

### GESVDA Conclusions

- Replaces SGEQRF with DGEMM → up to **16x speedup**
- Inhouse batched eigenvalue solver avoids kernel-launch bottleneck
- Good quality for singular values and vectors in common use cases
- **Shipped in CUDA 10.1** with batched API

---

## Part 3: Randomized SVD (SVDR)

### Motivation: Low-Rank Approximation

Given matrix A, find best rank-k approximation:
- **Fixed precision**: find smallest k such that ||A - A_k|| < ε
- **Fixed rank**: given k, minimize ||A - A_k||
- Truncated SVD gives optimal rank-k approximation [Eckart-Young-Mirsky]
- But full SVD costs O(n^3) — too expensive

### When Approximate is Enough

Randomized SVD computes top-k eigenpairs to *sufficient* accuracy:
- Data analytics / PCA / clustering: 1e-2 accuracy may be enough
- Physics simulations: 1e-2 may be useless

**Highlights:**
- Reduced time and space complexity
- Preserves sparsity of A
- One-pass or streaming algorithms (read A only once)

### Core Idea: Range Finder (Halko et al., SIAM Review 2011)

**Algorithm** (inputs: A, k; output: Q, B such that A ≈ QB):
```
1. Ω = SketchingMatrix(A, O(k))      -- random projection matrix
2. C = A Ω                            -- project A into low-dimensional space
3. Q = orth(C)                        -- orthogonalize
4. B = Q^T A                          -- compress A
```

**Complexity:**
- Step 1: O(m × k)
- Step 2: O(m × n × k)
- Step 3: O(m × k^2)
- Step 4: O(m × n × k)

**Error bound** (with q power iteration steps):
||A - QQ^H A|| ≤ λ_{k+1} + (1 + 4·sqrt(2m/(k-1)))^(1/(2q+1)) · λ_{k+1}

### Sketching Matrix Options

An m-by-O(k) matrix that captures the column space of A:
- **Gaussian projection** — easy to apply, expensive to construct
- **Subsampled Randomized Hadamard Transform (SRHT)**
- **Count sketch**
- **Leverage-score subsampling** (sparse cases)

### SVDR Algorithm

**Inputs:** A, k  
**Outputs:** Û, Σ_k, V_k^H

```
1. Ω = SketchingMatrix(A, O(k))
2. C = A Ω
3. [Q, R] = qr(C)
4. [Ũ, Σ̃, Ṽ] = svd(Q^H A)
5. Û = Q Ũ
6. A_k = Û Σ_k V_k^H
```

Error bound:
||A - Ã_k|| = ||A - Q̃ Σ̃_k Ṽ^H|| ≤ λ_{k+1} + (1 + 4·sqrt(2m/(k-1)))^(1/(2q+1)) · λ_{k+1}

### Error Metric for Numerical Experiments

Don't measure:
- ||U - Ũ_k|| < ε  or  ||V - Ṽ_k|| < ε

Instead measure relative residual:
- ||A - Ã_k|| = (1 + η) λ_{k+1}
- η close to 0 means near-optimal approximation

### Spectral Norm Estimator (Magdon-Ismail & Malik, arXiv:1104.2076)

To adaptively determine when the approximation is good enough:

p_j(A) = sqrt( ||(A^H A)^j ω̃||_2 / ||(A^H A)^{j-1} ω̃||_2 )

P[ p_j(A) ≥ ||A||_2 / 10 ] > 1 - 4·sqrt(n/(j-1)) · 100^{-j}

**Bottom line: 6 iterations estimate ||A||_2 within a factor of 10.**

### Test Cases for Accuracy Experiments

Accuracy depends on spectral gap. Three test cases with 100 eigenvalues:
- **Fast decay**: eigenvalues drop steeply after k (large gap)
- **S-shape**: flat plateau then steep drop at k=20
- **Slow decay**: gradual decay (small gap, hard case)

### Accuracy Results

**Fast decay** (k=20 threshold):
- 0 power iterations: good accuracy (η ~ 1e-7 to 1e-8) near k=20
- 1 power iteration: noisy improvement at larger k

**S-shape** (transition at k=20):
- 0 iters: accuracy collapses above k=20 (η → 1 when gap disappears)
- 1 iter: significantly better, but still degrades above k=20

**Slow decay**:
- 0 iters and 1 iter nearly identical — power iteration doesn't help
- Error η ~ 1e-2 throughout (no spectral gap to exploit)

**Key insight**: SVDR accuracy is controlled by the spectral gap, not by k alone.

### Speedup: RSVD over GESVD (Rank-10)

Speedup matrix (rows = #rows of A, cols = #cols of A):

- At 2048×2048: ~103x speedup
- At 4096×4096: ~82x speedup
- Smallest matrices (256×256): only ~2.5x (GEMM not efficient at small sizes)

### Speedup: RSVD over GESVD (Rank-20)

- At 2048×2048: ~43x speedup
- At 4096×4096: ~34x speedup
- Speedup scales inversely with k (more columns in sketch = less speedup vs GESVD)

---

## Part 4: SVDR Conclusions

**Concerns:**
- No error estimator for spectrum accuracy; theoretical bounds are very pessimistic in practice

**Results:**
- SVDR provides good accuracy for top-k eigenvalues
- Works out of the box for k < 10 in practice
- Good alternative for PCA — but **not a substitute for GESVD**
- Can get impressive performance if you know your data's spectral structure

**Status:** Internal research project at time of talk — not included in CUDA 10.1. (Feedback welcomed.)

---

## Summary Table: Which Algorithm to Use?

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Square matrix, n < 1024 | GESVDJ | Jacobi faster than QR iteration |
| Square matrix, n > 2048 | GESVD | QR iteration catches up |
| Tall skinny, full accuracy | GESVDA | Up to 16x faster via GEMM + Jacobi |
| Tall skinny, approximate (data analytics, PCA) | SVDR | Up to 100x faster for rank-10 |
| Need 1e-14 accuracy | GESVDA + LGEMM (fp128) | Only worth it for M > 100,000 |
