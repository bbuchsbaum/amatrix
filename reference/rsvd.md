# GPU-native randomized SVD (Halko et al. 2011)

Computes a truncated SVD via randomized projection entirely on the GPU.
All QR, matmul, and SVD steps stay on device; a single `mlx_eval`
materializes the results. Falls back to `irlba::svdr` on CPU if no GPU
backend with rsvd support is active.

## Usage

``` r
rsvd(x, k, n_oversamples = 10L, n_iter = 2L, ...)
```

## Arguments

- x:

  An `adgeMatrix` or plain numeric matrix.

- k:

  Number of singular values/vectors to compute.

- n_oversamples:

  Extra columns for the random projection (default 10). Increasing this
  improves accuracy at modest cost.

- n_iter:

  Number of power-iteration passes (default 2). More passes give better
  accuracy for matrices with slowly decaying spectra.

- ...:

  Ignored (for forward compatibility).

## Value

A list with components `u` (m x k), `d` (length-k singular values,
decreasing), and `v` (n x k).

## References

Halko, N., Martinsson, P. G., & Tropp, J. A. (2011). Finding structure
with randomness: Probabilistic algorithms for constructing approximate
matrix decompositions. *SIAM Review*, 53(2), 217-288.

## See also

[`block_lanczos`](https://bbuchsbaum.github.io/amatrix/reference/block_lanczos.md),
[`eigh`](https://bbuchsbaum.github.io/amatrix/reference/eigh.md)

## Examples

``` r
A <- matrix(rnorm(200), nrow = 20)
res <- rsvd(A, k = 3L)
length(res$d)
#> [1] 3
```
