# Doubly-stochastic scaling via Sinkhorn-Knopp iterations

Alternates row and column normalization until the matrix is
approximately doubly stochastic. When the chosen backend supports
resident broadcast and reduction kernels, the hot loop stays on device
via `resident_handle` and returns a deferred `adgeMatrix` bound to the
resident result.

## Usage

``` r
sinkhorn(
  A,
  max_iter = 200L,
  tol = 1e-08,
  check_every = 5L,
  eps = 1e-15,
  mode = "fast",
  backend = NULL,
  return_info = FALSE
)
```

## Arguments

- A:

  A dense numeric matrix or `adgeMatrix`. Sparse inputs are not yet
  supported in this surface.

- max_iter:

  Maximum number of Sinkhorn iterations.

- tol:

  Convergence tolerance on the maximum row/column sum error.

- check_every:

  Check convergence every `check_every` iterations.

- eps:

  Floor applied to row/column sums before division.

- mode:

  Execution mode used when coercing a plain matrix. Default `"fast"`
  allows accelerated backends to use lower precision.

- backend:

  Backend name used when coercing a plain matrix. Ignored when `A` is
  already an `adgeMatrix`.

- return_info:

  When `TRUE`, return convergence metadata alongside the scaled matrix.

## Value

By default, an `adgeMatrix`. With `return_info = TRUE`, a list
containing `result`, `iterations`, `converged`, `row_error`,
`col_error`, `backend`, and `method`.

## See also

[`dist_matrix`](https://bbuchsbaum.github.io/amatrix/reference/dist_matrix.md)

## Examples

``` r
A <- abs(matrix(rnorm(16), nrow = 4)) + 0.1
S <- sinkhorn(A, max_iter = 50L)
# Row sums should be close to 1
rowSums(as.matrix(S))
#> [1] 1 1 1 1
```
