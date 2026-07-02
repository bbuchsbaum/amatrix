# Backend-dispatched covariance matrix

Computes the (possibly weighted) sample or population covariance matrix
of `X`. Sparse inputs use a memory-efficient path; dense inputs use a
fused GPU kernel when available, otherwise fall back to CPU centering
followed by a GPU cross-product.

## Usage

``` r
covariance(X, center = TRUE, sample = TRUE, weights = NULL, block_size = NULL)
```

## Arguments

- X:

  Numeric matrix, `adgeMatrix`, or sparse `adgCMatrix` / `sparseMatrix`,
  shape `[n, p]`.

- center:

  Logical. When `TRUE` (default), columns are mean-centred before
  computing the cross-product.

- sample:

  Logical. When `TRUE` (default), divides by `n - 1` (sample
  covariance); when `FALSE`, divides by `n`.

- weights:

  Optional numeric vector of length `n` with non-negative observation
  weights. When supplied, a weighted covariance is computed.

- block_size:

  Optional integer. When set, the cross-product is computed in
  column-blocks of this size to limit peak memory. Ignored for sparse
  inputs.

## Value

An `adgeMatrix` of shape `[p, p]`.

## See also

[`many_lm`](https://bbuchsbaum.github.io/amatrix/reference/many_lm.md)

## Examples

``` r
X <- matrix(rnorm(200), nrow = 40)
C <- covariance(X)
dim(C)
#> [1] 5 5
```
