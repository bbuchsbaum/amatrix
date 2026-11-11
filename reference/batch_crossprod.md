# Batch crossproduct

Compute `t(A_b) %*% A_b` for each matrix in a batch.

## Usage

``` r
batch_crossprod(A)
```

## Arguments

- A:

  A list of numeric matrices, or a 3-D array `[n, p, B]`.

## Value

A list of `p x p` crossproduct matrices.

## See also

[`batch_chol`](https://bbuchsbaum.github.io/amatrix/reference/batch_chol.md)

## Examples

``` r
mats <- lapply(1:3, function(i) matrix(rnorm(12), 4, 3))
xtx <- batch_crossprod(mats)
dim(xtx[[1]])
#> [1] 3 3
```
