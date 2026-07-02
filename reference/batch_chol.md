# Batch Cholesky factorization

Factorize B symmetric positive-definite matrices in parallel. Each
matrix is dispatched through the same backend as
[`chol_factor`](https://bbuchsbaum.github.io/amatrix/reference/chol_factor.md),
so MLX GPU acceleration applies to every element when available.

## Usage

``` r
batch_chol(A)
```

## Arguments

- A:

  A list of square numeric matrices, or a 3-D array `[n, n, B]`.

## Value

A list of `amChol` objects, one per input matrix.

## See also

[`chol_factor`](https://bbuchsbaum.github.io/amatrix/reference/chol_factor.md),
[`batch_solve`](https://bbuchsbaum.github.io/amatrix/reference/batch_solve.md)
