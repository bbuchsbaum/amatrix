# Batch triangular solve

Solve B linear systems `A_b x_b = B_b` where each `A_b` is represented
by its Cholesky factor from
[`batch_chol`](https://bbuchsbaum.github.io/amatrix/reference/batch_chol.md).

## Usage

``` r
batch_solve(Ls, B)
```

## Arguments

- Ls:

  A list of `amChol` objects (output of `batch_chol`).

- B:

  A list of right-hand-side matrices/vectors, or a 3-D array
  `[n, k, B]`. Length / third dimension must match `Ls`.

## Value

A list of solution matrices (or vectors when each rhs is a vector).

## See also

[`batch_chol`](https://bbuchsbaum.github.io/amatrix/reference/batch_chol.md),
[`chol_solve`](https://bbuchsbaum.github.io/amatrix/reference/chol_solve.md)

## Examples

``` r
mats <- lapply(1:2, function(i) crossprod(matrix(rnorm(12), 4, 3)) + diag(3))
facs <- batch_chol(mats)
rhs <- lapply(1:2, function(i) rnorm(3))
sols <- batch_solve(facs, rhs)
all.equal(sols[[1]], solve(mats[[1]], rhs[[1]]))
#> [1] TRUE
```
