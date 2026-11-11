# Solve many right-hand-side batches with one Cholesky factor

This is the same operation as repeatedly calling
`chol_solve(factor, B[[i]])`, but it packs all RHS batches into one wide
solve and then splits the result. BLAS/GPU backends generally amortize
launch and dispatch overhead much better on one wide RHS than on many
small independent solves.

## Usage

``` r
chol_solve_batches(factor, B)
```

## Arguments

- factor:

  An `amChol` object from
  [`chol_factor`](https://bbuchsbaum.github.io/amatrix/reference/chol_factor.md).

- B:

  A list of RHS vectors/matrices, or a 3-D array `[n, k, batch]`.

## Value

A list of solution vectors/matrices.

## Examples

``` r
S <- crossprod(matrix(rnorm(20), 5, 4)) + diag(4)
fac <- chol_factor(adgeMatrix(S))
rhs <- list(rnorm(4), matrix(rnorm(8), 4, 2))
sols <- chol_solve_batches(fac, rhs)
all.equal(sols[[1]], solve(S, rhs[[1]]))
#> [1] TRUE
```
