# Stochastic trace estimator (Hutchinson)

Estimates \\\mathrm{tr}(A)\\ or \\\mathrm{tr}(A^{-1})\\ using
Hutchinson's method with `k` Rademacher probe vectors. Supply `solve_fn`
to estimate the trace of an inverse without forming it explicitly.

## Usage

``` r
trace_estim(A = NULL, k = 30L, seed = NULL, solve_fn = NULL, n = NULL)
```

## Arguments

- A:

  Square matrix or `aMatrix`; required when `solve_fn` is `NULL`.

- k:

  Integer number of Rademacher probe vectors. Default `30L`.

- seed:

  Optional integer random seed for reproducibility.

- solve_fn:

  Optional function `function(V)` that returns `A^{-1} %*% V`; use this
  to estimate \\\mathrm{tr}(A^{-1})\\ without materialising the inverse.

- n:

  Integer dimension of the matrix; required when `solve_fn` is supplied.

## Value

A single numeric scalar estimate of the trace.

## Examples

``` r
A <- crossprod(matrix(rnorm(25), 5, 5)) + 5 * diag(5)
trace_estim(A, k = 50L, seed = 1L)
#> [1] 52.56089
```
