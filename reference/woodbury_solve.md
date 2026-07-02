# Solve a linear system using the Woodbury matrix identity

Computes \\(A + UCV)^{-1} b\\ in \\O(nk^2 + k^3)\\ time using the
Woodbury matrix identity, avoiding an \\O(n^3)\\ refactorisation of the
updated matrix.

## Usage

``` r
woodbury_solve(A_factor, U, b, V = NULL, C_inv = NULL)
```

## Arguments

- A_factor:

  An `amChol` object from
  [`chol_factor()`](https://bbuchsbaum.github.io/amatrix/reference/chol_factor.md),
  or a square numeric matrix that is automatically Cholesky-factored.

- U:

  Numeric matrix of shape `[n, k]`; low-rank left factor.

- b:

  Numeric matrix of shape `[n, rhs]`; right-hand side(s).

- V:

  Numeric matrix of shape `[k, n]`; low-rank right factor. Defaults to
  `t(U)` (symmetric update).

- C_inv:

  Numeric matrix of shape `[k, k]`; inverse of the central factor \\C\\.
  Defaults to `diag(k)` (pure rank-k update with \\C = I\\).

## Value

Numeric matrix of shape `[n, rhs]`: the solution \\(A + UCV)^{-1} b\\.

## See also

[`woodbury_logdet`](https://bbuchsbaum.github.io/amatrix/reference/woodbury_logdet.md),
[`chol_factor`](https://bbuchsbaum.github.io/amatrix/reference/chol_factor.md)

## Examples

``` r
A <- crossprod(matrix(rnorm(25), 5)) + diag(5)
U <- matrix(rnorm(10), 5, 2)
b <- rnorm(5)
x <- woodbury_solve(A, U, b)
length(x)
#> [1] 5
```
