# Log-determinant via the Woodbury matrix determinant lemma

Computes \\\log\|A + UCV\|\\ using the matrix determinant lemma in
\\O(nk^2 + k^3)\\ time, reusing an existing Cholesky factor of \\A\\.

## Usage

``` r
woodbury_logdet(A_factor, U, V = NULL, C_inv = NULL)
```

## Arguments

- A_factor:

  An `amChol` object from
  [`chol_factor()`](https://bbuchsbaum.github.io/amatrix/reference/chol_factor.md),
  or a square numeric matrix that is automatically Cholesky-factored.

- U:

  Numeric matrix of shape `[n, k]`; low-rank left factor.

- V:

  Numeric matrix of shape `[k, n]`; low-rank right factor. Defaults to
  `t(U)` (symmetric update).

- C_inv:

  Numeric matrix of shape `[k, k]`; inverse of the central factor \\C\\.
  Defaults to `diag(k)`.

## Value

A length-1 numeric: \\\log\|A + UCV\|\\.

## See also

[`woodbury_solve`](https://bbuchsbaum.github.io/amatrix/reference/woodbury_solve.md),
[`chol_factor`](https://bbuchsbaum.github.io/amatrix/reference/chol_factor.md)

## Examples

``` r
A <- crossprod(matrix(rnorm(25), 5)) + diag(5)
U <- matrix(rnorm(10), 5, 2)
ld <- woodbury_logdet(A, U)
is.finite(ld)
#> [1] TRUE
```
