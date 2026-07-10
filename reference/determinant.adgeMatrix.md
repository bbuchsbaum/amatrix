# Determinant for adgeMatrix

Computes the determinant of an `adgeMatrix` with base matrix semantics.
The method materializes the current host value before delegating to
[`determinant`](https://rdrr.io/r/base/det.html) so
[`base::det()`](https://rdrr.io/r/base/det.html) works even when Matrix
is not attached on the search path.

## Usage

``` r
# S3 method for class 'adgeMatrix'
determinant(x, logarithm = TRUE, ...)
```

## Arguments

- x:

  An `adgeMatrix`.

- logarithm:

  Logical; if `TRUE`, return the logarithm of the modulus.

- ...:

  Further arguments passed to
  [`determinant`](https://rdrr.io/r/base/det.html).

## Value

A `"det"` object with `modulus` and `sign`.

## Examples

``` r
A <- adgeMatrix(crossprod(matrix(rnorm(9), 3, 3)) + 3 * diag(3))
det(A)
#> [1] 209.3798
```
