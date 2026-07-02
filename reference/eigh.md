# Symmetric eigendecomposition

Computes eigenvalues and eigenvectors of a real symmetric matrix by
dispatching to the active backend via
[`eigen`](https://rdrr.io/r/base/eigen.html) with `symmetric = TRUE`.

## Usage

``` r
eigh(x)
```

## Arguments

- x:

  A real symmetric numeric matrix, `adgeMatrix`, or other object
  accepted by [`eigen`](https://rdrr.io/r/base/eigen.html).

## Value

A list with components `values` (numeric vector of eigenvalues in
ascending order) and `vectors` (numeric matrix whose columns are the
corresponding eigenvectors).

## See also

[`rsvd`](https://bbuchsbaum.github.io/amatrix/reference/rsvd.md)

## Examples

``` r
S <- crossprod(matrix(rnorm(25), nrow = 5))
ev <- eigh(adgeMatrix(S))
length(ev$values)
#> [1] 5
```
