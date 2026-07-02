# Force materialization of an aMatrix to a host Matrix object

Downloads any GPU-resident data and returns a standard `Matrix`-package
object on the host. For `adgeMatrix` inputs the result is a `dgeMatrix`;
for `adgCMatrix` inputs the result is a `dgCMatrix`; for
`aTransposeView` the transposed dense host matrix is returned. Host-only
objects are returned unchanged.

## Usage

``` r
amatrix_materialize_host(x)
```

## Arguments

- x:

  An `aMatrix` object (`adgeMatrix`, `adgCMatrix`, or `aTransposeView`).

## Value

A `dgeMatrix`, `dgCMatrix`, or the original object if no materialization
is needed.

## See also

[`amatrix_residency_info`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_residency_info.md),
[`amatrix_gc`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_gc.md)

## Examples

``` r
m <- adgeMatrix(matrix(1:6, 2, 3))
host <- amatrix_materialize_host(m)
class(host)
#> [1] "dgeMatrix"
#> attr(,"package")
#> [1] "Matrix"
```
