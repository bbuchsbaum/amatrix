# Symmetrise a matrix

Returns `(x + t(x)) / 2`, enforcing exact symmetry. Handles both dense
`aMatrix` and sparse `adgCMatrix` inputs.

## Usage

``` r
sym(x)
```

## Arguments

- x:

  A square matrix or `aMatrix` object.

## Value

A symmetric matrix or `aMatrix` of the same class and dimensions as `x`.

## Examples

``` r
m <- matrix(c(1, 2, 3, 4), 2, 2)
sym(m)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 2 x 2 Matrix of class "adgeMatrix"
#>      [,1] [,2]
#> [1,]  1.0  2.5
#> [2,]  2.5  4.0
```
