# Inner product of two vectors or matrices

Computes the element-wise inner product `sum(x * y)`, equivalent to
`as.numeric(t(x) %*% y)` for vectors.

## Usage

``` r
dot(x, y)
```

## Arguments

- x:

  A numeric vector, matrix, or `aMatrix`.

- y:

  A numeric vector, matrix, or `aMatrix` conformable with `x`.

## Value

A single numeric scalar.

## Examples

``` r
dot(1:4, 4:1)
#> [1] 20
```
