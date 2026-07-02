# Matrix trace

Returns the trace (sum of diagonal elements) of a square matrix or
`aMatrix`.

## Usage

``` r
trace(x)
```

## Arguments

- x:

  A square matrix, sparse `sparseMatrix`, or `aMatrix`.

## Value

A single numeric scalar equal to the sum of diagonal elements.

## Examples

``` r
trace(diag(1:4))
#> [1] 10
```
