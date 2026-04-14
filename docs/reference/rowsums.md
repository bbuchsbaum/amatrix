# Row and column sums

Compute row or column sums of a matrix or `aMatrix`, dispatching to an
accelerated backend when one is available.

## Usage

``` r
rowsums(x, na.rm = FALSE, dims = 1L)

colsums(x, na.rm = FALSE, dims = 1L)
```

## Arguments

- x:

  A matrix or `aMatrix` object.

- na.rm:

  Logical; if `TRUE`, missing values are removed before summing. Default
  `FALSE`.

- dims:

  Integer; the number of dimensions to regard as rows (for `rowsums`) or
  columns (for `colsums`). Default `1L`.

## Value

A numeric vector of length `nrow(x)` (`rowsums`) or `ncol(x)`
(`colsums`).

## Examples

``` r
m <- adgeMatrix(matrix(1:12, 3, 4))
rowsums(m)
#> [1] 22 26 30
colsums(m)
#> [1]  6 15 24 33
```
