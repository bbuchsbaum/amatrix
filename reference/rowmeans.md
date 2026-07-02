# Row and column means

Compute row or column means of a matrix or `aMatrix`, dispatching to an
accelerated backend when one is available.

## Usage

``` r
rowmeans(x, na.rm = FALSE)

colmeans(x, na.rm = FALSE)
```

## Arguments

- x:

  A matrix or `aMatrix` object.

- na.rm:

  Logical; if `TRUE`, `NA` values are excluded before averaging. Default
  `FALSE`.

## Value

A numeric vector of length `nrow(x)` (`rowmeans`) or `ncol(x)`
(`colmeans`).

## Examples

``` r
m <- matrix(1:12, 3, 4)
rowmeans(m)
#> [1] 5.5 6.5 7.5
colmeans(m)
#> [1]  2  5  8 11
```
