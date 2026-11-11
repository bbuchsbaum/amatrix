# Backend-dispatched sweep

Apply a function to each row or column of a matrix, dispatching to the
preferred GPU backend when available.

## Usage

``` r
am_sweep(x, MARGIN, STATS, FUN = "-")
```

## Arguments

- x:

  A numeric matrix or `adgeMatrix`.

- MARGIN:

  1 for rows, 2 for columns.

- STATS:

  Numeric vector of statistics to apply.

- FUN:

  Operation: `"+"`, `"-"` (the default, as for
  [`sweep`](https://rdrr.io/r/base/sweep.html)), `"*"`, or `"/"`.

## Value

A matrix of the same dimensions as `x`.

## See also

[`am_sweep_inplace`](https://bbuchsbaum.github.io/amatrix/reference/am_sweep_inplace.md)

## Examples

``` r
x <- adgeMatrix(matrix(1:6, 2, 3))
# Subtract the column means
am_sweep(x, 2, colMeans(as.matrix(x)), FUN = "-")
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 2 x 3 Matrix of class "adgeMatrix"
#>      [,1] [,2] [,3]
#> [1,] -0.5 -0.5 -0.5
#> [2,]  0.5  0.5  0.5
```
