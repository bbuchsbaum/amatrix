# Scatter mean by group labels

Compute the mean of rows of `x` grouped by integer `labels`.

## Usage

``` r
am_scatter_mean(x, labels, K)
```

## Arguments

- x:

  A numeric matrix or `adgeMatrix`.

- labels:

  Integer vector of group labels (1-based).

- K:

  Number of groups.

## Value

A `K`-by-`ncol(x)` matrix of group means.

## Examples

``` r
x <- matrix(1:12, 6, 2)
labels <- c(1L, 2L, 1L, 3L, 2L, 1L)
# One row of column means per cluster
am_scatter_mean(x, labels, K = 3L)
#>          [,1]      [,2]
#> [1,] 3.333333  9.333333
#> [2,] 3.500000  9.500000
#> [3,] 4.000000 10.000000
```
