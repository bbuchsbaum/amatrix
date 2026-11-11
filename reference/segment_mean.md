# Segment mean by group labels

Compute the mean of rows of `x` grouped by integer `labels`, dispatching
to GPU when available.

## Usage

``` r
segment_mean(x, labels, K)
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

## See also

[`segment_sum`](https://bbuchsbaum.github.io/amatrix/reference/segment_sum.md),
[`am_scatter_mean`](https://bbuchsbaum.github.io/amatrix/reference/am_scatter_mean.md)

## Examples

``` r
x <- matrix(1:12, 6, 2)
labels <- c(1L, 2L, 1L, 3L, 2L, 1L)
segment_mean(x, labels, K = 3L)
#>          [,1]      [,2]
#> [1,] 3.333333  9.333333
#> [2,] 3.500000  9.500000
#> [3,] 4.000000 10.000000
```
