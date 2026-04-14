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
