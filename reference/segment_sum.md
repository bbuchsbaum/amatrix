# Segment sum by group labels

Sum rows of `x` grouped by integer `labels`, dispatching to GPU when
available.

## Usage

``` r
segment_sum(x, labels, K)
```

## Arguments

- x:

  A numeric matrix or `adgeMatrix`.

- labels:

  Integer vector of group labels (1-based).

- K:

  Number of groups.

## Value

A `K`-by-`ncol(x)` matrix of group sums.

## See also

[`segment_mean`](https://bbuchsbaum.github.io/amatrix/reference/segment_mean.md),
[`am_scatter_mean`](https://bbuchsbaum.github.io/amatrix/reference/am_scatter_mean.md)
