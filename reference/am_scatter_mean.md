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
