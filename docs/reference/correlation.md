# Compute a correlation matrix

Computes the sample correlation matrix of the columns of `X`, optionally
with observation weights and column-blocked covariance accumulation.

## Usage

``` r
correlation(X, center = TRUE, weights = NULL, block_size = NULL)
```

## Arguments

- X:

  Numeric matrix or `adgeMatrix` of shape `[n, p]`.

- center:

  Logical; when `TRUE` (default) column means are subtracted before
  computing covariances.

- weights:

  Numeric vector of length `n`, or `NULL` for unweighted correlation.

- block_size:

  Positive integer or `NULL`. When non-`NULL`, covariances are
  accumulated in blocks of this many columns to limit memory usage.

## Value

An `adgeMatrix` of shape `[p, p]`: the sample correlation matrix with
diagonal entries set to 1.

## See also

[`covariance`](https://bbuchsbaum.github.io/amatrix/reference/covariance.md)

## Examples

``` r
X <- matrix(rnorm(50), nrow = 10)
R <- correlation(X)
round(diag(as.matrix(R)), 6)
#> [1] 1 1 1 1 1
```
