# Column sums of a GPU-resident handle

Computes column sums of the matrix stored in a `resident_handle`, using
a GPU-resident reduction when the backend supports it. Falls back to
[`base::colSums`](https://rdrr.io/r/base/colSums.html) on the
materialized matrix when no resident reduction is available.

## Usage

``` r
rh_colSums(h)
```

## Arguments

- h:

  A `resident_handle`.

## Value

Numeric vector of length `ncol(h)`.

## See also

[`rh_rowSums`](https://bbuchsbaum.github.io/amatrix/reference/rh_rowSums.md),
[`am_sweep_inplace`](https://bbuchsbaum.github.io/amatrix/reference/am_sweep_inplace.md)

## Examples

``` r
# \donttest{
# requires a backend with residency support (e.g. MLX, OpenCL)
# }
```
