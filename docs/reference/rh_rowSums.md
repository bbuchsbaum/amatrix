# Row sums of a GPU-resident handle

Computes row sums of the matrix stored in a `resident_handle`, using a
GPU-resident reduction when the backend supports it to avoid a
round-trip download. Falls back to
[`base::rowSums`](https://rdrr.io/r/base/colSums.html) on the
materialized matrix when no resident reduction is available.

## Usage

``` r
rh_rowSums(h)
```

## Arguments

- h:

  A `resident_handle`.

## Value

Numeric vector of length `nrow(h)`.

## See also

[`rh_colSums`](https://bbuchsbaum.github.io/amatrix/reference/rh_colSums.md),
[`am_sweep_inplace`](https://bbuchsbaum.github.io/amatrix/reference/am_sweep_inplace.md)

## Examples

``` r
# \donttest{
# requires a backend with residency support (e.g. MLX, OpenCL)
# }
```
