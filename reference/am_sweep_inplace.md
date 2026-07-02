# In-place broadcast sweep on a resident handle

Applies a row-wise or column-wise arithmetic operation between the
resident matrix and a statistics vector, mutating the handle in place.
Equivalent to `sweep(as.matrix(h), MARGIN, STATS, FUN)` but avoids
downloading the matrix to host.

## Usage

``` r
am_sweep_inplace(h, MARGIN, STATS, FUN = "+")
```

## Arguments

- h:

  A `resident_handle`.

- MARGIN:

  Integer. `1L` to sweep across rows (one value per row), `2L` to sweep
  across columns.

- STATS:

  Numeric vector of length equal to the number of rows or columns
  selected by `MARGIN`.

- FUN:

  Character string. Arithmetic operator to apply: `"+"`, `"-"`, `"*"`,
  or `"/"`. Default `"+"`.

## Value

`h`, invisibly. The handle is modified in place; the underlying device
buffer is replaced with the sweep result.

## See also

[`resident_handle`](https://bbuchsbaum.github.io/amatrix/reference/resident_handle.md),
[`am_ewise_inplace`](https://bbuchsbaum.github.io/amatrix/reference/am_ewise_inplace.md)

## Examples

``` r
# \donttest{
# requires a backend with residency support (e.g. MLX, OpenCL)
# }
```
