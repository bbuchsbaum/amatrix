# Backend-dispatched sweep

Apply a function to each row or column of a matrix, dispatching to the
preferred GPU backend when available.

## Usage

``` r
am_sweep(x, MARGIN, STATS, FUN = "+")
```

## Arguments

- x:

  A numeric matrix or `adgeMatrix`.

- MARGIN:

  1 for rows, 2 for columns.

- STATS:

  Numeric vector of statistics to apply.

- FUN:

  Operation: `"+"`, `"-"`, `"*"`, or `"/"`.

## Value

A matrix of the same dimensions as `x`.

## See also

[`am_sweep_inplace`](https://bbuchsbaum.github.io/amatrix/reference/am_sweep_inplace.md)
