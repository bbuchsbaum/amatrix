# Query GPU residency state of an aMatrix object

Returns a single-row data.frame describing whether `x` is currently
uploaded to a GPU backend and, if so, which backend holds it and whether
that binding is still live (the device buffer still exists).

## Usage

``` r
amatrix_residency_info(x)
```

## Arguments

- x:

  An `aMatrix` object.

## Value

A data.frame with one row and columns:

- backend:

  Character. Backend name, or `NA` when not resident.

- resident_key:

  Character. Internal device buffer key, or `NA`.

- pinned_backend:

  Character. Backend name when the binding is confirmed live, otherwise
  `NA`.

- live:

  Logical. `TRUE` when the backend still holds the buffer identified by
  `resident_key`.

## See also

[`amatrix_materialize_host`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_materialize_host.md),
[`amatrix_memory_stats`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_memory_stats.md)

## Examples

``` r
m <- adgeMatrix(matrix(1:4, 2, 2))
amatrix_residency_info(m)
#>   backend resident_key pinned_backend  live
#> 1    <NA>         <NA>           <NA> FALSE
```
