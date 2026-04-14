# Query the precision modes supported by a registered backend

Returns the precision mode strings advertised by the named backend.
Valid values are `"strict"` (double precision) and `"fast"`
(single/mixed precision).

## Usage

``` r
amatrix_backend_precision_modes(name)
```

## Arguments

- name:

  Character string. Name of a registered backend.

## Value

Character vector of precision mode identifiers, a subset of
`c("strict", "fast")`.

## See also

[`amatrix_backend_capabilities`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_capabilities.md),
[`amatrix_backend_status`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_status.md)

## Examples

``` r
amatrix_backend_precision_modes("cpu")
#> [1] "strict" "fast"  
```
