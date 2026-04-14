# Get the session-level default precision mode

Returns the precision mode used when constructing new `aMatrix` objects
that do not specify their own precision.

## Usage

``` r
amatrix_default_precision()
```

## Value

Character string, either `"strict"` (double precision) or `"fast"`
(single/mixed precision).

## See also

[`amatrix_set_default_precision`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_set_default_precision.md),
[`amatrix_default_policy`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_default_policy.md)

## Examples

``` r
amatrix_default_precision()
#> [1] "strict"
```
