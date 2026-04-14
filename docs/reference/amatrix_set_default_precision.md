# Set the session-level default precision mode

Sets the precision mode applied to new `aMatrix` objects that do not
specify their own precision. Use `"strict"` for reproducible
double-precision results and `"fast"` for maximum GPU throughput with
single/mixed precision.

## Usage

``` r
amatrix_set_default_precision(precision)
```

## Arguments

- precision:

  Character string. Must be one of `"strict"` or `"fast"`.

## Value

Invisibly, `precision`.

## See also

[`amatrix_default_precision`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_default_precision.md),
[`amatrix_set_default_policy`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_set_default_policy.md)

## Examples

``` r
old <- amatrix_default_precision()
amatrix_set_default_precision("strict")
amatrix_set_default_precision(old) # restore
```
