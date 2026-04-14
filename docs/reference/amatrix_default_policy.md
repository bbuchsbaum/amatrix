# Get the session-level default dispatch policy

Returns the dispatch policy used when an `aMatrix` object does not
specify its own policy. The policy controls which backend is preferred
for operations on new matrices.

## Usage

``` r
amatrix_default_policy()
```

## Value

Character string, one of `"auto"`, `"cpu"`, `"mlx"`, `"metal"`,
`"arrayfire"`, or `"torch"`.

## See also

[`amatrix_set_default_policy`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_set_default_policy.md),
[`amatrix_default_precision`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_default_precision.md)

## Examples

``` r
amatrix_default_policy()
#> [1] "auto"
```
