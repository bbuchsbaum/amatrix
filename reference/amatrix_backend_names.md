# List names of all registered backends

Returns the names of every backend currently in the session registry.
When optional backends are enabled (the default), this also attempts to
auto-register any installed optional backend packages before returning
the list.

## Usage

``` r
amatrix_backend_names()
```

## Value

Character vector of registered backend names, sorted alphabetically.
Always includes at least `"cpu"`.

## See also

[`amatrix_backend_status`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_status.md),
[`amatrix_register_backend`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_register_backend.md)

## Examples

``` r
amatrix_backend_names()
#> [1] "cpu"
```
