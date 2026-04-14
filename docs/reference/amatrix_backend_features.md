# Query the features of a registered backend

Returns the unique feature strings advertised by the named backend, as
reported by its `features()` function. Features describe optional
capabilities such as sparse residency or deferred execution.

## Usage

``` r
amatrix_backend_features(name)
```

## Arguments

- name:

  Character string. Name of a registered backend.

## Value

Character vector of feature identifiers.

## See also

[`amatrix_backend_capabilities`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_capabilities.md),
[`amatrix_backend_status`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_status.md)

## Examples

``` r
amatrix_backend_features("cpu")
#> [1] "dense_f64"   "dense_f32"   "solve"       "chol"        "svd"        
#> [6] "sparse_spmm"
```
