# Query the capabilities of a registered backend

Returns the unique capability strings advertised by the named backend,
as reported by its
[`capabilities()`](https://rdrr.io/r/base/capabilities.html) function.

## Usage

``` r
amatrix_backend_capabilities(name)
```

## Arguments

- name:

  Character string. Name of a registered backend.

## Value

Character vector of capability identifiers (e.g. `"matmul"`, `"svd"`).

## See also

[`amatrix_backend_features`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_features.md),
[`amatrix_backend_status`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_status.md)

## Examples

``` r
amatrix_backend_capabilities("cpu")
#>  [1] "matmul"          "crossprod"       "tcrossprod"      "ewise"          
#>  [5] "broadcast_ewise" "argmax"          "scatter_mean"    "segment_sum"    
#>  [9] "segment_mean"    "rowSums"         "colSums"         "solve"          
#> [13] "chol"            "qr"              "svd"             "eigen"          
#> [17] "diag"           
```
