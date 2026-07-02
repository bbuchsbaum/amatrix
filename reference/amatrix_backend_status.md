# Summarise the status of registered backends

Returns a data.frame with one row per backend describing its
availability, supported precision modes, features, capabilities, and
whether it supports GPU residency.

## Usage

``` r
amatrix_backend_status(names = NULL)
```

## Arguments

- names:

  Character vector of backend names to query. When `NULL` (default) all
  registered backends are included, with optional backends
  auto-registered first if possible.

## Value

A data.frame with columns:

- name:

  Character. Backend identifier.

- available:

  Logical. Whether the backend reports itself as available on this
  machine.

- precision_modes:

  Character. Comma-separated precision modes (`"strict"`, `"fast"`).

- features:

  Character. Comma-separated feature strings.

- residency_capable:

  Logical. Whether the backend supports GPU-resident matrix storage.

- capabilities:

  Character. Comma-separated operation capability strings.

## See also

[`amatrix_backend_names`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_names.md),
[`amatrix_register_backend`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_register_backend.md)

## Examples

``` r
amatrix_backend_status()
#>        name available   health
#> 1 arrayfire     FALSE unprobed
#> 2       cpu      TRUE  healthy
#> 3     metal     FALSE unprobed
#> 4       mlx     FALSE unprobed
#> 5    opencl     FALSE unprobed
#>                                                                                      health_reason
#> 1                            namespace unavailable: there is no package called ‘amatrix.arrayfire’
#> 2                                                                                             <NA>
#> 3                                namespace unavailable: there is no package called ‘amatrix.metal’
#> 4                                  namespace unavailable: there is no package called ‘amatrix.mlx’
#> 5 probe disabled; call amatrix_use_gpu() or set AMATRIX_OPENCL_PROBE_GPU=1 for explicit probe runs
#>   precision_modes                                       features
#> 1            <NA>                                           <NA>
#> 2     strict,fast dense_f64,dense_f32,solve,chol,svd,sparse_spmm
#> 3            <NA>                                           <NA>
#> 4            <NA>                                           <NA>
#> 5            <NA>                                           <NA>
#>   residency_capable
#> 1             FALSE
#> 2             FALSE
#> 3             FALSE
#> 4             FALSE
#> 5             FALSE
#>                                                                                                                                  capabilities
#> 1                                                                                                                                        <NA>
#> 2 matmul,crossprod,tcrossprod,ewise,broadcast_ewise,argmax,scatter_mean,segment_sum,segment_mean,rowSums,colSums,solve,chol,qr,svd,eigen,diag
#> 3                                                                                                                                        <NA>
#> 4                                                                                                                                        <NA>
#> 5                                                                                                                                        <NA>
amatrix_backend_status("cpu")
#>   name available  health health_reason precision_modes
#> 1  cpu      TRUE healthy          <NA>     strict,fast
#>                                         features residency_capable
#> 1 dense_f64,dense_f32,solve,chol,svd,sparse_spmm             FALSE
#>                                                                                                                                  capabilities
#> 1 matmul,crossprod,tcrossprod,ewise,broadcast_ewise,argmax,scatter_mean,segment_sum,segment_mean,rowSums,colSums,solve,chol,qr,svd,eigen,diag
```
