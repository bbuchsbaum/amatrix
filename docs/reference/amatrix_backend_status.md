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
#>   name available precision_modes
#> 1  cpu      TRUE     strict,fast
#> 2  mlx     FALSE            fast
#>                                                                                                 features
#> 1                                                         dense_f64,dense_f32,solve,chol,svd,sparse_spmm
#> 2 dense_f32,resident_dense,unified_memory,custom_ops,qr,rsvd,chol_gpu,batched_trsm,eigen_sym,sparse_spmm
#>   residency_capable
#> 1             FALSE
#> 2              TRUE
#>                                                                                                                                                                   capabilities
#> 1                                  matmul,crossprod,tcrossprod,ewise,broadcast_ewise,argmax,scatter_mean,segment_sum,segment_mean,rowSums,colSums,solve,chol,qr,svd,eigen,diag
#> 2 matmul,crossprod,tcrossprod,ewise,broadcast_ewise,argmax,scatter_mean,segment_sum,segment_mean,addmm,rowSums,colSums,qr,svd,rsvd,chol,chol_gpu,batched_trsm,eigen,covariance
amatrix_backend_status("cpu")
#>   name available precision_modes                                       features
#> 1  cpu      TRUE     strict,fast dense_f64,dense_f32,solve,chol,svd,sparse_spmm
#>   residency_capable
#> 1             FALSE
#>                                                                                                                                  capabilities
#> 1 matmul,crossprod,tcrossprod,ewise,broadcast_ewise,argmax,scatter_mean,segment_sum,segment_mean,rowSums,colSums,solve,chol,qr,svd,eigen,diag
```
