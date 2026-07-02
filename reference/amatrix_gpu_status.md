# GPU backend status: why am I (not) on the GPU?

One row per known GPU backend with the state of every gate between
"installed" and "computing on the GPU": package installed, backend
registered, device available, health, and the registry's recorded reason
when something is off.

## Usage

``` r
amatrix_gpu_status()
```

## Value

A data frame with columns `backend`, `package`, `installed`,
`registered`, `available`, `health`, and `reason`.

## See also

[`amatrix_use_gpu`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_use_gpu.md),
[`amatrix_backend_status`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_status.md),
[`amatrix_explain`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_explain.md)

## Examples

``` r
amatrix_gpu_status()
#>     backend           package installed registered available   health
#> 1       mlx       amatrix.mlx     FALSE      FALSE     FALSE unprobed
#> 2     metal     amatrix.metal     FALSE      FALSE     FALSE unprobed
#> 3 arrayfire amatrix.arrayfire     FALSE      FALSE     FALSE unprobed
#> 4    opencl    amatrix.opencl     FALSE      FALSE     FALSE unprobed
#>                                    reason
#> 1       package amatrix.mlx not installed
#> 2     package amatrix.metal not installed
#> 3 package amatrix.arrayfire not installed
#> 4    package amatrix.opencl not installed
```
