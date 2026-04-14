# Compute the dispatch plan for a single operation

Evaluates each candidate backend in preference order and returns a
structured plan describing which backend was chosen and why each
candidate was accepted or rejected. The plan respects GPU residency,
precision compatibility, and calibration thresholds.

## Usage

``` r
amatrix_backend_plan(x, op, y = NULL)
```

## Arguments

- x:

  An `aMatrix` object.

- op:

  Character string naming the operation, e.g. `"matmul"`, `"crossprod"`,
  `"svd"`.

- y:

  Right-hand-side `aMatrix` or `NULL`. Used for binary operations such
  as `"matmul"` to check compatibility and calibration workload.

## Value

A named list with elements:

- op:

  Character. The requested operation.

- pinned_backend:

  Character or `NULL`. Backend to which `x` is currently GPU-resident.

- preferred:

  Character vector. Backends evaluated in order.

- requested_precision:

  Character. Precision mode of `x`.

- chosen:

  Character. Name of the chosen backend.

- chosen_path:

  Character. Either `"resident"` or `"cold"`.

- candidates:

  List of per-candidate evaluation records, each a named list with
  logical flags for `registered`, `available`, `precision_compatible`,
  `supported_cold`, `supported_resident`, `calibration_ok`, `supported`,
  and `chosen`.

## See also

[`amatrix_backend_matrix`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_matrix.md),
[`amatrix_explain`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_explain.md),
[`amatrix_execution_info`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_execution_info.md)

## Examples

``` r
m <- adgeMatrix(matrix(1:6, 2, 3))
amatrix_backend_plan(m, "matmul")
#> $op
#> [1] "matmul"
#> 
#> $pinned_backend
#> NULL
#> 
#> $preferred
#> [1] "cpu"  "auto"
#> 
#> $requested_precision
#> [1] "strict"
#> 
#> $chosen
#> [1] "cpu"
#> 
#> $chosen_path
#> [1] "cold"
#> 
#> $candidates
#> $candidates[[1]]
#> $candidates[[1]]$name
#> [1] "cpu"
#> 
#> $candidates[[1]]$registered
#> [1] TRUE
#> 
#> $candidates[[1]]$capabilities
#>  [1] "matmul"          "crossprod"       "tcrossprod"      "ewise"          
#>  [5] "broadcast_ewise" "argmax"          "scatter_mean"    "segment_sum"    
#>  [9] "segment_mean"    "rowSums"         "colSums"         "solve"          
#> [13] "chol"            "qr"              "svd"             "eigen"          
#> [17] "diag"           
#> 
#> $candidates[[1]]$features
#> [1] "dense_f64"   "dense_f32"   "solve"       "chol"        "svd"        
#> [6] "sparse_spmm"
#> 
#> $candidates[[1]]$precision_modes
#> [1] "strict" "fast"  
#> 
#> $candidates[[1]]$available
#> [1] TRUE
#> 
#> $candidates[[1]]$precision_compatible
#> [1] TRUE
#> 
#> $candidates[[1]]$resident_active
#> [1] FALSE
#> 
#> $candidates[[1]]$supported_cold
#> [1] TRUE
#> 
#> $candidates[[1]]$supported_resident
#> [1] FALSE
#> 
#> $candidates[[1]]$calibration_ok
#> [1] TRUE
#> 
#> $candidates[[1]]$supported
#> [1] TRUE
#> 
#> $candidates[[1]]$chosen_path
#> [1] "cold"
#> 
#> $candidates[[1]]$chosen
#> [1] TRUE
#> 
#> 
#> $candidates[[2]]
#> $candidates[[2]]$name
#> [1] "auto"
#> 
#> $candidates[[2]]$registered
#> [1] FALSE
#> 
#> $candidates[[2]]$capabilities
#> character(0)
#> 
#> $candidates[[2]]$features
#> character(0)
#> 
#> $candidates[[2]]$precision_modes
#> character(0)
#> 
#> $candidates[[2]]$available
#> [1] FALSE
#> 
#> $candidates[[2]]$precision_compatible
#> [1] FALSE
#> 
#> $candidates[[2]]$resident_active
#> [1] FALSE
#> 
#> $candidates[[2]]$supported_cold
#> [1] FALSE
#> 
#> $candidates[[2]]$supported_resident
#> [1] FALSE
#> 
#> $candidates[[2]]$calibration_ok
#> [1] TRUE
#> 
#> $candidates[[2]]$supported
#> [1] FALSE
#> 
#> $candidates[[2]]$chosen_path
#> [1] NA
#> 
#> $candidates[[2]]$chosen
#> [1] FALSE
#> 
#> 
#> 
```
