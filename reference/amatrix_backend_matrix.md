# Tabulate dispatch plans across multiple operations

Runs
[`amatrix_backend_plan`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_plan.md)
for each requested operation and returns the results as a single
data.frame, one row per operation. Useful for inspecting which backend
will be used across an entire workload.

## Usage

``` r
amatrix_backend_matrix(
  x,
  ops = c("matmul", "crossprod", "tcrossprod", "ewise", "rowSums", "colSums", "solve",
    "chol", "qr", "svd", "eigen", "diag"),
  y_map = list()
)
```

## Arguments

- x:

  An `aMatrix` object.

- ops:

  Character vector of operation names. Defaults to the twelve standard
  operations.

- y_map:

  Named list mapping operation names to right-hand-side objects. Use to
  supply a `y` argument for binary operations such as `"matmul"`.

## Value

A data.frame with one row per operation and columns:

- op:

  Character. Operation name.

- precision:

  Character. Precision mode.

- pinned_backend:

  Character. Backend to which `x` is GPU-resident, or `NA`.

- preferred:

  Character. Preference order string.

- chosen:

  Character. Selected backend.

- chosen_path:

  Character. `"resident"` or `"cold"`.

- resident_reuse:

  Logical. Whether the resident path is active.

- cpu_fallback:

  Logical. Whether CPU was chosen despite not being first preference.

- candidate_summary:

  Character. Compact flag string for all candidates.

## See also

[`amatrix_backend_plan`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_plan.md),
[`amatrix_execution_info`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_execution_info.md)

## Examples

``` r
m <- adgeMatrix(matrix(1:6, 2, 3))
amatrix_backend_matrix(m, ops = c("matmul", "crossprod"))
#>          op precision pinned_backend preferred chosen chosen_path
#> 1    matmul    strict           <NA>       cpu    cpu        cold
#> 2 crossprod    strict           <NA>       cpu    cpu        cold
#>   resident_reuse cpu_fallback candidate_summary
#> 1          FALSE        FALSE    cpu[RAP-C-KSX]
#> 2          FALSE        FALSE    cpu[RAP-C-KSX]
```
