# Explain dispatch decisions for an aMatrix operation

Prints a human-readable diagnostic showing which backend was chosen for
`op` on `x`, the accept/reject status of every candidate backend, and
actionable suggestions for improving performance.

## Usage

``` r
amatrix_explain(x, op, y = NULL)
```

## Arguments

- x:

  An `aMatrix` object.

- op:

  Character string. Operation to explain, e.g. `"matmul"`,
  `"crossprod"`, `"svd"`.

- y:

  Right-hand-side `aMatrix` or `NULL`. Supply for binary operations to
  include workload-specific advice.

## Value

Invisibly, the dispatch plan list from
[`amatrix_backend_plan`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_plan.md).
Called primarily for its printed output.

## See also

[`amatrix_backend_plan`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_plan.md),
[`amatrix_backend_matrix`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_matrix.md),
[`amatrix_execution_info`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_execution_info.md)

## Examples

``` r
m <- adgeMatrix(matrix(1:12, 3, 4))
amatrix_explain(m, "matmul")
#> ── amatrix dispatch: matmul ──────────────────────────────────────────── 
#>   object:    adgeMatrix [3 × 4]  precision=strict  preferred=cpu
#>   residency: host (not GPU-resident)
#> 
#> ── candidates ────────────────────────────────────────────────────────── 
#>   ► CHOSEN   cpu            reg avail prec cold calib   [cold]
#>     ......    auto           NO-reg NO-avail NO-prec NO-cold calib 
#> 
#> ── result ────────────────────────────────────────────────────────────── 
#>   chosen: cpu  via cold path (upload + compute)
#> 
#> ──────────────────────────────────────────────────────────────────────── 
```
