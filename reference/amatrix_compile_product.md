# Compile a reusable matrix-product plan

Prepares a fixed left operand for repeated products, choosing and
binding a resident accelerator backend when beneficial. The returned
object is a callable function, so repeated right-hand sides can be
applied without rethinking backend selection each time.

## Usage

``` r
amatrix_compile_product(
  x,
  op = c("matmul", "crossprod", "tcrossprod"),
  backend = "auto",
  precision = amatrix_default_precision(),
  policy = amatrix_default_policy()
)
```

## Arguments

- x:

  Fixed left operand.

- op:

  Product primitive: `"matmul"`, `"crossprod"`, or `"tcrossprod"`.

- backend:

  Backend name or `"auto"`.

- precision:

  Precision to use when wrapping base matrices.

- policy:

  Policy to use when wrapping base matrices.

## Value

A callable `am_product_plan`. Calling the plan with
`materialize = "matrix"` returns a base matrix directly, which is useful
for internal algorithms that immediately need host matrix data.

## Examples

``` r
x <- adgeMatrix(matrix(rnorm(20), 5, 4))
plan <- amatrix_compile_product(x, op = "matmul")
# Reuse the plan for several right-hand sides
plan(matrix(rnorm(8), 4, 2))
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 5 x 2 Matrix of class "adgeMatrix"
#>            [,1]       [,2]
#> [1,] -2.7054843  2.9514750
#> [2,]  2.5681976 -1.9322226
#> [3,] -0.3880687 -1.9392116
#> [4,]  0.6225021 -0.6963171
#> [5,]  3.9459324 -3.4105414
plan(matrix(rnorm(4), 4, 1))
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 5 x 1 Matrix of class "adgeMatrix"
#>             [,1]
#> [1,]  0.17194633
#> [2,]  0.07256624
#> [3,] -0.21162045
#> [4,] -0.02776609
#> [5,]  0.05398815
```
