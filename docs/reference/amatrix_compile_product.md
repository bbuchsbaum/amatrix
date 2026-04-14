# Compile a Reusable Matrix-Product Plan

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

A callable `am_product_plan`.
