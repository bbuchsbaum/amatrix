# Prepare Operands for a Repeated Matrix Product

Converts inputs to `amatrix` wrappers when needed, chooses a
residency-capable accelerator backend in automatic mode, and binds the
operands so repeated products reuse the resident fast path.

## Usage

``` r
amatrix_prepare_operands(
  x,
  y,
  op = c("matmul", "crossprod", "tcrossprod"),
  backend = "auto",
  precision = amatrix_default_precision(),
  policy = amatrix_default_policy()
)
```

## Arguments

- x:

  Left operand.

- y:

  Right operand.

- op:

  Product primitive: `"matmul"`, `"crossprod"`, or `"tcrossprod"`.

- backend:

  Backend name or `"auto"`.

- precision:

  Precision to use when wrapping base matrices.

- policy:

  Policy to use when wrapping base matrices.

## Value

A list with elements `x`, `y`, and `backend`.
