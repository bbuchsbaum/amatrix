# Element-wise operations

Apply an element-wise arithmetic operation to one or two matrices,
dispatching to the preferred GPU backend when available.

## Usage

``` r
ewise(op, e1, e2 = NULL)
```

## Arguments

- op:

  Character string: `"+"`, `"-"`, `"*"`, or `"/"`.

- e1:

  A numeric matrix or `adgeMatrix`.

- e2:

  A numeric matrix, `adgeMatrix`, or `NULL` for unary ops.

## Value

A matrix of the same dimensions as `e1`.
