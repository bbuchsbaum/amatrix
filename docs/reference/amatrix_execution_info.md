# Collect full dispatch information for an aMatrix object

Returns a snapshot of the dispatch state for an `aMatrix`, including
residency, preferred backend, policy, precision, and the per-operation
dispatch matrix for a set of operations.

## Usage

``` r
amatrix_execution_info(
  x,
  ops = c("matmul", "crossprod", "tcrossprod", "ewise", "rowSums", "colSums"),
  y_map = list()
)
```

## Arguments

- x:

  An `aMatrix` object.

- ops:

  Character vector of operation names to include in the dispatch matrix.
  Default covers the six core operations.

- y_map:

  Named list mapping operation names to right-hand-side objects used
  when planning binary operations such as `"matmul"`.

## Value

A named list with elements:

- object_id:

  Character. Internal object identifier.

- preferred_backend:

  Character. Preferred backend slot value.

- pinned_backend:

  Character or `NULL`. Backend to which the object is currently
  GPU-resident.

- policy:

  Character. Dispatch policy slot value.

- precision:

  Character. Precision mode (`"strict"` or `"fast"`).

- residency:

  data.frame. Output of
  [`amatrix_residency_info`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_residency_info.md).

- plans:

  data.frame. Output of
  [`amatrix_backend_matrix`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_matrix.md).

## See also

[`amatrix_backend_plan`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_plan.md),
[`amatrix_backend_matrix`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_matrix.md),
[`amatrix_explain`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_explain.md)
