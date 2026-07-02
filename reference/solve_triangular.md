# Solve a triangular linear system

Solves `R %*% x = B` (or `t(R) %*% x = B` when `lower = TRUE`) for `x`,
where `R` is a triangular matrix. Dispatches to a GPU backend when `R`
is an
[`amChol`](https://bbuchsbaum.github.io/amatrix/reference/amChol-class.md)
or `adgeMatrix` with a live resident key and a capable backend.

## Usage

``` r
solve_triangular(R, B, lower = FALSE)
```

## Arguments

- R:

  An
  [`amChol`](https://bbuchsbaum.github.io/amatrix/reference/amChol-class.md),
  `adgeMatrix`, or numeric matrix holding the triangular factor. Upper
  triangular by default.

- B:

  Numeric vector or matrix; the right-hand side.

- lower:

  Logical scalar; `FALSE` (default) treats `R` as upper triangular,
  `TRUE` treats it as lower triangular.

## Value

Numeric vector or matrix `x` satisfying `R %*% x == B` (or its transpose
variant). Returns a vector when `B` is a vector or single-column matrix.

## Examples

``` r
R <- chol(crossprod(matrix(rnorm(16), 4, 4)) + diag(4))
b <- rnorm(4)
x <- solve_triangular(R, b)
```
