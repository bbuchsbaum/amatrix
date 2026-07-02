# Subset a lazy Kronecker product

Extracts elements of a
[`KronMatrix`](https://bbuchsbaum.github.io/amatrix/reference/KronMatrix-class.md)
using standard matrix indexing (`K[i, j]`, `K[i, ]`, `K[, j]`, or linear
`K[i]`). Without this method `K[i, j]` fails with “object of type 'S4'
is not subsettable”.

## Usage

``` r
# S4 method for class 'KronMatrix,ANY,ANY,ANY'
x[i, j, ..., drop = TRUE]
```

## Arguments

- x:

  A
  [`KronMatrix`](https://bbuchsbaum.github.io/amatrix/reference/KronMatrix-class.md).

- i, j:

  Row and column subscripts, following base matrix semantics.

- ...:

  Unused.

- drop:

  Logical; drop dimensions when the result has a single row or column.
  Defaults to `TRUE`, as for base matrices.

## Value

A numeric vector or matrix, exactly as base matrix indexing of
`as.matrix(x)` would return.

## Details

**Note:** this is a materialize-on-subset implementation. The full \\(mp
\times nq)\\ Kronecker product is formed via
[`kronecker`](https://rdrr.io/r/base/kronecker.html) before indexing, so
subsetting does not preserve the memory advantage of the lazy
representation. It is intended for convenient inspection rather than
large-scale extraction.

## Examples

``` r
K <- kron_matrix(matrix(1:4, 2, 2), diag(2))
K[1, ]
#> [1] 1 0 3 0
K[, 2]
#> [1] 0 1 0 2
K[2, 3]
#> [1] 0
```
