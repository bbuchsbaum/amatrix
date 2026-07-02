# Kronecker product of backend-aware matrices

S4 methods for [`kronecker`](https://rdrr.io/r/base/kronecker.html) and
the `%x%` operator that keep the result as an amatrix. Without them,
`kronecker(A, B)` and `A %x% B` dispatch to the Matrix methods for the
parent `dgeMatrix` / `dgCMatrix` classes and silently demote to a plain
(non-amatrix) result, discarding backend-dispatch metadata.

## Usage

``` r
# S4 method for class 'adgeMatrix,adgeMatrix'
kronecker(X, Y, FUN = "*", make.dimnames = FALSE, ...)

# S4 method for class 'adgeMatrix,adgCMatrix'
kronecker(X, Y, FUN = "*", make.dimnames = FALSE, ...)

# S4 method for class 'adgCMatrix,adgeMatrix'
kronecker(X, Y, FUN = "*", make.dimnames = FALSE, ...)

# S4 method for class 'adgCMatrix,adgCMatrix'
kronecker(X, Y, FUN = "*", make.dimnames = FALSE, ...)

# S4 method for class 'adgeMatrix,matrix'
kronecker(X, Y, FUN = "*", make.dimnames = FALSE, ...)

# S4 method for class 'matrix,adgeMatrix'
kronecker(X, Y, FUN = "*", make.dimnames = FALSE, ...)

# S4 method for class 'adgCMatrix,matrix'
kronecker(X, Y, FUN = "*", make.dimnames = FALSE, ...)

# S4 method for class 'matrix,adgCMatrix'
kronecker(X, Y, FUN = "*", make.dimnames = FALSE, ...)
```

## Arguments

- X, Y:

  Kronecker factors. At least one is an
  [`aMatrix`](https://bbuchsbaum.github.io/amatrix/reference/aMatrix-class.md)
  subclass; the other may be an amatrix, a base `matrix`, or a Matrix
  object.

- FUN:

  Function (or its name) applied to the outer products; passed to the
  underlying Matrix method. Defaults to `"*"`.

- make.dimnames:

  Logical; construct dimnames from the factors. Passed to the underlying
  method.

- ...:

  Further arguments passed to the underlying method.

## Value

An
[`adgeMatrix`](https://bbuchsbaum.github.io/amatrix/reference/adgeMatrix-class.md)
(dense) or
[`adgCMatrix`](https://bbuchsbaum.github.io/amatrix/reference/adgCMatrix-class.md)
(sparse).

## Details

The product itself is computed by Matrix's own Kronecker methods on the
materialized host contents, so values are identical to
[`base::kronecker()`](https://rdrr.io/r/base/kronecker.html) on the
dense contents. The result is re-wrapped as an
[`adgCMatrix`](https://bbuchsbaum.github.io/amatrix/reference/adgCMatrix-class.md)
when it is sparse and an
[`adgeMatrix`](https://bbuchsbaum.github.io/amatrix/reference/adgeMatrix-class.md)
otherwise, inheriting the preferred backend, policy, and precision of
the first amatrix operand.

## Examples

``` r
A <- adgeMatrix(matrix(1:4, 2, 2))
B <- adgeMatrix(diag(2))
kronecker(A, B)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    3    0
#> [2,]    0    1    0    3
#> [3,]    2    0    4    0
#> [4,]    0    2    0    4
A %x% B
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    3    0
#> [2,]    0    1    0    3
#> [3,]    2    0    4    0
#> [4,]    0    2    0    4
```
