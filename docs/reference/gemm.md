# Generalised matrix multiply (BLAS DGEMM interface)

Computes `alpha * op(A) %*% op(B) + beta * C`, where `op(X) = t(X)` when
the corresponding `trans` flag is `TRUE`. Routes internally to the most
efficient resident operation for the chosen transpose combination.

## Usage

``` r
gemm(A, B, C = NULL, alpha = 1, beta = 1, transA = FALSE, transB = FALSE)
```

## Arguments

- A:

  A matrix or `aMatrix`.

- B:

  A matrix or `aMatrix`.

- C:

  Optional matrix or `aMatrix` to add after scaling; `NULL` omits the
  addition term.

- alpha:

  Numeric scalar multiplier for `op(A) %*% op(B)`. Default `1.0`.

- beta:

  Numeric scalar multiplier for `C`. Default `1.0`.

- transA:

  Logical; transpose `A` before multiplying. Default `FALSE`.

- transB:

  Logical; transpose `B` before multiplying. Default `FALSE`.

## Value

A matrix of dimensions `nrow(op(A))` by `ncol(op(B))`.

## Examples

``` r
A <- adgeMatrix(matrix(1:6, 2, 3))
B <- adgeMatrix(matrix(1:6, 2, 3))
gemm(A, B, transA = TRUE)          # t(A) %*% B
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 3 x 3 Matrix of class "adgeMatrix"
#>      [,1] [,2] [,3]
#> [1,]    5   11   17
#> [2,]   11   25   39
#> [3,]   17   39   61
gemm(A, B, transB = TRUE)          # A %*% t(B)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 2 x 2 Matrix of class "adgeMatrix"
#>      [,1] [,2]
#> [1,]   35   44
#> [2,]   44   56
```
