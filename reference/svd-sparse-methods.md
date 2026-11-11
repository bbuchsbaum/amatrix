# Singular value decomposition for adgCMatrix

Compute the singular value decomposition of a sparse `adgCMatrix`,
dispatching through the amatrix backend.

## Usage

``` r
# S4 method for class 'adgCMatrix'
svd(x, nu = min(dim(x)), nv = min(dim(x)), LINPACK = FALSE, ...)
```

## Arguments

- x:

  An `adgCMatrix`.

- nu:

  Number of left singular vectors to compute.

- nv:

  Number of right singular vectors to compute.

- LINPACK:

  Ignored; retained for signature compatibility.

- ...:

  Further arguments passed to the backend.

## Value

A list with components `d`, `u`, and `v`.

## Examples

``` r
S <- Matrix::sparseMatrix(i = c(1, 2, 3, 1), j = c(1, 2, 3, 3),
                          x = c(2, 3, 4, 1))
A <- adgCMatrix(S)
svd(A)$d
#> [1] 4.159415 3.000000 1.923347
```
