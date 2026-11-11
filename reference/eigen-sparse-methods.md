# Eigendecomposition for adgCMatrix

Compute eigenvalues and eigenvectors of a sparse `adgCMatrix`,
dispatching through the amatrix backend.

## Usage

``` r
# S4 method for class 'adgCMatrix'
eigen(x, symmetric, only.values = FALSE, EISPACK = FALSE)
```

## Arguments

- x:

  An `adgCMatrix`.

- symmetric:

  Logical; whether `x` is symmetric. Auto-detected when missing.

- only.values:

  Logical; if `TRUE` only eigenvalues are returned.

- EISPACK:

  Ignored; retained for signature compatibility.

## Value

A list with components `values` and `vectors`.

## Examples

``` r
S <- Matrix::sparseMatrix(i = c(1, 2, 3, 1, 2), j = c(1, 2, 3, 2, 1),
                          x = c(2, 3, 4, 1, 1))
A <- adgCMatrix(S)
eigen(A, symmetric = TRUE)$values
#> [1] 4.000000 3.618034 1.381966
```
