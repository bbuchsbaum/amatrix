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
