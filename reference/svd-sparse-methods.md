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
