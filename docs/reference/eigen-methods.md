# Eigendecomposition for adgeMatrix

Compute eigenvalues and eigenvectors of an `adgeMatrix`, dispatching
through the amatrix backend for symmetric matrices when GPU acceleration
is available. A fallback method for plain `matrix` preserves base R
behaviour after the generic is promoted to S4.

## Usage

``` r
# S4 method for class 'adgeMatrix'
eigen(x, symmetric, only.values = FALSE, EISPACK = FALSE)
```

## Arguments

- x:

  An `adgeMatrix` or plain `matrix`.

- symmetric:

  Logical indicating whether `x` is symmetric. When missing, symmetry is
  detected automatically from the host copy.

- only.values:

  Logical; if `TRUE` only eigenvalues are returned.

- EISPACK:

  Ignored; retained for signature compatibility.

## Value

A list with components `values` (numeric vector) and `vectors` (matrix,
omitted when `only.values = TRUE`).

## Examples

``` r
S <- adgeMatrix(crossprod(matrix(rnorm(9), 3, 3)))
ev <- eigen(S, symmetric = TRUE)
length(ev$values)
#> [1] 3
```
