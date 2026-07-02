# Singular value decomposition for adgeMatrix

Compute the singular value decomposition of an `adgeMatrix`, dispatching
through the amatrix backend. A fallback method for plain `matrix`
objects is also provided to preserve base R behaviour after the generic
is promoted to S4.

## Usage

``` r
# S4 method for class 'adgeMatrix'
svd(x, nu = min(dim(x)), nv = min(dim(x)), LINPACK = FALSE, ...)
```

## Arguments

- x:

  An `adgeMatrix` or plain `matrix`.

- nu:

  Number of left singular vectors to compute.

- nv:

  Number of right singular vectors to compute.

- LINPACK:

  Ignored; retained for signature compatibility.

- ...:

  Further arguments passed to the backend.

## Value

A list with components `d` (singular values), `u` (left singular
vectors, `nrow(x)` by `nu`), and `v` (right singular vectors, `ncol(x)`
by `nv`).

## Examples

``` r
A <- adgeMatrix(matrix(rnorm(12), 4, 3))
s <- svd(A)
length(s$d)
#> [1] 3
```
