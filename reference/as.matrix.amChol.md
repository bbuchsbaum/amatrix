# Materialize the dense upper Cholesky factor of an amChol

Returns the factor as a base R matrix regardless of where it currently
lives: the host copy when present, otherwise the resident/deferred
factor object is downloaded and cached.

## Usage

``` r
# S3 method for class 'amChol'
as.matrix(x, ...)
```

## Arguments

- x:

  An
  [`amChol`](https://bbuchsbaum.github.io/amatrix/reference/amChol-class.md)
  object.

- ...:

  Ignored.

## Value

A base R numeric matrix containing the upper-triangular factor.
