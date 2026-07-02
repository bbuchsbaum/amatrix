# Internal S3 methods for amatrix helper classes

These S3 methods implement standard base generics (`as.matrix`, `dim`,
`nrow`, `ncol`) for internal amatrix helper classes (`KronMatrix`,
`resident_handle`). They are not part of the public user-facing API —
use the generics directly. This help page exists only to satisfy R CMD
check.

## Usage

``` r
as.matrix.KronMatrix(x, ...)

as.matrix.resident_handle(x, ...)

# S3 method for class 'resident_handle'
dim(x)

nrow.resident_handle(x)

ncol.resident_handle(x)
```

## Arguments

- x:

  A `KronMatrix` or `resident_handle` object.

- ...:

  Additional arguments passed to base methods.

## Value

For `as.matrix` methods: a plain R `matrix`. For `dim`, `nrow`, `ncol`:
an integer (or length-2 integer vector for `dim`).
