# Coerce amatrix objects to base R types

Convert `adgeMatrix`, `adgCMatrix`, or `aTransposeView` objects to base
R `matrix`, numeric vector, or array by materializing the host copy.

## Usage

``` r
# S4 method for class 'adgeMatrix'
as.matrix(x, ...)

# S3 method for class 'adgeMatrix'
as.matrix(x, ...)

# S4 method for class 'adgCMatrix'
as.matrix(x, ...)

# S3 method for class 'adgCMatrix'
as.matrix(x, ...)

# S4 method for class 'aTransposeView'
as.matrix(x, ...)

# S3 method for class 'aTransposeView'
as.matrix(x, ...)

# S4 method for class 'amChol'
as.matrix(x, ...)

# S4 method for class 'KronMatrix'
as.matrix(x, ...)

# S4 method for class 'adgeMatrix'
as.numeric(x, ...)

# S4 method for class 'adgeMatrix'
as.vector(x, mode = "any")

# S4 method for class 'adgeMatrix'
as.array(x, ...)

# S4 method for class 'adgCMatrix'
as.array(x, ...)
```

## Arguments

- x:

  An `adgeMatrix`, `adgCMatrix`, or `aTransposeView`.

- ...:

  Further arguments passed to the corresponding base R coercion
  function.

- mode:

  Storage mode string passed to `as.vector`.

## Value

A plain R `matrix`, numeric vector, or `array` containing the
materialized host data.

## Examples

``` r
A <- adgeMatrix(matrix(1:6, 2, 3))
as.matrix(A)
#>      [,1] [,2] [,3]
#> [1,]    1    3    5
#> [2,]    2    4    6
```
