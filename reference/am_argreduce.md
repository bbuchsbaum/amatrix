# Row and column argmax/argmin

Return the index of the maximum or minimum value in each row or column.

## Usage

``` r
am_rowargmax(x)

am_rowargmin(x)

am_colargmax(x)

am_colargmin(x)
```

## Arguments

- x:

  A numeric matrix or `adgeMatrix`.

## Value

An integer vector of indices.

## Examples

``` r
x <- adgeMatrix(matrix(c(3, 1, 2, 9, 7, 8), 2, 3))
am_rowargmax(x)
#> [1] 3 2
am_rowargmin(x)
#> [1] 2 1
am_colargmax(x)
#> [1] 1 2 2
am_colargmin(x)
#> [1] 2 1 1
```
