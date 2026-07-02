# Row and column summary methods for adgeMatrix

Compute row or column sums and means for an `adgeMatrix`, dispatching
through the amatrix backend when GPU acceleration is available.

## Usage

``` r
# S4 method for class 'adgeMatrix'
rowSums(x, na.rm = FALSE, dims = 1L)

# S4 method for class 'adgeMatrix'
colSums(x, na.rm = FALSE, dims = 1L)

# S4 method for class 'adgeMatrix'
rowMeans(x, na.rm = FALSE, dims = 1L)

# S4 method for class 'adgeMatrix'
colMeans(x, na.rm = FALSE, dims = 1L)
```

## Arguments

- x:

  An `adgeMatrix`.

- na.rm:

  Logical; if `TRUE`, `NA` values are ignored.

- dims:

  Integer; dimensions to sum over (passed to the backend).

## Value

A numeric vector of length equal to the number of rows or columns.

## Examples

``` r
A <- adgeMatrix(matrix(1:12, 3, 4))
rowSums(A)
#> [1] 22 26 30
colMeans(A)
#> [1]  2  5  8 11
```
