# Row and column summary methods for adgCMatrix

Compute row or column sums and means for a sparse `adgCMatrix`,
dispatching through the amatrix backend when GPU acceleration is
available.

## Usage

``` r
# S4 method for class 'adgCMatrix'
rowSums(x, na.rm = FALSE, dims = 1L)

# S4 method for class 'adgCMatrix'
colSums(x, na.rm = FALSE, dims = 1L)

# S4 method for class 'adgCMatrix'
rowMeans(x, na.rm = FALSE, dims = 1L)

# S4 method for class 'adgCMatrix'
colMeans(x, na.rm = FALSE, dims = 1L)
```

## Arguments

- x:

  An `adgCMatrix`.

- na.rm:

  Logical; if `TRUE`, `NA` values are ignored.

- dims:

  Integer; dimensions to sum over (passed to the backend).

## Value

A numeric vector of length equal to the number of rows or columns.

## Examples

``` r
sp <- as(matrix(c(1, 0, 2, 0, 3, 0), 2, 3), "dgCMatrix")
A  <- adgCMatrix(sp)
rowSums(A)
#> [1] 6 0
colMeans(A)
#> [1] 0.5 1.0 1.5
```
