# Coerce an object to adgCMatrix

Converts a sparse or dense matrix-like object to an `adgCMatrix` with
the requested backend metadata.

## Usage

``` r
as_adgCMatrix(
  x,
  mode = NULL,
  backend = NULL,
  preferred_backend = NULL,
  policy = NULL,
  precision = NULL
)
```

## Arguments

- x:

  A `dgCMatrix`, other `sparseMatrix`, or base R `matrix`.

- mode:

  Single string shortcut; see
  [`adgCMatrix`](https://bbuchsbaum.github.io/amatrix/reference/adgCMatrix.md).

- backend:

  Alias for `preferred_backend`.

- preferred_backend:

  Single string; preferred compute backend.

- policy:

  Single string dispatch policy.

- precision:

  Single string; `"strict"` or `"fast"`.

## Value

An `adgCMatrix`.

## Examples

``` r
m <- Matrix::sparseMatrix(i = c(1, 2, 3), j = c(1, 3, 2), x = c(1, 2, 3))
A <- as_adgCMatrix(m)
class(A)
#> [1] "adgCMatrix"
#> attr(,"package")
#> [1] "amatrix"
A %*% c(1, 1, 1)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 3 x 1 Matrix of class "adgeMatrix"
#>      [,1]
#> [1,]    1
#> [2,]    2
#> [3,]    3
```
