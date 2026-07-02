# QR decomposition of an amatrix object

Computes the QR decomposition of a matrix or `aMatrix`, routing to a
backend-specific implementation when available.

## Usage

``` r
am_qr(x, ...)
```

## Arguments

- x:

  A matrix or `aMatrix` object.

- ...:

  Additional arguments passed to the underlying QR routine.

## Value

An object of class `amDenseQR` (or a wrapped sparse QR for `adgCMatrix`
input) containing the factorisation components.

## Examples

``` r
m <- adgeMatrix(matrix(rnorm(12), 4, 3))
qr_obj <- am_qr(m)
qr.R(qr_obj)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 3 x 3 Matrix of class "adgeMatrix"
#>          [,1]      [,2]      [,3]
#> [1,] 2.822338  1.369303 0.5724748
#> [2,] 0.000000 -1.791752 0.2272415
#> [3,] 0.000000  0.000000 0.6799935
```
