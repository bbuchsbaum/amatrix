# QR downdate after removing one row

Updates a QR factorization to reflect the removal of a single row from
the original matrix. The current implementation refits from the reduced
matrix; a Givens-rotation path is planned for a future release.

## Usage

``` r
qr_downdate(qr_factor, row_idx, X = NULL)
```

## Arguments

- qr_factor:

  An `amQR` factor object returned by
  [`am_qr()`](https://bbuchsbaum.github.io/amatrix/reference/am_qr.md),
  or any object for which a method is defined.

- row_idx:

  Positive integer index of the row to remove.

- X:

  The original numeric matrix or `adgeMatrix` used to compute
  `qr_factor`. Required for `amQR` factors because they do not store the
  source matrix.

## Value

An updated `amQR` factor with row `row_idx` excluded.

## See also

[`am_qr`](https://bbuchsbaum.github.io/amatrix/reference/am_qr.md),
[`qr_info`](https://bbuchsbaum.github.io/amatrix/reference/qr_info.md),
[`lm_loo_cv`](https://bbuchsbaum.github.io/amatrix/reference/lm_loo_cv.md)

## Examples

``` r
if (FALSE) { # \dontrun{
X <- adgeMatrix(matrix(rnorm(40), nrow = 8))
qf <- am_qr(X)
qf2 <- qr_downdate(qf, row_idx = 3L, X = X)
qr_info(qf2)$dim
} # }
```
