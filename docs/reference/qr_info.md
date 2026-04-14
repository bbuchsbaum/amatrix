# Inspect an amQR factorization object

Returns a named list of metadata fields describing an `amQR` factor
produced by
[`am_qr()`](https://bbuchsbaum.github.io/amatrix/reference/am_qr.md).

## Usage

``` r
qr_info(qr)
```

## Arguments

- qr:

  An `amQR` object.

## Value

A named list with the following elements:

- rank:

  Integer effective rank of the factored matrix.

- dim:

  Integer vector of length 2: `c(nrow, ncol)` of the source matrix.

- thin:

  Logical; `TRUE` when a thin (economy) QR was computed.

- pivoted:

  Logical; `TRUE` when column pivoting was used.

- pivot:

  Integer permutation vector, or `NULL` when unpivoted.

- representation:

  Character string describing the internal storage format.

- backend_ops:

  Character string naming the backend that owns any resident buffers, or
  `NULL`.

- backend:

  Character string: the preferred backend.

- precision:

  Character string: `"strict"` or `"fast"`.

- method:

  Character string: QR algorithm used.

- compact_factor_available:

  Logical.

- compact_factor_source:

  Character string or `NULL`.

- compact_factor_materialized:

  Logical.

- q_materialized:

  Logical.

- r_materialized:

  Logical.

## See also

[`am_qr`](https://bbuchsbaum.github.io/amatrix/reference/am_qr.md),
[`qr_downdate`](https://bbuchsbaum.github.io/amatrix/reference/qr_downdate.md)

## Examples

``` r
X <- adgeMatrix(matrix(rnorm(30), nrow = 6))
qf <- am_qr(X)
info <- qr_info(qf)
info$rank
#> [1] 5
```
