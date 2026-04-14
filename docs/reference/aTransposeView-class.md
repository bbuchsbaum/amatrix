# Lazy transpose view of an adgeMatrix

`aTransposeView` is a zero-copy structural view representing the
transpose of an `adgeMatrix`. It carries no independent dense host
storage; the underlying data lives in `source`. The transposed matrix is
materialized on demand via
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) or
[`amatrix_materialize_host()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_materialize_host.md).

## Slots

- `source`:

  The originating `adgeMatrix`; kept alive by this reference.

- `Dim`:

  Integer vector of length 2 giving the transposed dimensions
  `c(ncol_src, nrow_src)`.

- `Dimnames`:

  List of length 2 with transposed dimnames.
