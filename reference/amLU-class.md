# LU factorization result for general square matrices

Stores a square matrix together with its pivoted LU factorization
(LAPACK `DGETRF`), computed once by
[`lu_factor`](https://bbuchsbaum.github.io/amatrix/reference/lu_factor.md)
and reused by every
[`lu_solve`](https://bbuchsbaum.github.io/amatrix/reference/lu_solve.md)
call (LAPACK `DGETRS`).

## Slots

- `A`:

  Numeric square matrix; the original matrix passed to
  [`lu_factor`](https://bbuchsbaum.github.io/amatrix/reference/lu_factor.md).

- `lu`:

  The Matrix `denseLU` factorization of `A`, or `NULL` for an empty
  object.

- `source_id`:

  Character string; the `object_id` of the source `adgeMatrix`, or `NA`
  for base matrices.

- `precision`:

  Character string; `"strict"` or `"fast"`, or `NA` for base matrices.

- `backend`:

  Character string; the preferred backend of the source object, or `NA`
  for base matrices.

## See also

[`lu_factor`](https://bbuchsbaum.github.io/amatrix/reference/lu_factor.md),
[`lu_solve`](https://bbuchsbaum.github.io/amatrix/reference/lu_solve.md)
