# LU factorization result for general square matrices

Stores the original square matrix for use with LAPACK's `DGESV` routine.
Unlike
[`amChol`](https://bbuchsbaum.github.io/amatrix/reference/amChol-class.md),
which caches the explicit triangular factor, `amLU` retains `A` and
delegates factorization to
[`base::solve`](https://rdrr.io/r/base/solve.html) at solve time.

## Slots

- `A`:

  Numeric square matrix; the original matrix passed to
  [`lu_factor`](https://bbuchsbaum.github.io/amatrix/reference/lu_factor.md).

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
