# Cholesky factorization result

Stores the upper-triangular Cholesky factor `R` of a symmetric
positive-definite `adgeMatrix`, as returned by
[`chol_factor`](https://bbuchsbaum.github.io/amatrix/reference/chol_factor.md).
When the factor is resident on a GPU backend the host-side `@factor`
slot may be an empty zero-row matrix; use
[`chol_solve`](https://bbuchsbaum.github.io/amatrix/reference/chol_solve.md)
rather than accessing slots directly.

## Slots

- `factor`:

  Numeric matrix; the upper-triangular factor `R` such that `t(R) %*% R`
  equals the source matrix. May be `matrix(numeric(0), 0, 0)` when the
  factor lives only on the device.

- `factor_obj`:

  Either an `adgeMatrix` holding the GPU-resident factor, or `NULL` for
  CPU-only factors.

- `source_id`:

  Character string; the `object_id` of the source `adgeMatrix`.

- `precision`:

  Character string; `"strict"` or `"fast"`.

- `backend`:

  Character string; the backend that computed the factorization.

## See also

[`chol_factor`](https://bbuchsbaum.github.io/amatrix/reference/chol_factor.md),
[`chol_solve`](https://bbuchsbaum.github.io/amatrix/reference/chol_solve.md),
[`chol_logdet`](https://bbuchsbaum.github.io/amatrix/reference/chol_logdet.md)
