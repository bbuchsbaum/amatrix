# Truncated SVD factorization result

Stores the rank-`k` truncated singular value decomposition of an
`aMatrix` as returned by
[`svd_factor`](https://bbuchsbaum.github.io/amatrix/reference/svd_factor.md).
Left and right singular vectors are kept as base R matrices so they are
always host-accessible; optional `adgeMatrix` copies of `t(u)` and `v`
may be present for GPU-accelerated projection.

## Slots

- `u`:

  Numeric matrix of left singular vectors; `nrow(u)` equals the number
  of rows of the source matrix and `ncol(u)` equals `k`.

- `d`:

  Numeric vector of singular values of length `k`, in descending order.

- `v`:

  Numeric matrix of right singular vectors; `nrow(v)` equals the number
  of columns of the source matrix and `ncol(v)` equals `k`.

- `k`:

  Integer; the requested rank.

- `method`:

  Character string; one of `"exact"`, `"rsvd"`, or `"subspace"`.

- `engine`:

  Character string identifying the low-level solver used (e.g.,
  `"exact_svd"`, `"irlba_svdr"`, `"backend_rsvd"`).

- `source_id`:

  Character string; the `object_id` of the source `aMatrix`.

- `precision`:

  Character string; `"strict"` or `"fast"`.

- `backend`:

  Character string; the backend that computed the SVD.

- `ut_am`:

  Either an `adgeMatrix` holding `t(u)` for GPU matrix-multiply routing,
  or `NULL` on CPU paths.

- `v_am`:

  Either an `adgeMatrix` holding `v` for GPU matrix-multiply routing, or
  `NULL` on CPU paths.

## See also

[`svd_factor`](https://bbuchsbaum.github.io/amatrix/reference/svd_factor.md),
[`svd_project`](https://bbuchsbaum.github.io/amatrix/reference/svd_project.md),
[`svd_reconstruct`](https://bbuchsbaum.github.io/amatrix/reference/svd_reconstruct.md)
