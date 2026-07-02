# GPU-accelerated truncated SVD via irlba

Wraps `irlba::irlba()` with an `adgeMatrix` input so that every Lanczos
matrix-vector product routes through the amatrix GPU dispatch path. The
matrix `A` is kept resident on device; consecutive matvecs in the
Lanczos loop avoid host round-trips.

## Usage

``` r
irlba(
  A,
  nv = 5,
  nu = nv,
  ...,
  mode = "fast",
  backend = NULL,
  implementation = c("compat", "block"),
  block_size = NULL,
  n_steps = NULL
)
```

## Arguments

- A:

  A matrix, `adgeMatrix`, or `adgCMatrix`. Plain matrices are coerced
  via `adgeMatrix(A, mode=mode, backend=backend)`.

- nv:

  Number of right singular vectors.

- nu:

  Number of left singular vectors. Defaults to `nv`.

- ...:

  Additional arguments forwarded to `irlba::irlba()`. `fastpath` is
  always forced to `FALSE` — the C fastpath bypasses S4 dispatch and
  cannot be GPU-accelerated.

- mode:

  Execution mode passed to
  [`adgeMatrix()`](https://bbuchsbaum.github.io/amatrix/reference/adgeMatrix.md)
  when coercing. `"fast"` permits float32 and enables GPU routing.
  Ignored if `A` is already an `adgeMatrix` or `adgCMatrix`.

- backend:

  Backend name (e.g. `"mlx"`, `"arrayfire"`). Ignored if `A` is already
  an amatrix object.

- implementation:

  Lanczos implementation to use. `"compat"` preserves the current
  `irlba::irlba()` wrapper behavior. `"block"` routes to
  [`block_lanczos`](https://bbuchsbaum.github.io/amatrix/reference/block_lanczos.md)
  for a GEMM-oriented approximation.

- block_size:

  Block size passed to `am_block_lanczos()` when
  `implementation = "block"`. Defaults to a small MLX-friendly block
  size derived from the requested rank.

- n_steps:

  Number of block Krylov steps passed to `am_block_lanczos()` when
  `implementation = "block"`.

## Value

Same structure as `irlba::irlba()`: a list with components `d`, `u`,
`v`, `iter`, `mprod`.

## Details

The hot loop in irlba is two matrix-vector products per Lanczos step:
`A %*% v` and `w %*% A`. Both route through `am_matmul()` when `A` is an
`adgeMatrix`, giving GPU acceleration on the dominant cost.
Orthogonalization, `svd(B)`, and convergence tests remain on CPU where
they belong (the subspace dimension `work` is always small).

Do **not** pass `mult=` — it is deprecated in irlba and forces a
non-standard dispatch path. Pass an `adgeMatrix` instead.

## See also

[`adgeMatrix`](https://bbuchsbaum.github.io/amatrix/reference/adgeMatrix.md),
[`svd_factor`](https://bbuchsbaum.github.io/amatrix/reference/svd_factor.md)
