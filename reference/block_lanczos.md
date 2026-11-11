# Block Lanczos SVD via block Krylov iteration

Computes a truncated SVD using a block Lanczos bidiagonalization. Each
Krylov step issues one GPU GEMM per block rather than sequential GEMVs,
significantly reducing kernel-launch overhead on accelerated backends.

## Usage

``` r
block_lanczos(
  A,
  nv = 5L,
  nu = nv,
  block_size = NULL,
  n_steps = NULL,
  mode = "fast",
  backend = NULL,
  tol = 1e-05
)

block_svd(
  A,
  k,
  block_size = NULL,
  n_steps = NULL,
  mode = "fast",
  backend = NULL,
  tol = 1e-05
)
```

## Arguments

- A:

  Numeric matrix, `adgeMatrix`, or `adgCMatrix`. Plain matrices are
  coerced to `adgeMatrix` using `mode` and `backend`.

- nv:

  Number of right singular vectors to return.

- nu:

  Number of left singular vectors to return. Defaults to `nv`.

- block_size:

  Integer block width for the Krylov iteration. When `NULL` (default), a
  size is chosen automatically based on `nv` and `nu`.

- n_steps:

  Number of Krylov steps. When `NULL` (default), a starting value is
  chosen automatically and doubled until the residuals meet `tol`. An
  explicit value is used as given.

- mode:

  Execution mode passed to
  [`adgeMatrix()`](https://bbuchsbaum.github.io/amatrix/reference/adgeMatrix.md)
  when coercing plain matrices.

- backend:

  Backend name passed to
  [`adgeMatrix()`](https://bbuchsbaum.github.io/amatrix/reference/adgeMatrix.md)
  when coercing plain matrices. Ignored when `A` is already an
  `aMatrix`.

- tol:

  Convergence tolerance on the residuals \\\\A v_i - \sigma_i u_i\\ /
  \sigma_1\\. A warning of class `amatrix_not_converged` is raised when
  it is not met.

- k:

  Number of singular values/vectors. Alias for `nv = nu = k`.

## Value

A named list with components:

- u:

  Numeric matrix `[m, nu]`: left singular vectors.

- d:

  Numeric vector of length `max(nu, nv)`: singular values in decreasing
  order.

- v:

  Numeric matrix `[n, nv]`: right singular vectors.

- iter:

  Integer: number of Krylov steps performed.

- mprod:

  Integer: total matrix-vector products issued.

- residuals:

  Numeric vector: relative residual of each returned triplet.

- converged:

  Logical: whether all residuals met `tol`.

## See also

[`rsvd`](https://bbuchsbaum.github.io/amatrix/reference/rsvd.md),
`block_svd`

## Examples

``` r
A <- matrix(rnorm(200), nrow = 20)
res <- block_lanczos(A, nv = 3L)
length(res$d)
#> [1] 3
```
