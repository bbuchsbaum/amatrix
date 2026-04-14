# GPU-native truncated SVD via Lanczos bidiagonalization

Implements Golub-Kahan Lanczos bidiagonalization directly in ArrayFire
C, keeping all matvecs and CGS2 reorthogonalization on the GPU. Only
2\*work scalars and the final basis matrices cross PCIe per restart; no
per-step host transfers.

## Usage

``` r
irlba_native(
  A,
  nv = 5L,
  nu = nv,
  tol = sqrt(.Machine$double.eps),
  maxit = 100L,
  work = max(nv + 20L, 3L * nv),
  v0 = NULL,
  mode = "fast",
  backend = NULL
)
```

## Arguments

- A:

  A matrix or `adgeMatrix`. Coerced if necessary.

- nv:

  Number of singular values/vectors to compute.

- nu:

  Number of left singular vectors (default = `nv`).

- tol:

  Convergence tolerance.

- maxit:

  Maximum number of restarts.

- work:

  Size of the Lanczos subspace per restart. Larger values converge in
  fewer restarts at the cost of more memory and work per restart.
  Default is `max(nv + 20L, 3L * nv)`.

- v0:

  Optional starting vector (length `ncol(A)`).

- mode, backend:

  Passed to
  [`adgeMatrix()`](https://bbuchsbaum.github.io/amatrix/reference/adgeMatrix.md)
  when coercing.

## Value

A list with components `d`, `u`, `v`, `iter`, `mprod`, compatible with
[`irlba::irlba()`](https://rdrr.io/pkg/irlba/man/irlba.html).

## Details

Compared to `am_irlba`, which routes each Lanczos matvec through S4
dispatch, this function:

- eliminates S4 overhead on the hot path

- replaces k sequential GEMVs for reorthogonalization with one GEMM

- uploads A once and never re-uploads it across restarts

## See also

[`irlba`](https://bbuchsbaum.github.io/amatrix/reference/irlba.md),
[`adgeMatrix`](https://bbuchsbaum.github.io/amatrix/reference/adgeMatrix.md)
