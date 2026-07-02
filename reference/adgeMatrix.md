# Create a backend-aware dense matrix

Converts a base R matrix or `Matrix::dgeMatrix` to an `adgeMatrix` with
the specified backend, policy, and precision. This is the primary
user-facing constructor for dense amatrix objects.

## Usage

``` r
adgeMatrix(
  x,
  mode = NULL,
  backend = NULL,
  preferred_backend = NULL,
  policy = NULL,
  precision = NULL
)
```

## Arguments

- x:

  A base R `matrix`, `dgeMatrix`, or any `denseMatrix` coercible to
  `dgeMatrix`.

- mode:

  Single string shortcut accepted by `.amatrix_resolve_mode()`; used to
  set `backend`, `policy`, and `precision` together. Pass `NULL` to use
  the individual arguments instead. In particular, `mode = "fast"`
  requests reduced precision and prefers an available fast-capable
  accelerator automatically, with CPU fallback when none is available.

- backend:

  Alias for `preferred_backend`; ignored when `preferred_backend` is
  non-`NULL`.

- preferred_backend:

  Single string naming the preferred compute backend, e.g. `"cpu"`,
  `"mlx"`, or `"metal"`.

- policy:

  Single string; one of `"auto"`, `"cpu"`, `"mlx"`, `"metal"`,
  `"arrayfire"`, `"opencl"`.

- precision:

  Single string; `"strict"` for full double-precision accuracy or
  `"fast"` to allow reduced precision on GPU backends.

## Value

An `adgeMatrix` with the data from `x` and the requested backend
metadata.

## Examples

``` r
m <- matrix(1:6, nrow = 2)
A <- adgeMatrix(m)
A
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 2 x 3 Matrix of class "adgeMatrix"
#>      [,1] [,2] [,3]
#> [1,]    1    3    5
#> [2,]    2    4    6
```
