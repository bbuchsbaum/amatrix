# Create a backend-aware sparse matrix

Converts a sparse or dense matrix to an `adgCMatrix` with the specified
backend, policy, and precision. This is the primary user-facing
constructor for sparse amatrix objects.

## Usage

``` r
adgCMatrix(
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

  A `dgCMatrix`, other `sparseMatrix`, or base R `matrix`.

- mode:

  Single string shortcut passed to `.amatrix_resolve_mode()`. Pass
  `NULL` to use the individual arguments. `mode = "fast"` prefers an
  available fast-capable accelerator automatically, with CPU fallback.

- backend:

  Alias for `preferred_backend`; ignored when `preferred_backend` is
  non-`NULL`.

- preferred_backend:

  Single string naming the preferred compute backend.

- policy:

  Single string; one of `"auto"`, `"cpu"`, `"mlx"`, `"metal"`,
  `"arrayfire"`, `"torch"`.

- precision:

  Single string; `"strict"` or `"fast"`.

## Value

An `adgCMatrix` with the data from `x` and the requested backend
metadata.

## Examples

``` r
m <- matrix(c(1, 0, 0, 2), nrow = 2)
S <- adgCMatrix(m)
S
#> An amatrix sparse matrix [cpu|policy=auto|precision=strict]
#> 2 x 2 sparse Matrix of class "adgCMatrix"
#>         
#> [1,] 1 .
#> [2,] . 2
```
