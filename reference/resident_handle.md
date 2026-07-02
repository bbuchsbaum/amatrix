# Create a mutable GPU-resident handle

Wraps an `adgeMatrix` or plain matrix in a lightweight mutable
environment that holds a GPU-resident buffer key. Unlike `adgeMatrix`,
the handle can be updated in place, making it suitable for iterative
algorithms that would otherwise incur per-step S4 object allocation
overhead. The handle owns its resident key and releases the device
buffer when garbage collected.

## Usage

``` r
resident_handle(x, backend = NULL)
```

## Arguments

- x:

  An `adgeMatrix` or plain `matrix`. If `x` is already GPU-resident on
  `backend`, the existing device buffer is reused without re-uploading.

- backend:

  Character string. Name of the backend to use. Defaults to
  `x@preferred_backend` for `adgeMatrix` inputs and `"cpu"` for plain
  matrices. The backend must support GPU residency.

## Value

A `resident_handle` environment with fields `backend_name`,
`resident_key`, `dim`, `dimnames`, `policy`, `precision`, and `active`.

## See also

[`am_sweep_inplace`](https://bbuchsbaum.github.io/amatrix/reference/am_sweep_inplace.md),
[`rh_rowSums`](https://bbuchsbaum.github.io/amatrix/reference/rh_rowSums.md),
[`rh_colSums`](https://bbuchsbaum.github.io/amatrix/reference/rh_colSums.md)

## Examples

``` r
# \donttest{
m <- adgeMatrix(matrix(runif(12), 3, 4), preferred_backend = "cpu")
# resident_handle requires a backend with residency support (e.g. MLX, OpenCL)
# }
```
