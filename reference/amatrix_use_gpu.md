# Enable GPU acceleration for this session

Finds, enables, and health-checks an installed GPU backend, then adopts
it as the session default for `"fast"`-precision work. On Apple Silicon
with amatrix.mlx installed this is usually unnecessary: MLX probing is
on by default and activates on first use. Call this for the opt-in
backends (amatrix.opencl, amatrix.arrayfire, amatrix.metal), to force a
specific backend, or to get an explicit confirmation line.

## Usage

``` r
amatrix_use_gpu(backend = NULL, quiet = FALSE)
```

## Arguments

- backend:

  Optional backend name (`"mlx"`, `"metal"`, `"arrayfire"`, `"opencl"`).
  Default `NULL` tries the automatic preference order and adopts the
  first healthy one.

- quiet:

  Logical; suppress the status messages. Default `FALSE`.

## Value

Invisibly, the name of the enabled backend, or `FALSE` if no GPU backend
could be enabled.

## Details

GPU backends compute in float32 (`"fast"` precision, conformance
tolerance ~1e-4); `"strict"` float64 work always stays on the CPU
reference backend regardless of this setting.

## See also

[`amatrix_gpu_status`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_gpu_status.md),
[`amatrix_backend_status`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_status.md),
[`amatrix_set_default_precision`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_set_default_precision.md)

## Examples

``` r
status <- amatrix_gpu_status()
if (interactive()) amatrix_use_gpu()
```
