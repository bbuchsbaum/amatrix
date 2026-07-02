# Run a canary health probe against a registered backend

Executes a small matmul round-trip against the named backend and
compares the result to the base R reference. On success the backend is
marked `healthy`; on failure it is marked `unhealthy:<reason>`.
Subsequent calls to
[`amatrix_backend_status()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_status.md)
reflect the recorded health.

## Usage

``` r
amatrix_backend_health_probe(name, tol = NULL)
```

## Arguments

- name:

  Character string. Name of a registered backend.

- tol:

  Numeric. Residual tolerance for the probe, default `1e-8` (float64) or
  `1e-4` (if the backend only supports fast precision).

## Value

Invisibly, the health record as a list with elements `status`, `reason`,
`timestamp`.

## Details

The probe is intentionally tiny (10x10 double-precision matmul) so it
completes in milliseconds even on cold GPU. It is not a benchmark; it is
a liveness check.

## See also

[`amatrix_backend_status`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_status.md),
[`amatrix_fallback_log`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_fallback_log.md)

## Examples

``` r
amatrix_backend_health_probe("cpu")
```
