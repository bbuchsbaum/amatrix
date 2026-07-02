# Return the amatrix backend fallback log

The fallback log records every runtime fall-through from a preferred
backend to the CPU reference path. A non-empty log after a clean
conformance run is a stop-ship condition: it means a backend claimed
support for an op it cannot actually execute, so the result silently
came from a different backend than the one that was requested.

## Usage

``` r
amatrix_fallback_log()
```

## Value

A data.frame with columns `timestamp`, `op`, `from_backend`,
`to_backend`, `reason`. Zero rows means no fallbacks have been recorded.

## See also

[`amatrix_fallback_log_reset`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_fallback_log_reset.md),
[`amatrix_backend_health_probe`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_health_probe.md)

## Examples

``` r
amatrix_fallback_log()
#> [1] timestamp    op           from_backend to_backend   reason      
#> <0 rows> (or 0-length row.names)
amatrix_fallback_log_reset()
```
