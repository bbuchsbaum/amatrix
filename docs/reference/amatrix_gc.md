# Free stale GPU residency entries and optionally flush the model cache

Scans the residency registry and removes entries whose backend no longer
reports the associated device buffer as present. Such stale entries
arise when a backend is unloaded or the device is reset between
sessions. Optionally flushes all cached matrix factors (QR, Cholesky,
SVD) from the session model cache.

## Usage

``` r
amatrix_gc(cache = FALSE)
```

## Arguments

- cache:

  Logical. If `TRUE`, also remove all model-cache entries. Default
  `FALSE`.

## Value

Invisibly, a list with two integer elements:

- dead_entries:

  Number of stale residency slots removed.

- cache_entries_cleared:

  Number of model-cache entries flushed (0 when `cache = FALSE`).

## See also

[`amatrix_memory_stats`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_memory_stats.md),
[`amatrix_residency_info`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_residency_info.md)

## Examples

``` r
amatrix_gc()
amatrix_gc(cache = TRUE)
```
