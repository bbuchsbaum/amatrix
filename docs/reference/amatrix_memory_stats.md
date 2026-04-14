# Report GPU residency and model cache usage

Returns a snapshot of GPU device memory usage and model cache occupancy
for the current session. Backends that expose a `memory_usage()` method
contribute device byte counts; backends without that method show `NA`
for byte fields.

## Usage

``` r
amatrix_memory_stats()
```

## Value

An object of class `amatrix_memory_stats`, which is a list with two
components:

- residency:

  data.frame with one row per registered backend and columns `backend`
  (character), `resident_objects` (integer count of GPU-resident R
  objects), `bytes_used` (numeric, device bytes in use, or `NA`), and
  `bytes_total` (numeric, total device capacity, or `NA`).

- model_cache:

  List with `n_entries` (integer, number of cached matrix factors) and
  `max_size` (integer or `Inf`, the cache size limit).

## See also

[`amatrix_gc`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_gc.md),
[`amatrix_residency_info`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_residency_info.md)

## Examples

``` r
stats <- amatrix_memory_stats()
print(stats)
#> -- amatrix memory stats ----------------------------------------
#>   model cache: 0 entries (max: unlimited)
#>   residency:
#>     cpu           0 resident object(s)
#>     mlx           0 resident object(s)
#> ----------------------------------------------------------------
```
