# Get or set the model cache maximum size

`amatrix_cache_max_size` returns the current limit.
`amatrix_set_cache_max_size` changes the limit and immediately evicts
the least-recently-used entries if the cache exceeds the new bound. When
`max_size` is `Inf` (the default) the cache grows without bound.

## Usage

``` r
amatrix_cache_max_size()

amatrix_set_cache_max_size(max_size)
```

## Arguments

- max_size:

  Positive numeric scalar or `Inf`; the maximum number of factorizations
  to retain in the model cache.

## Value

`amatrix_cache_max_size` returns a length-1 numeric giving the current
limit. `amatrix_set_cache_max_size` returns the new limit invisibly.

## Examples

``` r
old <- amatrix_cache_max_size()
amatrix_set_cache_max_size(10)
amatrix_cache_max_size()
#> [1] 10
amatrix_set_cache_max_size(old)
```
