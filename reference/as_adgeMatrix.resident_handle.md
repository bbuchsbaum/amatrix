# Convert a resident handle back to an adgeMatrix

Creates an adgeMatrix with the resident key still bound. By default the
GPU data is materialized to a host copy. If `defer_host = TRUE`, the
host copy is not materialized until first host access.

## Usage

``` r
as_adgeMatrix.resident_handle(h, ..., defer_host = FALSE)
```

## Arguments

- h:

  A `resident_handle`.

- ...:

  Reserved for future use.

- defer_host:

  When `TRUE`, return a deferred-host `adgeMatrix` that materializes
  lazily. Deferred-host objects are not process-serializable unless
  materialized before persistence; after
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html)/[`readRDS()`](https://rdrr.io/r/base/readRDS.html)
  they fail cleanly instead of returning sentinel data.

## Value

An `adgeMatrix`.
