# Convert a resident handle back to an adgeMatrix

Materialises the GPU data and creates an adgeMatrix with the resident
key still bound. The handle becomes inert after this call.

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
  lazily.

## Value

An `adgeMatrix`.
