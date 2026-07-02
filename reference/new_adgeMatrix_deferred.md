# Construct a deferred adgeMatrix with GPU-only storage

Creates an `adgeMatrix` whose host `@x` slot holds a `NaN` sentinel
vector. The true data lives only on the device until the first host
access, which triggers a transparent download.

## Usage

``` r
new_adgeMatrix_deferred(
  dim,
  dimnames = list(NULL, NULL),
  preferred_backend = .amatrix_default_preferred_backend(policy, precision),
  policy = amatrix_default_policy(),
  precision = amatrix_default_precision(),
  src_id = ""
)
```

## Arguments

- dim:

  Integer vector of length 2 giving `c(nrow, ncol)`.

- dimnames:

  List of length 2 with row and column names, or `list(NULL, NULL)`.

- preferred_backend:

  Single string naming the preferred backend.

- policy:

  Single string; backend dispatch policy.

- precision:

  Single string; `"strict"` or `"fast"`.

- src_id:

  String recording the source object identifier.

## Value

An `adgeMatrix` with `finalizer_env$host_deferred` set to `TRUE`.

## Details

Deferred objects are intentionally not process-serializable: after a
serialization boundary such as
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html)/[`readRDS()`](https://rdrr.io/r/base/readRDS.html)
or
[`serialize()`](https://rdrr.io/r/base/serialize.html)/[`unserialize()`](https://rdrr.io/r/base/serialize.html),
the device resident key is no longer valid unless the host copy was
materialized before persistence. Coercion or printing of such a dead
deferred object fails with a clean error rather than returning sentinel
data.
