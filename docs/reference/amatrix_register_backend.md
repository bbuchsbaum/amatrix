# Register a backend with the amatrix dispatch system

Adds a named backend to the session backend registry. The backend must
be a named list containing all required callable fields. Once
registered, the backend is available for dispatch by any `aMatrix`
object whose `preferred_backend` or `policy` slot names it.

## Usage

``` r
amatrix_register_backend(name, backend, overwrite = FALSE)
```

## Arguments

- name:

  Character string. Unique identifier for the backend (e.g. `"mlx"`,
  `"opencl"`).

- backend:

  Named list implementing the backend contract. Required fields:
  `capabilities`, `features`, `precision_modes` (each a zero-argument
  function returning a character vector), `available` (zero-argument
  logical function), `supports`, `matmul`, `crossprod`, `tcrossprod`,
  `ewise`, `rowSums`, `colSums`.

- overwrite:

  Logical. Allow replacement of an existing registration with the same
  `name`. Default `FALSE`.

## Value

Invisibly, `name`.

## See also

[`amatrix_backend_names`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_names.md),
[`amatrix_backend_status`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_status.md)

## Examples

``` r
# Minimal no-op backend for illustration only
noop <- list(
  capabilities   = function() character(),
  features       = function() character(),
  precision_modes = function() "strict",
  available      = function() FALSE,
  supports       = function(op, x, y = NULL) FALSE,
  matmul         = function(x, y) x,
  crossprod      = function(x, y = NULL) x,
  tcrossprod     = function(x, y = NULL) x,
  ewise          = function(x, y, op) x,
  rowSums        = function(x) numeric(nrow(x)),
  colSums        = function(x) numeric(ncol(x))
)
amatrix_register_backend("noop_test", noop, overwrite = TRUE)
```
