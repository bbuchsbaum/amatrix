# Warm up GPU backends to eliminate cold-start latency

Pre-compiles GPU kernels by running tiny dummy operations through each
requested backend. Call once before timed work to pay JIT compilation
costs upfront. Errors are silently swallowed; warming never alters
numerical state.

## Usage

``` r
amatrix_warm(
  backend = NULL,
  ops = c("matmul", "crossprod", "qr", "chol"),
  size = c(64L, 64L),
  quiet = FALSE
)
```

## Arguments

- backend:

  Character vector of backend names to warm, or `NULL` to warm all
  non-CPU backends currently registered.

- ops:

  Character vector of operation names to trigger. Recognised values:
  `"matmul"`, `"crossprod"`, `"tcrossprod"`, `"qr"`, `"chol"`, `"svd"`,
  `"solve"`.

- size:

  Integer vector of length 2 giving the dimensions `c(nrow, ncol)` of
  the dummy matrices used during warming.

- quiet:

  Logical; suppress progress messages when `TRUE`.

## Value

An invisible named list, one entry per backend, each a list with
elements `warmed` (logical) and `elapsed_ms` (numeric milliseconds, or
`NA` when unavailable).

## See also

[`amatrix_backend_names`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_names.md)

## Examples

``` r
# \donttest{
results <- amatrix_warm(quiet = TRUE)
# }
```
