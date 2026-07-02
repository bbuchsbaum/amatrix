# Coerce an object to adgeMatrix

Converts a matrix-like object or a `resident_handle` to an `adgeMatrix`.
When `x` is a `resident_handle`, ownership of the GPU-resident buffer is
transferred to the new `adgeMatrix`.

## Usage

``` r
as_adgeMatrix(
  x,
  mode = NULL,
  backend = NULL,
  preferred_backend = NULL,
  policy = NULL,
  precision = NULL
)
```

## Arguments

- x:

  A `resident_handle`, base R `matrix`, `dgeMatrix`, or any
  `denseMatrix`.

- mode:

  Single string shortcut; see
  [`adgeMatrix`](https://bbuchsbaum.github.io/amatrix/reference/adgeMatrix.md).

- backend:

  Alias for `preferred_backend`.

- preferred_backend:

  Single string; preferred compute backend.

- policy:

  Single string dispatch policy.

- precision:

  Single string; `"strict"` or `"fast"`.

## Value

An `adgeMatrix`.

## Examples

``` r
m <- matrix(1:6, nrow = 2)
A <- as_adgeMatrix(m)
dim(A)
#> [1] 2 3
```
