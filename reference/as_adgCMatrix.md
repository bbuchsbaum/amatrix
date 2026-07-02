# Coerce an object to adgCMatrix

Converts a sparse or dense matrix-like object to an `adgCMatrix` with
the requested backend metadata.

## Usage

``` r
as_adgCMatrix(
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

  A `dgCMatrix`, other `sparseMatrix`, or base R `matrix`.

- mode:

  Single string shortcut; see
  [`adgCMatrix`](https://bbuchsbaum.github.io/amatrix/reference/adgCMatrix.md).

- backend:

  Alias for `preferred_backend`.

- preferred_backend:

  Single string; preferred compute backend.

- policy:

  Single string dispatch policy.

- precision:

  Single string; `"strict"` or `"fast"`.

## Value

An `adgCMatrix`.
