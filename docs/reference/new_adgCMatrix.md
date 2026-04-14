# Construct an adgCMatrix from a sparse or dense matrix

Wraps a `Matrix::dgCMatrix` or any sparse or dense matrix in an
`adgCMatrix`, attaching backend-dispatch metadata.

## Usage

``` r
new_adgCMatrix(
  x,
  preferred_backend = .amatrix_default_preferred_backend(policy, precision),
  policy = amatrix_default_policy(),
  precision = amatrix_default_precision()
)
```

## Arguments

- x:

  A `dgCMatrix`, other `sparseMatrix`, or base R `matrix` to convert.

- preferred_backend:

  Single string; the preferred compute backend. Defaults to `"cpu"`.

- policy:

  Single string; backend dispatch policy.

- precision:

  Single string; `"strict"` or `"fast"`.

## Value

An `adgCMatrix` object.
