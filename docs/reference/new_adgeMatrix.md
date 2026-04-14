# Construct an adgeMatrix from a matrix or dgeMatrix

Wraps a base R matrix or `Matrix::dgeMatrix` in an `adgeMatrix`,
attaching backend-dispatch metadata.

## Usage

``` r
new_adgeMatrix(
  x,
  preferred_backend = .amatrix_default_preferred_backend(policy, precision),
  policy = amatrix_default_policy(),
  precision = amatrix_default_precision(),
  src_id = ""
)
```

## Arguments

- x:

  A base R `matrix`, `dgeMatrix`, or any `denseMatrix` coercible to
  `dgeMatrix`.

- preferred_backend:

  Single string; the preferred compute backend. Defaults to `"cpu"`.

- policy:

  Single string; backend dispatch policy. Defaults to
  [`amatrix_default_policy()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_default_policy.md).

- precision:

  Single string; either `"strict"` or `"fast"`. Defaults to
  [`amatrix_default_precision()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_default_precision.md).

- src_id:

  String recording the source object identifier. Pass `""` (default) for
  new objects.

## Value

An `adgeMatrix` object with the same data as `x`.
