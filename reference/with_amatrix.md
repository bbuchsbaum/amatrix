# Evaluate code with temporary amatrix defaults

Temporarily overrides the session-default dispatch policy and/or
precision mode for the duration of `code`, then restores the previous
values on exit, even when `code` errors.

## Usage

``` r
with_amatrix(policy = NULL, precision = NULL, code)
```

## Arguments

- policy:

  Optional temporary policy. Must be one of `"auto"`, `"cpu"`, `"mlx"`,
  `"metal"`, or `"arrayfire"`.

- precision:

  Optional temporary precision. Must be either `"strict"` or `"fast"`.

- code:

  Expression to evaluate under the temporary defaults.

## Value

The result of evaluating `code`.

## See also

[`adgeMatrix`](https://bbuchsbaum.github.io/amatrix/reference/adgeMatrix.md),
[`amatrix_set_default_policy`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_set_default_policy.md),
[`amatrix_set_default_precision`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_set_default_precision.md)

## Examples

``` r
with_amatrix(policy = "auto", precision = "fast", {
  adgeMatrix(matrix(1:4, nrow = 2))
})
#> An amatrix dense matrix [cpu|policy=auto|precision=fast]
#> 2 x 2 Matrix of class "adgeMatrix"
#>      [,1] [,2]
#> [1,]    1    3
#> [2,]    2    4
```
