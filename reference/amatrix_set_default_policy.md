# Set the session-level default dispatch policy

Sets the dispatch policy applied to new `aMatrix` objects that do not
specify their own policy. The change affects all subsequent matrix
constructions in the current session.

## Usage

``` r
amatrix_set_default_policy(policy)
```

## Arguments

- policy:

  Character string. Must be one of `"auto"`, `"cpu"`, `"mlx"`,
  `"metal"`, `"arrayfire"`, or `"opencl"`.

## Value

Invisibly, `policy`.

## See also

[`amatrix_default_policy`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_default_policy.md),
[`amatrix_set_default_precision`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_set_default_precision.md)

## Examples

``` r
old <- amatrix_default_policy()
amatrix_set_default_policy("auto")
amatrix_set_default_policy(old) # restore
```
