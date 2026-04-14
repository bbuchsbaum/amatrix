# In-place elementwise operation on a resident handle

Applies an elementwise arithmetic operation between the handle's
resident matrix and either a scalar or another resident handle,
replacing the handle's device buffer with the result.

## Usage

``` r
am_ewise_inplace(h, rhs, op)
```

## Arguments

- h:

  A `resident_handle`.

- rhs:

  A length-1 numeric scalar, or a `resident_handle` with identical
  dimensions to `h`.

- op:

  Character string. Arithmetic operator: `"+"`, `"-"`, `"*"`, or `"/"`.

## Value

`h`, invisibly. The handle is modified in place.

## See also

[`am_sweep_inplace`](https://bbuchsbaum.github.io/amatrix/reference/am_sweep_inplace.md),
[`resident_handle`](https://bbuchsbaum.github.io/amatrix/reference/resident_handle.md)

## Examples

``` r
# \donttest{
# requires a backend with residency support (e.g. MLX, OpenCL)
# }
```
