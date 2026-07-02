# Low-level backend dispatch for a single operation

Resolves the best available backend for `op` on `x`, attempts the
GPU-resident path when applicable, and falls back to the cold path
(materializing `x` to host) if needed. If the chosen backend does not
implement `method`, the `fallback` function is called instead.

## Usage

``` r
amatrix_dispatch_op(x, op, method = op, y = NULL, args = list(), fallback)
```

## Arguments

- x:

  An `aMatrix` object.

- op:

  Character string. Operation key used for backend selection (e.g.
  `"matmul"`, `"svd"`).

- method:

  Character string. Name of the backend list element to call. Defaults
  to `op`; override when the backend method name differs from the
  operation key.

- y:

  Right-hand-side `aMatrix` or `NULL`. Passed to the backend method and
  used during backend selection.

- args:

  Named list of additional arguments forwarded to the backend method on
  the cold path.

- fallback:

  Zero-argument function called when the chosen backend does not
  implement `method`.

## Value

The result of the backend method, or the result of `fallback()` if the
method is unavailable.

## See also

[`amatrix_backend_plan`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_plan.md),
[`amatrix_materialize_host`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_materialize_host.md)
