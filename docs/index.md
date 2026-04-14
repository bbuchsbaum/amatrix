# amatrix

`amatrix` keeps matrix-heavy R code in the `Matrix` idiom while adding
backend-aware dispatch and optional accelerator execution. The main
target is repeated linear algebra on the same design matrix, especially
many-response regression where factorization reuse matters.

## What problem does it solve?

If your workflow already starts with matrices and ends in `%*%`,
[`crossprod()`](https://rdrr.io/r/base/crossprod.html),
[`qr()`](https://rdrr.io/r/base/qr.html), or repeated least-squares
fits, `amatrix` gives you:

- Matrix-compatible dense and sparse objects
- predictable CPU fallback when no accelerator backend is available
- optional fast backends such as MLX on Apple Silicon
- shared-design workflows such as
  [`many_lm()`](https://bbuchsbaum.github.io/amatrix/reference/many_lm.md)

You do not need to adopt a new tensor framework or rewrite the analysis
around a separate object model.

## Quick start

``` r
library(amatrix)

set.seed(1)
X <- matrix(rnorm(120 * 6), nrow = 120, ncol = 6)
Y <- matrix(rnorm(120 * 8), nrow = 120, ncol = 8)

# Safe default: strict CPU semantics
X_am <- adgeMatrix(X)
fit <- many_lm(X_am, Y, method = "qr", cache = TRUE)

dim(coef(fit))
fit$sigma2[1:3]
```

On a machine with an available accelerator backend, the workflow stays
the same and only the constructor metadata changes.

``` r
X_fast <- adgeMatrix(X, mode = "fast")
fit_fast <- many_lm(X_fast, Y, method = "qr", cache = TRUE)
```

## How does backend selection work?

`amatrix` separates the matrix object from the execution backend. You
can inspect what is available and what will be chosen for a specific
operation.

``` r
amatrix_backend_status()
amatrix_backend_plan(X_am, "qr")
amatrix_explain(X_am, "qr")
```

If no accelerator backend is available, the same code still runs on CPU.

For library code, the default per-object path is to wrap at the boundary
with `mode = "fast"` and then keep the rest of the code generic:

``` r
X_am <- as_adgeMatrix(X, mode = "fast")
fit <- many_lm(X_am, Y, method = "qr", cache = TRUE)
```

That asks `amatrix` to prefer an available fast-capable accelerator
automatically and otherwise fall back to CPU. You do not need to
hardcode `"mlx"` or set session-global defaults first.

If a caller wants to flip defaults for one local block instead of one
object, use
[`with_amatrix()`](https://bbuchsbaum.github.io/amatrix/reference/with_amatrix.md):

``` r
with_amatrix(policy = "auto", precision = "fast", {
  X_am <- as_adgeMatrix(X)
  fit <- many_lm(X_am, Y, method = "qr", cache = TRUE)
})
```

Persistent session setters such as
[`amatrix_set_default_policy()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_set_default_policy.md)
and
[`amatrix_set_default_precision()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_set_default_precision.md)
are still available for power users who genuinely want session-wide
behavior, for example from `.Rprofile`.

For especially hot repeated paths, you can still bind resident storage
explicitly with
[`amatrix_bind_resident()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_bind_resident.md),
but that should be an optimization step rather than the default
package-author workflow.

## When is it worth using?

- one design matrix and many response columns
- repeated fits where QR or other factors can be reused
- workflows that need Matrix-compatible objects instead of a separate
  tensor API
- Apple Silicon systems where MLX can accelerate the hot path

## Modes

- `adgeMatrix(x)` uses the conservative default path and is safe on
  CPU-only machines.
- `adgeMatrix(x, mode = "exact")` pins execution to strict CPU
  semantics.
- `adgeMatrix(x, mode = "fast")` requests reduced-precision execution
  and prefers an available fast-capable accelerator automatically, with
  CPU fallback when none is available.

## Installation

``` r
# Core package
pak::pkg_install("amatrix")

# Optional MLX backend for Apple Silicon
pak::pkg_install(c("amatrix", "amatrix.mlx"))
```

Other accelerator backends can be installed separately when needed.

## Start here

- [`vignette("amatrix")`](https://bbuchsbaum.github.io/amatrix/articles/amatrix.md)
  for the getting-started workflow
- [`?adgeMatrix`](https://bbuchsbaum.github.io/amatrix/reference/adgeMatrix.md)
  for matrix constructors
- [`?many_lm`](https://bbuchsbaum.github.io/amatrix/reference/many_lm.md)
  for the flagship shared-design regression path
- [`?amatrix_backend_status`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_status.md)
  for backend availability and capability checks
