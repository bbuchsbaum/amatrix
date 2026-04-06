# amatrix

GPU acceleration for R matrix workloads, with Matrix-compatible objects and no framework lock-in.

## What it is

`amatrix` is an acceleration layer, not a framework. Your existing Matrix-based R code gets faster. You do not adopt a new computational model.

```r
# Before
X <- matrix(design, nrow = n)
fit <- lm.fit(X, y)

# After — one constructor change, same downstream code
library(amatrix)
library(amatrix.models)

X <- adgeMatrix(design, mode = "fast", backend = "mlx")
fit <- many_lm(X, Y_many, method = "qr", cache = TRUE)
```

## What it is not

`amatrix` is not trying to become a bag of `gpuFoo()` wrappers.

The design target is:

- one Matrix-compatible object model
- one backend-planning layer
- one factor-first runtime
- a few benchmark-backed workflows that users actually feel

That means the package should win as infrastructure, not as a collection of one-off GPU demos. A short design note on this point lives in [docs/legacy-r-gpu-lessons.md](/Users/bbuchsbaum/code/amatrix/docs/legacy-r-gpu-lessons.md).

## Why

The primary value is GPU speed on Apple Silicon via [MLX](https://github.com/ml-explore/mlx). The design is optimized for workloads that are common in statistical computing but are not served by existing GPU frameworks:

- one shared design matrix, many response columns
- repeated least-squares fits with the same `X`
- f64-correct numerics (no silent float32 downcast under the default mode)
- Matrix-compatible S4 objects that work with existing R packages

`torch` for R cannot win this space: torch tensors are not Matrix subclasses, torch is float32-first, and torch has no equivalent of a cached shared-X QR factorization across many response columns.

The package is also intentionally conservative about where acceleration is promised:

- large dense repeated workloads are the main target
- benchmark-backed workflow wins matter more than isolated kernel claims
- unsupported or low-value shapes should remain boring CPU fallbacks

## Modes

```r
adgeMatrix(x)                                   # default balanced: strict f64, CPU
adgeMatrix(x, mode = "exact")                   # strict f64, CPU-pinned
adgeMatrix(x, mode = "fast", backend = "mlx")   # fast: f32-permitted, MLX GPU
```

Current contract:

- `exact` is strict CPU semantics.
- `balanced` is currently conservative strict-CPU behavior unless you supply an explicit `backend =`.
- `fast` permits float32-oriented backend execution, but you still need `backend = "mlx"` or another supported backend to leave CPU.

## Flagship workload

One design matrix, many response columns, GPU-accelerated QR:

```r
X <- adgeMatrix(design, mode = "fast", backend = "mlx")
fit <- many_lm(X, Y_many, method = "qr", cache = TRUE)

coef(fit)       # coefficients for all response columns
fit$rss         # residual sum of squares
fit$sigma2      # per-column sigma^2
```

Benchmark (this machine, 1024×128 design):

| RHS columns | CPU cached QR | MLX resident |
|-------------|--------------|--------------|
| 8           | 0.0016 s     | 0.0020 s     |
| 32          | 0.0058 s     | 0.0018 s     |
| 128         | 0.0215 s     | 0.0025 s     |

MLX dominates once the number of right-hand sides grows.

## Packages

| Package | Role |
|---------|------|
| `amatrix` | Core: classes, constructors, backend registry, kernels |
| `amatrix.mlx` | Apple Silicon MLX backend |
| `amatrix.arrayfire` | Portable dense backend (ArrayFire) |
| `amatrix.models` | Model-core fitters: `many_lm`, `array_lm`, `wls_fit`, `covariance` |

## Installation

```r
# Core package (no backend required)
pak::pkg_install("amatrix")

# With MLX backend (Apple Silicon)
pak::pkg_install(c("amatrix", "amatrix.mlx", "amatrix.models"))
```

## CPU mode

`adgeMatrix(x)` without a backend is a valid entry point: it opts into the model surface (`many_lm`, shared-X cache) and makes the object GPU-ready when a backend is added. The package works on CPU. But the reason to adopt it is GPU acceleration.
