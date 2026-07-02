# amatrix

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/bbuchsbaum/amatrix/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/amatrix/actions/workflows/R-CMD-check.yaml)
[![amatrix status badge](https://bbuchsbaum.r-universe.dev/badges/amatrix)](https://bbuchsbaum.r-universe.dev/amatrix)
<!-- badges: end -->

`amatrix` keeps matrix-heavy R code in the `Matrix` idiom while adding backend-aware dispatch and optional accelerator execution. The main target is repeated linear algebra on the same design matrix, especially many-response regression where factorization reuse matters.

## What problem does it solve?

If your workflow already starts with matrices and ends in `%*%`, `crossprod()`, `qr()`, or repeated least-squares fits, `amatrix` gives you:

- Matrix-compatible dense and sparse objects
- predictable CPU fallback when no accelerator backend is available
- optional fast backends such as MLX on Apple Silicon
- shared-design workflows such as `many_lm()`

You do not need to adopt a new tensor framework or rewrite the analysis around a separate object model.

## Quick start

```r
library(amatrix)

set.seed(1)
X <- matrix(rnorm(120 * 6), nrow = 120, ncol = 6)
Y <- matrix(rnorm(120 * 8), nrow = 120, ncol = 8)

# Safe default: strict CPU semantics
X_am <- adgeMatrix(X)
fit <- many_lm(X_am, Y, method = "qr", cache = TRUE,
               include_residuals = TRUE)

dim(coef(fit))
fit$sigma2[1:3]  # residual variances (needs include_residuals = TRUE)
```

On a machine with an available accelerator backend, the workflow stays the same and only the constructor metadata changes.

```r
X_fast <- adgeMatrix(X, mode = "fast")
fit_fast <- many_lm(X_fast, Y, method = "qr", cache = TRUE)
```

## How does backend selection work?

`amatrix` separates the matrix object from the execution backend. You can inspect what is available and what will be chosen for a specific operation.

```r
amatrix_backend_status()
amatrix_backend_plan(X_am, "qr")
amatrix_explain(X_am, "qr")
```

If no accelerator backend is available, the same code still runs on CPU.

For library code, the default per-object path is to wrap at the boundary with `mode = "fast"` and then keep the rest of the code generic:

```r
X_am <- as_adgeMatrix(X, mode = "fast")
fit <- many_lm(X_am, Y, method = "qr", cache = TRUE)
```

That asks `amatrix` to prefer an available fast-capable accelerator automatically and otherwise fall back to CPU. You do not need to hardcode `"mlx"` or set session-global defaults first.

If a caller wants to flip defaults for one local block instead of one object, use `with_amatrix()`:

```r
with_amatrix(policy = "auto", precision = "fast", {
  X_am <- as_adgeMatrix(X)
  fit <- many_lm(X_am, Y, method = "qr", cache = TRUE)
})
```

Persistent session setters such as `amatrix_set_default_policy()` and `amatrix_set_default_precision()` are still available for power users who genuinely want session-wide behavior, for example from `.Rprofile`.

For especially hot repeated paths, you can still bind resident storage explicitly with `amatrix_bind_resident()`, but that should be an optimization step rather than the default package-author workflow.

## When is it worth using?

- one design matrix and many response columns
- repeated fits where QR or other factors can be reused
- workflows that need Matrix-compatible objects instead of a separate tensor API
- Apple Silicon systems where MLX can accelerate the hot path

## Modes

- `adgeMatrix(x)` uses the conservative default path and is safe on CPU-only machines.
- `adgeMatrix(x, mode = "exact")` pins execution to strict CPU semantics.
- `adgeMatrix(x, mode = "fast")` requests reduced-precision execution and prefers an available fast-capable accelerator automatically, with CPU fallback when none is available.

## Backend tiers

amatrix makes explicit support claims per backend. Tier labels track what the
quality gates have actually proven on the current release, not aspiration.

| Tier             | Backend(s)                     | What it means                                                                                          |
|------------------|--------------------------------|--------------------------------------------------------------------------------------------------------|
| **Authoritative**| `cpu`                          | Reference of record for correctness. Always available. CPU failures are stop-ship.                    |
| **Supported**    | `mlx` (Apple Silicon)          | Default Apple Silicon fast path when installed and healthy. Conformance-green on dense/model workloads. |
| **Explicit probe** | `arrayfire`, `opencl`, `metal` | Installed and loadable, but enabled only on explicit request — one call: `amatrix_use_gpu()`. Supported surface is backend-specific and evidence-backed. |
| **Experimental** | other registered backends      | Limited or local-only coverage. Health probe may route away silently. Not a general beta fast-path claim. |

The authoritative definitions, supported-op subsets, and gate evidence live in
`planning_docs/quality-tracking.md §8` and
`planning_docs/backend-certification.md`. At runtime, the same information
(plus current health state from the first-use canary probe) is available via:

```r
amatrix_backend_status()
```

Optional backends ship as separate packages: `amatrix.mlx`,
`amatrix.arrayfire`, `amatrix.opencl`, `amatrix.metal`. Install only the
ones you need — the core package works on pure CPU with no accelerator
dependencies.

## Installation

Core `amatrix` is CPU-only and has no accelerator dependencies. Install it from
the `bbuchsbaum` R-universe:

```r
install.packages(
  "amatrix",
  repos = c("https://bbuchsbaum.r-universe.dev", getOption("repos"))
)
```

GPU execution comes from optional sister packages — install the one that matches
your platform from the same repository:

```r
# Apple Silicon GPU via MLX
install.packages(
  "amatrix.mlx",
  repos = c("https://bbuchsbaum.r-universe.dev", getOption("repos"))
)
```

| Backend package     | Platform                     | System prerequisite                                                       |
|:--------------------|:-----------------------------|:-------------------------------------------------------------------------|
| `amatrix.mlx`       | macOS arm64 (Apple Silicon)  | Homebrew `mlx-c`, or `MLX_C_PREFIX`; builds a mock bridge without it      |
| `amatrix.metal`     | macOS                        | Xcode Command Line Tools                                                  |
| `amatrix.arrayfire` | any OS with ArrayFire >= 3.8 | ArrayFire runtime (on arm64 macOS it pins to the CPU runtime by default)  |
| `amatrix.opencl`    | unix with OpenCL + CLBlast   | OpenCL driver and CLBlast; enable with `amatrix_use_gpu()`                |

Install only the backends you need; the core package works on pure CPU without
any of them.

## GPU in one line

On **Apple Silicon** the GPU is zero-config: install `amatrix.mlx` and MLX
probing is on by default in every launch mode. The first operation that needs
the GPU probes it in a disposable child process and routes eligible work to MLX
automatically — there is nothing to call.

```r
library(amatrix)
#> amatrix GPU backends: mlx ready (activates on first use). See amatrix_gpu_status().

X <- adgeMatrix(matrix(rnorm(4096 * 4096), 4096), mode = "fast")
Z <- X %*% X   # runs on the GPU
```

**Everywhere else** (OpenCL, ArrayFire, Metal) the backends are opt-in, so turn
one on with a single call that enables, health-checks, and adopts the best
installed backend:

```r
amatrix_use_gpu()
#> amatrix: GPU enabled - opencl backend (float32 'fast' precision, ~1e-4 vs
#> float64; 'strict' float64 stays on CPU). amatrix_gpu_status() for details.
```

To see why you are (or are not) on the GPU, call `amatrix_gpu_status()`. GPU work
runs in float32 (`mode = "fast"`, conformance tolerance ~`1e-4`); strict float64
always stays on the CPU. See `vignette("gpu")` for the full story.

## Start here

- `vignette("amatrix")` for the getting-started workflow
- `vignette("gpu")` for going from install to a matmul on the GPU
- `vignette("performance")` for when amatrix is fast (and when it isn't)
- `?adgeMatrix` for matrix constructors
- `?many_lm` for the flagship shared-design regression path
- `?amatrix_backend_status` for backend availability and health
- `?amatrix_benchmark_report` for cold/warm timings and calibrated thresholds
- `?amatrix_fallback_log` for dispatch fallback telemetry
