# GPU Acceleration: Zero to First Matmul

``` r

library(amatrix)
```

Core `amatrix` is CPU-first and always works with no accelerator
dependencies. GPU execution comes from optional sister packages that
plug in at runtime: `amatrix.mlx` (Apple Silicon via MLX),
`amatrix.metal` (macOS), `amatrix.arrayfire` (ArrayFire), and
`amatrix.opencl` (OpenCL + CLBlast). This vignette takes you from a
fresh install to a matmul running on the GPU, then explains the three
things that surprise people first: precision, routing, and residency.

Because the GPU path depends on hardware that is not present on every
machine (and never on CRAN’s build servers), the GPU-specific chunks
below are shown with `eval = FALSE` and their console output written out
inline. Copy them into a session with a backend installed and they run
as printed.

## How do I install a backend?

The backends ship from the `bbuchsbaum` R-universe. Install core
`amatrix` on its own for pure-CPU work, and add exactly one accelerator
package for your platform.

``` r

# Core package (CPU-only, no accelerator dependencies)
install.packages("amatrix", repos = c("https://bbuchsbaum.r-universe.dev", "https://cloud.r-project.org"))

# Apple Silicon GPU via MLX
install.packages("amatrix.mlx", repos = c("https://bbuchsbaum.r-universe.dev", "https://cloud.r-project.org"))
```

Pick the backend package that matches your machine:

| Backend package | Platform | System prerequisite |
|:---|:---|:---|
| `amatrix.mlx` | macOS arm64 (Apple Silicon) | Homebrew `mlx-c`, or `MLX_C_PREFIX`; builds a mock bridge without it |
| `amatrix.metal` | macOS | Xcode Command Line Tools |
| `amatrix.arrayfire` | any OS with ArrayFire \>= 3.8 | ArrayFire runtime (on arm64 macOS it pins to the CPU runtime by default) |
| `amatrix.opencl` | unix with OpenCL + CLBlast | OpenCL driver and CLBlast; GPU probe is env-gated (see below) |

## Does it just work?

On Apple Silicon, yes. Install `amatrix.mlx` and MLX probing is on by
default in every launch mode – interactive, `Rscript -e`, and
`Rscript file.R`. The first time a session needs the GPU, `amatrix` runs
a one-shot probe in a disposable child process (so a bad driver crashes
the child, never your session) and then routes eligible work to MLX
automatically. There is nothing to call.

You will see a one-line note when the package attaches, confirming what
is available:

``` r

library(amatrix)
#> amatrix GPU backends: mlx ready (activates on first use). See amatrix_gpu_status().
```

That line is suppressible the usual way with
`suppressPackageStartupMessages(library(amatrix))`, and the automatic
probe can be turned off with `options(amatrix.auto_probe = FALSE)` or
the environment variable `AMATRIX_MLX_PROBE_GPU=0`. For a permanently
quiet attach, set `options(amatrix.quiet_startup = TRUE)` in your
`.Rprofile` (or export `AMATRIX_QUIET=1`);
`options(amatrix.optional_backends = FALSE)` goes further and turns the
optional GPU backends off entirely, and per-backend
`options(amatrix.disable_mlx = TRUE)` drops just one. None of this
affects correctness — the note is purely informational and attaching the
package does no GPU probing.

Everywhere else – OpenCL, ArrayFire, Metal – the backends are opt-in, so
you turn them on with one call.
[`amatrix_use_gpu()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_use_gpu.md)
walks the preferred order, enables and health-probes the first installed
backend, sets the fast-precision default, and prints one confirmation
line:

``` r

amatrix_use_gpu()
#> amatrix: GPU enabled - opencl backend (float32 'fast' precision, ~1e-4 vs
#> float64; 'strict' float64 stays on CPU). amatrix_gpu_status() for details.
```

[`amatrix_use_gpu()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_use_gpu.md)
works on Apple Silicon too if you want the explicit confirmation line,
or pass `backend = "mlx"` to force a specific one. It returns the
enabled backend name invisibly, or `FALSE` if none could be brought up.

## Am I actually on the GPU?

Two functions answer this.
[`amatrix_gpu_status()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_gpu_status.md)
is the “why am I (not) on the GPU?” table – one row per backend, walking
every gate between *installed* and *computing on the device*:

``` r

amatrix_gpu_status()
#>     backend       package installed registered available  health
#> 1       mlx   amatrix.mlx      TRUE       TRUE      TRUE healthy
#> 2     metal amatrix.metal     FALSE      FALSE     FALSE unprobed
#> 3 arrayfire amatrix.arr..     FALSE      FALSE     FALSE unprobed
#> 4    opencl amatrix.opencl    FALSE      FALSE     FALSE unprobed
#>                                       reason
#> 1                                       <NA>
#> 2              package amatrix.metal not installed
#> 3          package amatrix.arrayfire not installed
#> 4             package amatrix.opencl not installed
```

Read it left to right: a backend has to be `installed`, then
`registered`, then report `available`, then pass a `health` probe before
work lands on it. The `reason` column names the first gate that failed –
“package not installed”, “opt-in backend; enable with
amatrix_use_gpu()”, or a driver-specific message from the probe.

For a single operation,
[`amatrix_explain()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_explain.md)
shows the actual dispatch decision: which backend was chosen, why every
candidate was accepted or rejected, and what to do about it.

``` r

X <- adgeMatrix(matrix(rnorm(4096 * 4096), 4096), mode = "fast")
amatrix_explain(X, "matmul")
#> ── amatrix dispatch: matmul ─────────────────────────────────────────
#>   object:    adgeMatrix [4096 × 4096]  precision=fast  preferred=mlx
#>   residency: host (not GPU-resident)
#>
#> ── candidates ───────────────────────────────────────────────────────
#>   ► CHOSEN   mlx            reg avail prec cold calib  [cold]
#>     ......   cpu            reg avail prec cold calib
#>
#> ── result ───────────────────────────────────────────────────────────
#>   chosen: mlx  via cold path (upload + compute)
#> ─────────────────────────────────────────────────────────────────────
```

The `mode = "fast"` on the constructor is what makes the object eligible
for the GPU at all – more on that next.
[`amatrix_explain()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_explain.md)
returns the underlying plan invisibly;
`amatrix_backend_plan(X, "matmul")` gives you the same decision as a
plain list if you want to test against it in code.

## What does “fast” precision actually mean?

GPU backends in `amatrix` compute in **float32 only**. That is the
deliberate contract: `"fast"` precision runs on the GPU in single
precision, and `"strict"` precision (float64) **always** computes on the
CPU reference backend, no matter which GPU is enabled. There is no
float64 GPU path to accidentally fall into.

You opt an object into the GPU-eligible path when you build it:

``` r

X_fast <- adgeMatrix(A, mode = "fast")          # per object
amatrix_set_default_precision("fast")            # or session-wide
```

The practical consequence is that GPU results match base R to about
`1e-4` relative error, not to machine epsilon. This is the cross-backend
conformance tolerance the package tests against. So this is expected and
correct:

``` r

A <- matrix(rnorm(2048 * 2048), 2048)
gpu <- as.matrix(adgeMatrix(A, mode = "fast") %*% adgeMatrix(A, mode = "fast"))
cpu <- A %*% A

max(abs(gpu - cpu)) / max(abs(cpu))
#> [1] 7.2e-05

isTRUE(all.equal(gpu, cpu))
#> [1] "Mean relative difference: 3.1e-05"   # float32, as designed

isTRUE(all.equal(gpu, cpu, tolerance = 1e-4))
#> [1] TRUE
```

If you need bit-for-bit agreement with base R, use `mode = "exact"` (or
leave precision at the `"strict"` default) and the work stays on the
CPU. Reach for `"fast"` when a `1e-4` tolerance is acceptable for the
numerics you are doing – which, for most large matmuls, factorizations,
and iterative solvers, it is.

## Why did my small matrix stay on the CPU?

Because that is faster, and `amatrix` knows it. Moving data to the GPU
costs an upload; you only win that cost back when the compute is large
enough to amortize it. Each backend declares calibrated size thresholds
(for example, MLX declines `crossprod` below a 2048 leading dimension)
and simply refuses the op below them, so small matrices route to the
tuned CPU BLAS.

So seeing `cpu` chosen for a small input is the system working, not
failing:

``` r

small <- adgeMatrix(matrix(rnorm(64 * 64), 64), mode = "fast")
amatrix_explain(small, "crossprod")
#> ...
#> ── result ───────────────────────────────────────────────────────────
#>   chosen: cpu  via cold path (upload + compute)
#>
#> ── suggestions ──────────────────────────────────────────────────────
#> * 'mlx' skipped by calibration: 4,096 elements < threshold 4,194,304.
#>   → Use a larger matrix, or re-run amatrix_calibrate() if the
#>     threshold seems wrong for your hardware.
#> ─────────────────────────────────────────────────────────────────────
```

[`amatrix_explain()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_explain.md)
spells out the reason (“skipped by calibration”) so you never have to
guess whether a CPU result means the GPU is broken. It isn’t – the input
was just too small to be worth it. The `performance` vignette
([`vignette("performance")`](https://bbuchsbaum.github.io/amatrix/articles/performance.md))
covers where the crossover lies and how to recalibrate for your
hardware.

## What happens when the GPU fails?

It falls back to the CPU and keeps going. Every runtime GPU failure
degrades to the CPU reference path, marks the backend’s health, and
records an entry you can inspect:

``` r

amatrix_fallback_log()
#>             timestamp     op from_backend to_backend                     reason
#> 1 2026-07-01 10:14:22 matmul          mlx        cpu   kernel launch failed: ...
```

A non-empty fallback log means a backend claimed an op it could not
actually execute – worth reporting, but your results are still correct
because the CPU produced them. Clear it with
[`amatrix_fallback_log_reset()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_fallback_log_reset.md)
when you want to isolate a fresh run.

A more severe failure – a crash during the first-use probe – is
contained in the child process and disables that backend for the rest of
the session with an actionable reason. You will see it in
[`amatrix_gpu_status()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_gpu_status.md):

``` r

amatrix_gpu_status()[, c("backend", "available", "health", "reason")]
#>   backend available   health                              reason
#> 1     mlx     FALSE unhealthy  isolated GPU probe crashed; disabled for session
```

The troubleshooting loop is always the same: read the `reason` column.
It routes you to the fix – install the package, call
[`amatrix_use_gpu()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_use_gpu.md),
set the probe env var, use a bigger matrix, or file the driver crash.

## The one sharp edge: resident handles alias the GPU buffer

Ordinary `amatrix` objects are copy-on-write safe, exactly like base R
matrices. Assigning and then modifying a copy never touches the
original:

``` r

Y <- X
Y[1, 1] <- 999
X[1, 1]        # unchanged
#> [1] 0.418
```

[`resident_handle()`](https://bbuchsbaum.github.io/amatrix/reference/resident_handle.md)
is the deliberate exception. It is a mutable wrapper used by iterative
algorithms to avoid per-step object allocation, and its in-place
operators
([`am_sweep_inplace()`](https://bbuchsbaum.github.io/amatrix/reference/am_sweep_inplace.md),
[`am_ewise_inplace()`](https://bbuchsbaum.github.io/amatrix/reference/am_ewise_inplace.md))
mutate the device buffer directly. When you build a handle from a source
that is **already GPU-resident**, the handle *reuses that same buffer*
rather than re-uploading – so an in-place op through the handle also
mutates the source object and any aliases pointing at that buffer:

``` r

X <- adgeMatrix(A, mode = "fast")
X <- amatrix_bind_resident(X, "mlx")   # X now owns a resident buffer

h <- resident_handle(X)                # reuses X's buffer -- shared, not copied
am_ewise_inplace(h, 2, "*")            # doubles the buffer in place

# X sees the change too -- reference semantics, by design
as.matrix(X)[1, 1] == 2 * A[1, 1]
#> [1] TRUE
```

This is intentional: sharing the buffer is the whole point of a resident
handle on a hot path. If you want an independent copy to mutate, build
the handle from a **non-resident** source (a plain matrix, or an
`adgeMatrix` that has not been bound resident), in which case
[`resident_handle()`](https://bbuchsbaum.github.io/amatrix/reference/resident_handle.md)
uploads a fresh buffer that it owns exclusively. When in doubt, keep
in-place handle work off of objects you still intend to read.

## Where next?

- [`vignette("amatrix")`](https://bbuchsbaum.github.io/amatrix/articles/amatrix.md)
  – the getting-started workflow
  ([`adgeMatrix()`](https://bbuchsbaum.github.io/amatrix/reference/adgeMatrix.md),
  [`many_lm()`](https://bbuchsbaum.github.io/amatrix/reference/many_lm.md))
- [`vignette("performance")`](https://bbuchsbaum.github.io/amatrix/articles/performance.md)
  – when the GPU actually wins, and how to calibrate
- [`?amatrix_use_gpu`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_use_gpu.md)
  and
  [`?amatrix_gpu_status`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_gpu_status.md)
  – enablement and diagnostics
- [`?amatrix_explain`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_explain.md)
  – per-operation dispatch reasoning
- [`?resident_handle`](https://bbuchsbaum.github.io/amatrix/reference/resident_handle.md)
  – the mutable device-resident path for iterative code
