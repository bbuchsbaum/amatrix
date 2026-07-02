# When is amatrix fast?

amatrix exists because some linear-algebra workloads are dramatically
faster on a GPU than on a CPU. It also exists because **most workloads
are not**. This vignette is about telling those two apart on your
machine.

The short version: **a GPU helps when you’re amortizing upload cost over
a large amount of compute, and when the compute is GPU-shaped (dense
matmul, big factorizations, iterative solvers).** For everything else —
small matrices, per-row operations, low arithmetic intensity,
unpredictable shapes — the CPU is either faster or not worth the routing
complexity.

amatrix tries very hard to route to the CPU when the GPU wouldn’t win.
If you’re ever unsure whether a given call is routing where you want,
the dispatcher is introspectable: ask it.

## The speed contract

amatrix makes two promises about performance:

1.  **Never worse than CPU by more than 10% on a calibrated op.** When
    calibration data says a backend is slower than the CPU for a given
    (op, size) cell, the dispatcher routes to the CPU instead. You don’t
    need to guess — you can check.

2.  **Regressions are stop-ship.** The nightly benchmark gate fails when
    any op slows down by more than 20% against the recorded baseline on
    the reference machine, so a performance regression blocks a release
    the same way a correctness failure does.

Neither promise is magic — both require a baseline and a calibration on
your hardware. Running them is a one-time step per machine.

## Calibrate once per machine

Performance is machine-specific. A baseline from a MacBook M3 is not
comparable to one from a Linux GPU box, and the crossover point where a
GPU starts beating a CPU depends on the CPU/GPU you happen to own. The
first time you use amatrix on a machine — or any time the hardware
changes — run:

``` r

library(amatrix)

# Calibrate the currently-registered non-CPU backends.
amatrix_calibrate()
```

This measures every calibration op on every registered backend at a set
of canonical sizes, derives the minimum workload size at which each
backend beats the CPU, and persists the result to the user cache
directory. The cache is tagged with a `sys_hash` covering your OS,
machine, and R version; if you move the cache to a different machine,
amatrix detects the mismatch and refuses to use the stale data.

Inspect the current calibration at any time:

``` r

info <- amatrix_calibration_info()
str(info$thresholds, max.level = 2)
```

## Read the benchmark report

[`amatrix_benchmark_report()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_benchmark_report.md)
surfaces the cold-vs-warm timings from the machine-local baseline plus
the calibrated thresholds in one call:

``` r

rep <- amatrix_benchmark_report()

# Cold vs warm per op × backend, with warm/cold ratio.
head(rep$baseline)
#>            op     size     backend  cold_ms    warm_ms  warm_vs_cold_ratio  speedup_vs_cpu
#> 1  covariance  1024x128  cpu       11.19      11.24     1.00                1
#> 2  matmul      1024x128  cpu        5.17       5.17     1.00                1
#> ...
#> 10 many_lm     1024x128  cpu       10.90       1.13     9.67                1

# Calibrated thresholds: for each (backend, op) pair, the minimum workload
# where the backend beats CPU.
head(rep$calibration)
```

`warm_vs_cold_ratio` is the signal you care about most:

- **≈ 1.0** means the op has no warm-up overhead. First call and tenth
  call cost the same.
- **\> 2.0** means there’s a significant cold-start cost — typically the
  first upload of the matrix to the backend, or the first call pulling
  in a JIT-compiled kernel. The warm runs are fast; the cold run is not.
- For ops like `many_lm`, the ratio is close to 10× because the cache
  stores the QR factorization on the first call and reuses it on
  subsequent calls.

## When the GPU wins

From smallest to largest, here’s what actually happens when you call an
op on an `adgeMatrix`:

1.  **Small matrices (roughly \< 1000 elements).** Upload latency alone
    dwarfs the compute. The dispatcher routes to CPU. A 32×32 matmul on
    MLX is slower than on base BLAS, every time.

2.  **Medium matrices (10³–10⁵ elements).** Depends on the op. Dense
    matmul and factorizations start winning on the GPU around 512×512.
    Per-row reductions (`rowSums`, `colSums`) may still be CPU-bound
    because they’re memory-bound, not compute-bound.

3.  **Large matrices (\> 10⁵ elements, dense).** The GPU wins decisively
    for `%*%`, `crossprod`, `chol`, `svd`, `qr`. Factor caching
    (`chol_factor`, `svd_factor`) amplifies this: the first call pays
    the upload + factorize cost, subsequent calls reuse the cached
    factor at memory speed.

4.  **Very large matrices (memory-pressure regime, 4K×4K+).** The
    baseline tracks these explicitly because cache effects dominate. The
    CPU’s speedup-vs-GPU ratio gets worse as n grows.

5.  **Iterative algorithms** (`rsvd`, `block_lanczos`, `sinkhorn`,
    `irlba`). These do many sequential kernel invocations on resident
    data. If the input is already GPU-resident, they amortize the upload
    cost across hundreds of kernel launches. If the input is cold every
    time, you’re paying for upload over and over.

## When the CPU wins (and you should let it)

- **1×1, 1×n, n×1 matrices.** Always CPU. Degenerate shapes never route
  to a GPU backend; the package’s adversarial-input test suite pins this
  behavior.
- **Element-wise ops on tiny matrices.** The dispatch overhead exceeds
  the compute cost.
- **Ops that allocate more than they compute.** If the GPU spends more
  time copying bytes than running flops, the CPU wins.
- **Backends you don’t have.** Optional backends (`amatrix.mlx`,
  `amatrix.arrayfire`, `amatrix.opencl`, `amatrix.metal`) are separate
  packages; if they’re not installed, the dispatcher routes to CPU. The
  health probe in
  [`amatrix_backend_status()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_status.md)
  tells you what’s registered and which backends passed the canary
  check.

## Ask the dispatcher

If you’re debugging a performance surprise, don’t guess — ask the
planner what it did:

``` r

x <- adgeMatrix(matrix(rnorm(1e6), 1000, 1000))
amatrix_backend_plan(x, "matmul", y = diag(1000))
```

The plan shows which backends were considered, which were rejected and
why, and which one was chosen. The human-readable equivalent:

``` r

amatrix_explain(x %*% diag(1000))
```

## Residency and fallback telemetry

Two other surfaces help you understand runtime behavior:

- [`amatrix_backend_status()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_status.md)
  reports per-backend `health` (`unprobed` / `healthy` /
  `unhealthy:<reason>`), precision modes, feature flags, and
  capabilities. A backend that crashes during a health probe is marked
  `unhealthy` and never routed to again this session.

- [`amatrix_fallback_log()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_fallback_log.md)
  records every time dispatch was forced to fall back from a GPU backend
  to CPU. After a clean conformance run it should be empty. If it isn’t,
  something claimed support for an op it couldn’t actually execute.

``` r

amatrix_backend_status()
amatrix_fallback_log()
```

## Regenerating the baseline

The baseline you compare against is machine-local. For day-to-day
inspection you rarely touch the raw data —
[`amatrix_benchmark_report()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_benchmark_report.md)
reads and summarizes it for you. To *measure* a fresh baseline (a first
run, or after a hardware change), use the benchmark harness that ships
in the package repository,
[`tools/benchmark-regression.R`](https://github.com/bbuchsbaum/amatrix/blob/main/tools/benchmark-regression.R):

``` r

# First run or after a hardware change (writes a new baseline):
# Rscript tools/benchmark-regression.R --update

# Subsequent runs — compare to the saved baseline:
# Rscript tools/benchmark-regression.R
```

Regression against the baseline is informational on your laptop; the
authoritative gate runs in the nightly workflow on a reference machine.

## Honest defaults

amatrix ships with calibrated-conservative defaults: when no calibration
exists, the dispatcher routes small workloads to CPU even if a GPU is
available. This matches the Track 4 Speed Contract’s “never silently
worse than CPU” rule. Once you’ve calibrated, the thresholds sharpen and
the GPU gets its fair share of the work.

The bottom line: **you do not need to guess where amatrix is routing.**
Between
[`amatrix_backend_plan()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_plan.md),
[`amatrix_benchmark_report()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_benchmark_report.md),
and
[`amatrix_fallback_log()`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_fallback_log.md),
the runtime state is always introspectable.
