# amatrix.demo

**How a package author accelerates an algorithm with `amatrix`, in one code
change.**

This internal demo package implements logistic regression via IRLS — the
algorithm behind `glm(family = binomial())` — and shows the full adoption
path for `amatrix`: write the algorithm once, run it anywhere, accelerate it
by changing only the *input*.

## The algorithm, before and after

A textbook IRLS iteration forms the weighted normal equations:

```r
# before: base R
xtwx <- crossprod(X, w * X)          # X'WX   — the bottleneck
xtwz <- crossprod(X, w * z)          # X'Wz
beta <- solve(xtwx, xtwz)
```

The `amatrix.demo` version replaces the two hot products with `amatrix`
kernels (the "hot-kernel swap" from
`planning_docs/package-author-migration.md`):

```r
# after: amatrix kernels — see R/logit.R
xtwx <- amatrix::crossprod_weighted(X, w)
xtwz <- amatrix::xty_weighted(X, w, z)
beta <- solve(xtwx, xtwz)
```

That is the entire integration. There is no GPU code, no backend detection,
no branching in the package. The kernels accept plain base matrices and
compute on the CPU reference path.

One honest caveat: this forms the normal equations, which squares the
condition number of the weighted design; `glm.fit` instead solves the
weighted least-squares problem by pivoted QR, which also handles aliased
columns. For well-conditioned designs (the common case, and the fast case)
the results agree to high precision — the tests verify `1e-6` agreement —
and the demo fails loudly, not silently, on singular systems.

## Acceleration is the caller's one-liner

```r
library(amatrix.demo)

fit_cpu <- logit_fit(X, y)                    # plain matrix: CPU

X_gpu   <- amatrix::adgeMatrix(X, mode = "fast")  # wrap once at the boundary
fit_gpu <- logit_fit(X_gpu, y)                # same function: GPU
```

`mode = "fast"` picks a healthy installed backend (MLX on Apple Silicon,
OpenCL/ArrayFire elsewhere) and falls back to CPU honestly when none is
present — `logit_fit()` never needs to know.

## See it run

```r
logit_demo_benchmark()          # times CPU vs every available backend
# logit_fit() on n=20000, p=400 (3 reps, median)
#
#   base matrix (cpu):          736.0 ms
#   adgeMatrix (mlx):           247.0 ms   (3.0x vs cpu, max |coef diff| 1.2e-06)
#
#   Same logit_fit() code in every row — only the input type changed.
```

(Measured on an Apple Silicon laptop, 2026-07; run it on your hardware. On
GPU backends `precision = "fast"` computes in float32, hence coefficient
agreement at ~1e-6 rather than 1e-10.)

Note the comparison is `logit_fit()` against itself on the two input types,
never against `glm()`: normal equations are already algorithmically cheaper
per iteration than glm's pivoted QR, and that difference is not amatrix's
doing. `glm.fit` serves only as the correctness oracle in the tests, where
both solvers converge to the same MLE.

## Correctness is tested, not asserted

`tests/testthat/test-logit.R` verifies:

- coefficients and deviance match `stats::glm.fit()` to `1e-6`;
- a CPU `adgeMatrix` input reproduces the base-matrix result to `1e-10`;
- a `fast`-precision MLX input agrees to float32 tolerance (`1e-3`),
  skipped with a reason when MLX is not available;
- input validation and the intercept-only edge case.

## Running from the monorepo

```r
pkgload::load_all(".")                     # core amatrix, repo root
pkgload::load_all("demopkg/amatrix.demo")
logit_demo_benchmark()
```
