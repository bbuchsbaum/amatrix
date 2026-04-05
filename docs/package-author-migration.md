# Package-Author Migration Notes

This note is for authors who want GPU-capable acceleration from `amatrix` without rewriting their package around a new class system.

The intended promise is:

- small, explicit code changes
- trustworthy CPU fallback
- no backend-specific logic in your high-level code

## Adoption Levels

There are three practical adoption levels.

### 1. Constructor Swap

If your code is already Matrix-aware and generic enough, sometimes the only change needed is to construct `amatrix` objects at the boundary:

```r
library(amatrix)

X <- adgeMatrix(X_host, preferred_backend = "auto", precision = "strict")
Y <- adgeMatrix(Y_host, preferred_backend = "auto", precision = "strict")
```

This is the closest thing to near-zero adaptation. It is worth trying first, but it is not the universal contract.

### 2. Hot-Kernel Swap

If your package has a few clear matrix bottlenecks, replace only those calls.

Before:

```r
XtX <- crossprod(X)
XtY <- crossprod(X, Y)
beta <- solve(XtX, XtY)
```

After:

```r
XtX <- am_crossprod(X)
XtY <- am_crossprod(X, Y)
beta <- am_solve(XtX, XtY)
```

Why this is the preferred package-author path:

- it avoids depending on perfect dispatch
- it keeps backend choice and precision policy centralized
- it leaves the rest of your package structure alone

### 3. Workflow Swap

If the bottleneck is a repeated-fit pattern, adopt a flagship helper instead of manually looping.

Preferred first workload swap:

- `many_lm(..., method = "qr")`

That is currently the strongest public path in the repo for shared-`X`, many-`Y` work.

Before:

```r
fits <- lapply(response_list, function(y) {
  beta <- solve(crossprod(X), crossprod(X, y))
  list(coefficients = beta)
})
```

After:

```r
library(amatrix.models)

fit <- many_lm(X, Y_many, method = "qr", cache = TRUE)
coef(fit)
fit$rss
```

This is the most likely path to memorable speedups because it reuses the shared-`X` cache and keeps the workload on the model-core surface.

If the workload has many right-hand sides, this is also the path where the current MLX resident QR implementation is starting to outperform cached base QR on this machine.

## Similarity Workloads

For covariance or correlation code, the intended migration is equally small.

Before:

```r
S <- stats::cov(X)
R <- stats::cor(X)
```

After:

```r
library(amatrix.models)

S <- covariance(X)
R <- correlation(X)
```

The heavy second-moment step is then routed through `am_crossprod()` internally.

Weighted covariance uses the same shape:

```r
S_w <- covariance(X, weights = w)
```

If the matrix is wide enough that you want to bound the width of each second-moment multiply, use:

```r
S_block <- covariance(X, block_size = 256L)
```

## Precision Policy

Default behavior should remain honest:

```r
X <- adgeMatrix(X_host, precision = "strict")
```

Use `precision = "fast"` only when float32-oriented backend execution is acceptable for the workload.

```r
X <- adgeMatrix(X_host, preferred_backend = "mlx", precision = "fast")
```

Do not write backend-specific branches in your package code unless there is a proven need. Let the backend planner decide.

## What Not To Do

Do not:

- call backend-native bridges directly
- scatter `if (backend == "mlx")` conditionals through your algorithms
- depend on resident execution being present
- assume zero-adaptation is guaranteed for all callers

## Rule Of Thumb

If you can change 3 to 10 lines in a hot path and keep the rest of the code intact, you are using `amatrix` the way it is intended to be adopted.
