# amatrix.models

`amatrix.models` is the model-facing surface built on top of the `amatrix` substrate.

Start here:

```r
library(amatrix)
library(amatrix.models)

X <- adgeMatrix(design, mode = "fast", backend = "mlx")
fit <- many_lm(X, Y_many, method = "qr", cache = TRUE)

coef(fit)
fit$rss
fit$sigma2
```

Why this is the current flagship:

- one shared design matrix `X`
- many response columns `Y`
- QR-backed least squares
- automatic shared-`X` cache reuse
- resident MLX QR path dominates once the number of right-hand sides grows

Without a GPU backend, the package still works and caches the QR factorization of `X` across all response columns — but the primary reason to adopt it is GPU acceleration via `mode = "fast", backend = "mlx"`.

Other exported surfaces:

- `lm_fit(...)` — single least-squares fit
- `wls_fit(...)` — single weighted fit
- `ridge_fit(...)` — ridge regression
- `array_lm(...)` — array-shaped responses
- `covariance(...)` / `correlation(...)` — similarity structure

If you are evaluating the package for adoption, try `many_lm(..., method = "qr")` first.
