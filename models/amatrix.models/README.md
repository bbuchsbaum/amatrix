# amatrix.models

`amatrix.models` is the model-facing surface built on top of the `amatrix` substrate.

Start here:

```r
library(amatrix)
library(amatrix.models)

X <- adgeMatrix(design, preferred_backend = "mlx", precision = "fast")
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
- resident MLX QR path that now pays off when the number of right-hand sides grows

Other exported surfaces:

- `lm_fit(...)`
- `wls_fit(...)`
- `ridge_fit(...)`
- `array_lm(...)`
- `covariance(...)`
- `correlation(...)`

If you are evaluating the package for adoption, try `many_lm(..., method = "qr")` first.
