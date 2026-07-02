# GPU-accelerated pairwise kernel matrix

Computes the pairwise kernel matrix between rows of `X` and `Y`. The
expensive am_tcrossprod is GPU-dispatched; element-wise transforms (exp,
sqrt, pow) run on CPU.

## Usage

``` r
kernel_matrix(
  X,
  Y = NULL,
  kernel = c("linear", "rbf", "polynomial", "cosine", "laplacian"),
  sigma = 1,
  degree = 2L,
  coef = 0,
  preferred_backend = NULL,
  zero_diag = FALSE
)
```

## Arguments

- X:

  Numeric matrix or `adgeMatrix`, shape \[m, p\].

- Y:

  Numeric matrix or `adgeMatrix`, shape \[n, p\], or `NULL`.

- kernel:

  Kernel type string (see Details).

- sigma:

  Bandwidth for `"rbf"` and `"laplacian"`.

- degree:

  Polynomial degree for `"polynomial"`.

- coef:

  Constant term for `"polynomial"`: \\(\mathrm{coef} + x^\top
  y)^{\mathrm{degree}}\\.

- preferred_backend:

  Optional backend name to override the default dispatch (e.g., `"mlx"`,
  `"opencl"`).

- zero_diag:

  When `TRUE` and `Y` is `NULL`, set the diagonal of the kernel matrix
  to zero.

## Value

Numeric matrix \[m, n\] of kernel values.

## Details

Kernels:

- linear:

  \\k(x,y) = x^\top y\\

- rbf:

  \\k(x,y) = \exp(-\\x - y\\^2 / (2 \sigma^2))\\

- polynomial:

  \\k(x,y) = (\mathrm{coef} + x^\top y)^{\mathrm{degree}}\\

- cosine:

  \\k(x,y) = (x^\top y) / (\\x\\ \\y\\)\\

- laplacian:

  \\k(x,y) = \exp(-\\x - y\\ / \sigma)\\

## See also

[`dist_matrix`](https://bbuchsbaum.github.io/amatrix/reference/dist_matrix.md)
