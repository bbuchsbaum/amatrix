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

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(20), 5, 4)
K <- kernel_matrix(X, kernel = "rbf", sigma = 1.5)
dim(K)
#> [1] 5 5
# Cross-kernel between two sets of rows
kernel_matrix(X, X[1:2, ], kernel = "linear")
#>             [,1]        [,2]
#> [1,]  3.35311408  0.07512096
#> [2,]  0.07512096  0.42355183
#> [3,] -1.06387913 -0.05104410
#> [4,] -4.85682215 -0.30306633
#> [5,]  1.71810337  0.34058802
```
