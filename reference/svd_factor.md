# Compute a truncated SVD of an aMatrix

Computes a rank-`k` truncated singular value decomposition of `X`,
dispatching to the most suitable algorithm and backend based on matrix
size, requested rank, and available hardware.

## Usage

``` r
svd_factor(
  X,
  k = min(dim(X)),
  method = c("auto", "exact", "rsvd", "subspace"),
  n_oversamples = .amatrix_svd_factor_default_oversamples(k),
  n_iter = 2L
)
```

## Arguments

- X:

  An `aMatrix` (typically `adgeMatrix` or `adgCMatrix`).

- k:

  Positive integer; the number of singular values and vectors to
  compute. Defaults to `min(dim(X))`.

- method:

  One of `"auto"` (default), `"exact"`, `"rsvd"`, or `"subspace"`.
  `"auto"` selects the algorithm based on matrix dimensions and the
  available backend.

- n_oversamples:

  Non-negative integer; extra random vectors used to stabilize
  randomized algorithms. Ignored for `method = "exact"`.

- n_iter:

  Non-negative integer; number of power iterations for the subspace
  method. Ignored for `method = "exact"`.

## Value

An
[`amSVD`](https://bbuchsbaum.github.io/amatrix/reference/amSVD-class.md)
object containing slots `u`, `d`, `v`, and metadata.

## Examples

``` r
m <- matrix(rnorm(30), nrow = 6)
A <- adgeMatrix(m)
fac <- svd_factor(A, k = 3L)
fac
#> amSVD [6x5 -> rank 3 | strict | exact/exact_svd@cpu | source: 20260710043614.171005-4c29c517:am:118]
#>   d[1:min(3,k)]: 3.6559, 2.912, 2.0014
```
