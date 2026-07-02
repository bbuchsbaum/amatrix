# GPU-accelerated pairwise distance matrix

Computes the pairwise distance matrix between rows of `X` and `Y`. The
dominant cost (row inner-products via am_tcrossprod) is dispatched to
the active GPU backend (ArrayFire or MLX); norm computation and final
transforms run on CPU where they are O(mp + np) — negligible versus the
O(mnp) GEMM.

## Usage

``` r
dist_matrix(
  X,
  Y = NULL,
  method = c("euclidean", "sqeuclidean", "cosine"),
  tile_size = NULL
)
```

## Arguments

- X:

  Numeric matrix or `adgeMatrix`, shape \[m, p\].

- Y:

  Numeric matrix or `adgeMatrix`, shape \[n, p\], or `NULL` to compute
  pairwise distances within `X` (returns \[m, m\] matrix).

- method:

  One of `"euclidean"` (default), `"sqeuclidean"`, or `"cosine"`.

- tile_size:

  Integer row-block size for tiled computation, or `NULL` (default) to
  auto-tile when `nrow(X) > 50000` (self-distance only). Set explicitly
  to process any size in row-blocks; useful when GPU memory is limited.
  Not supported for `method = "cosine"`.

## Value

Numeric matrix \[m, n\] of pairwise distances.

## See also

[`kernel_matrix`](https://bbuchsbaum.github.io/amatrix/reference/kernel_matrix.md)

## Examples

``` r
X <- matrix(rnorm(30), nrow = 6)
D <- dist_matrix(X)
dim(D)
#> [1] 6 6
```
