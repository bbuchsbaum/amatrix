# Nearest-centroid assignment via fused squared-distance computation

Computes \\D\[i,k\] = \\x_i\\^2 - 2 x_i^\top c_k + \\c_k\\^2\\ and
returns \\\arg\min_k D\[i,k\]\\ for each row \\i\\, 1-indexed. GPU path
avoids host round-trips by chaining resident operations.

## Usage

``` r
pairwise_sqdist_argmin(X, Ct, x_norms = NULL, c_norms = NULL)
```

## Arguments

- X:

  n×p `adgeMatrix` or plain matrix (query points).

- Ct:

  p×K numeric matrix (centroids, transposed — columns are centroids).

- x_norms:

  Optional n-vector of precomputed \\\\x_i\\^2\\. Computed if `NULL`.

- c_norms:

  Optional K-vector of precomputed \\\\c_k\\^2\\. Computed if `NULL`.

## Value

Integer vector of length n, 1-indexed nearest centroid per row.
