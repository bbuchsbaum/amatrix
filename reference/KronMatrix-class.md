# Lazy Kronecker product of two matrices

`KronMatrix` stores the two factor matrices `A` (m x n) and `B` (p x q)
without forming the full (mp x nq) Kronecker product. Matrix-vector and
matrix-matrix products are evaluated using the vec-permutation identity
`(A x B) vec(X) = vec(B X t(A))`, keeping memory use at `O(mn + pq)`
rather than `O(mnpq)`.

## Slots

- `A`:

  Numeric matrix; the left factor of the Kronecker product.

- `B`:

  Numeric matrix; the right factor of the Kronecker product.

## See also

[`kron_matrix`](https://bbuchsbaum.github.io/amatrix/reference/kron_matrix.md)
