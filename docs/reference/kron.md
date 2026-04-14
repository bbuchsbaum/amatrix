# Eager Kronecker product

Computes `A ⊗ B` and returns the result as an `adgeMatrix`. Accepts
plain matrices or any `aMatrix` subclass. For a lazy variant that avoids
forming the full product see
[`kron_matrix`](https://bbuchsbaum.github.io/amatrix/reference/kron_matrix.md).

## Usage

``` r
kron(A, B)
```

## Arguments

- A, B:

  Matrices or `aMatrix` objects.

## Value

An `adgeMatrix` of dimension `(nrow(A)*nrow(B)) x (ncol(A)*ncol(B))`.

## See also

[`kron_matrix`](https://bbuchsbaum.github.io/amatrix/reference/kron_matrix.md)
