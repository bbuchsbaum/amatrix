# Dense general matrix with backend-dispatch metadata

`adgeMatrix` extends both `aMatrix` and `Matrix::dgeMatrix`, adding
backend-dispatch slots to a column-major dense double-precision matrix.
All arithmetic generics dispatch through the amatrix backend system
rather than directly to BLAS.

## See also

[`adgeMatrix`](https://bbuchsbaum.github.io/amatrix/reference/adgeMatrix.md)
for the user-facing constructor,
[`adgCMatrix`](https://bbuchsbaum.github.io/amatrix/reference/adgCMatrix.md)
for the sparse counterpart
