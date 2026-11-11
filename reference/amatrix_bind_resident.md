# Bind an amatrix object to resident backend storage

Upload a dense or sparse matrix to a residency-capable backend and
return the corresponding `aMatrix` object with a live resident binding.
This is primarily useful for repeated GPU work where paying the upload
cost once is preferable to relying on cold-path dispatch.

## Usage

``` r
amatrix_bind_resident(x, backend = NULL, op = NULL, y = NULL)
```

## Arguments

- x:

  An `adgeMatrix`, `adgCMatrix`, base matrix, or sparse Matrix object.

- backend:

  Backend name, `"auto"`, or `NULL`. When left `NULL` or set to
  `"auto"`, `amatrix` picks the first residency-capable accelerator
  backend that supports the requested resident operation.

- op:

  Optional operation name such as `"matmul"` used when selecting an
  automatic resident backend.

- y:

  Optional rhs object used when checking resident-op support for
  automatic backend selection.

## Value

An `adgeMatrix` or `adgCMatrix` with a live resident binding on
`backend`. When no suitable accelerator backend is available in
automatic mode, returns `x` unchanged.

## Examples

``` r
x <- adgeMatrix(matrix(rnorm(12), 3, 4))
# With only the CPU backend there is no device to bind to, so x is
# returned unchanged; with a GPU backend it becomes device-resident.
y <- amatrix_bind_resident(x)
amatrix_residency_info(y)$live
#> [1] FALSE
```
