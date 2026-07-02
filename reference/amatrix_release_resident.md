# Release GPU-resident data held by an amatrix object

Frees any device-resident buffer associated with `x` and drops its
residency-registry binding, leaving the host copy as the authoritative
storage. This gives long-lived GPU pipelines explicit control over
device memory instead of waiting for garbage collection to reclaim
resident handles.

## Usage

``` r
amatrix_release_resident(x)
```

## Arguments

- x:

  An
  [`aMatrix`](https://bbuchsbaum.github.io/amatrix/reference/aMatrix-class.md)
  object (for example an
  [`adgeMatrix`](https://bbuchsbaum.github.io/amatrix/reference/adgeMatrix-class.md)
  or
  [`adgCMatrix`](https://bbuchsbaum.github.io/amatrix/reference/adgCMatrix-class.md)).
  Non-amatrix inputs are ignored.

## Value

Invisibly, `TRUE` if a resident binding was released, and `FALSE`
otherwise (including the CPU-only no-op case).

## Details

The object remains fully usable afterwards: its data is served from the
host copy and is re-uploaded to the device on the next GPU operation if
needed. On CPU-only sessions, or for any object that currently holds no
device buffer, this is a safe no-op.

## Examples

``` r
A <- adgeMatrix(matrix(1:6, 2, 3))
# On a CPU-only session there is no device buffer, so this is a no-op:
released <- amatrix_release_resident(A)
released
#> [1] FALSE
```
