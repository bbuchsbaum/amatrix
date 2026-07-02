# Scalar extraction from an adgeMatrix

`[[` extracts a single element of an
[`adgeMatrix`](https://bbuchsbaum.github.io/amatrix/reference/adgeMatrix-class.md),
matching base matrix semantics: `A[[i]]` selects the `i`-th element in
column-major order, and `A[[i, j]]` selects the element in row `i`,
column `j`. Without this method `A[[1]]` fails with “this S4 class is
not subsettable” (as it does for a bare `dgeMatrix`).

## Usage

``` r
# S4 method for class 'adgeMatrix'
x[[i, j, ...]]
```

## Arguments

- x:

  An
  [`adgeMatrix`](https://bbuchsbaum.github.io/amatrix/reference/adgeMatrix-class.md).

- i, j:

  Element subscripts, following base matrix `[[` semantics.

- ...:

  Unused.

## Value

A length-one numeric value.

## Details

The host copy is materialized before extraction (transparently
downloading from the device for GPU-resident objects), so the returned
value always reflects the current contents.

## Examples

``` r
A <- adgeMatrix(matrix(1:6, 2, 3))
A[[1]]
#> [1] 1
A[[2, 3]]
#> [1] 6
```
