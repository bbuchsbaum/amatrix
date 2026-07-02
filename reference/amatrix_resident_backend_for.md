# Choose a residency-capable accelerator backend for a hot path

Returns the first available non-CPU backend that can keep `x` resident
for the requested operation. This is intended for package authors who
want repeated work to stay on the fastest available accelerator without
hardcoding backend names such as `"metal"` or `"mlx"`.

## Usage

``` r
amatrix_resident_backend_for(x, op = NULL, y = NULL)
```

## Arguments

- x:

  An `aMatrix`.

- op:

  Optional operation name such as `"matmul"`.

- y:

  Optional rhs object used when checking resident-op support.

## Value

A backend name, or `NULL` when no residency-capable accelerator is
available for the requested operation.
