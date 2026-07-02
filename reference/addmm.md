# Scaled matrix multiply with optional bias: alpha\*(A%\*%B) + beta\*C

Scaled matrix multiply with optional bias: alpha\*(A%\*%B) + beta\*C

## Usage

``` r
addmm(A, B, C = NULL, alpha = 1, beta = 1)
```

## Arguments

- A:

  n×p `adgeMatrix` or plain matrix.

- B:

  p×k numeric matrix.

- C:

  n×k numeric matrix or `NULL` (treated as zeros).

- alpha:

  Scalar multiplier for `A%*%B` (default 1).

- beta:

  Scalar multiplier for `C` (default 1).

## Value

`adgeMatrix` if A is resident, otherwise plain matrix.
