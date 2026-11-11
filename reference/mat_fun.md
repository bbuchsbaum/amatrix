# Matrix functions via symmetric eigendecomposition

Apply an elementwise function to the eigenvalues of a symmetric positive
definite matrix and reconstruct the result: \\f(X) = Q \\
\mathrm{diag}(f(\lambda)) \\ Q^T\\.

## Usage

``` r
mat_sqrt(X)

mat_pow(X, p)

mat_log(X)
```

## Arguments

- X:

  Symmetric positive definite numeric matrix or `adgeMatrix` of shape
  `[p, p]`.

- p:

  Numeric scalar exponent (used by `mat_pow` only).

## Value

An `adgeMatrix` of shape `[p, p]`: the matrix function applied to `X`.

## Examples

``` r
S <- crossprod(matrix(rnorm(16), 4)) + diag(4)
mat_sqrt(S)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>             [,1]         [,2]       [,3]         [,4]
#> [1,]  2.29537050  0.400923037 -0.5592246 -0.021803150
#> [2,]  0.40092304  1.174122387 -0.2950168 -0.004055295
#> [3,] -0.55922458 -0.295016753  3.1322639  0.385857389
#> [4,] -0.02180315 -0.004055295  0.3858574  1.132102003
mat_log(S)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>             [,1]        [,2]       [,3]       [,4]
#> [1,]  1.57127632  0.46007852 -0.4034181 0.02223645
#> [2,]  0.46007852  0.22574621 -0.2560563 0.03101484
#> [3,] -0.40341811 -0.25605632  2.2030453 0.40402871
#> [4,]  0.02223645  0.03101484  0.4040287 0.19011614
mat_pow(S, -1)
#> An amatrix dense matrix [cpu|policy=auto|precision=strict]
#> 4 x 4 Matrix of class "adgeMatrix"
#>             [,1]        [,2]        [,3]        [,4]
#> [1,]  0.25602453 -0.19738959  0.05424844 -0.02881547
#> [2,] -0.19738959  0.86630872  0.07187720 -0.04374878
#> [3,]  0.05424844  0.07187720  0.14812557 -0.15400348
#> [4,] -0.02881547 -0.04374878 -0.15400348  0.86739827
```
