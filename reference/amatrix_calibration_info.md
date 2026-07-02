# Retrieve the current calibration state

Returns the calibration data stored in the current session. If no
calibration has been run yet, the function attempts to load a previously
persisted calibration from the user cache directory. Returns `NULL` when
no calibration is available.

## Usage

``` r
amatrix_calibration_info()
```

## Value

A list as returned by
[`amatrix_calibrate`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_calibrate.md),
or `NULL` if no calibration data exists for this session.

## See also

[`amatrix_calibrate`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_calibrate.md)

## Examples

``` r
cal <- amatrix_calibration_info()
is.null(cal) # TRUE when no calibration has been run
#> [1] TRUE
```
