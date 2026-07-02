# Clear the amatrix backend fallback log

Resets the fallback log to empty. Typically called at the start of a
test block to isolate the assertion that the log is empty after a clean
run.

## Usage

``` r
amatrix_fallback_log_reset()
```

## Value

Invisibly, `NULL`.

## See also

[`amatrix_fallback_log`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_fallback_log.md)
