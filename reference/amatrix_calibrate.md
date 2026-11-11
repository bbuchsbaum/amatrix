# Calibrate GPU dispatch thresholds for this machine

Runs micro-benchmarks for each (op, backend, size) combination and
derives the minimum matrix element count at which each GPU backend
reliably outperforms CPU. Results are stored in the current session and
optionally persisted to disk for reuse in future sessions.

## Usage

``` r
amatrix_calibrate(
  backend = NULL,
  ops = c("gemm", "gemv", "crossprod", "rowSums", "colSums", "qr", "chol", "solve",
    "svd"),
  sizes = list(c(64L, 32L), c(128L, 64L), c(256L, 128L), c(512L, 256L), c(1024L, 512L)),
  sparse_densities = c(0.01, 0.05, 0.2),
  n_reps = 10L,
  margin = 0.1,
  persist = FALSE,
  quiet = FALSE
)
```

## Arguments

- backend:

  Character vector of backend names to benchmark. Defaults to all
  registered non-CPU backends that report `available = TRUE`.

- ops:

  Character vector of operations to benchmark. Supported values:
  `"matmul"` (alias for `"gemm"`), `"gemm"`, `"gemv"`, `"spmv"`,
  `"spmm"`, `"crossprod"`, `"rowSums"`, `"colSums"`, `"qr"`, `"chol"`,
  `"solve"`, `"svd"`.

- sizes:

  List of integer vectors of length 2 giving `c(nrow, ncol)` test matrix
  dimensions.

- sparse_densities:

  Numeric vector of target fill densities used when benchmarking
  `"spmv"` and `"spmm"`.

- n_reps:

  Positive integer. Number of timed repetitions per benchmark cell,
  after 2 warm-up repetitions.

- margin:

  Non-negative numeric less than 1. Fraction by which GPU median time
  must beat CPU to count as a GPU win (default `0.10` means GPU must be
  at least 10% faster).

- persist:

  Logical. If `TRUE`, save the calibration to the user cache directory
  (`tools::R_user_dir("amatrix", "cache")`) so future sessions load it
  automatically. The default `FALSE` keeps it in the current session
  only; nothing is written to disk unless you ask for it.

- quiet:

  Logical. Suppress progress messages.

## Value

Invisibly, a list with elements `version`, `calibrated_at` (POSIXct),
`thresholds` (nested list keyed by backend then op), and `results`
(data.frame of all benchmark measurements).

## See also

[`amatrix_calibration_info`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_calibration_info.md),
[`amatrix_backend_plan`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_backend_plan.md)

## Examples

``` r
# Benchmarks every registered non-CPU backend. With only the CPU backend
# registered there is nothing to measure and this returns immediately.
cal <- amatrix_calibrate(quiet = TRUE)
names(cal)
#> [1] "version"       "calibrated_at" "sys_hash"      "thresholds"   
#> [5] "results"      
```
