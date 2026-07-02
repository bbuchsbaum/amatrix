# Report amatrix benchmark status across ops and backends

Reads the machine-local baseline in `tools/baseline.csv` (if present)
and the cached calibration in the user cache directory, and returns a
structured data.frame surfacing per-op cold vs warm timings and the
currently-calibrated dispatch thresholds.

## Usage

``` r
amatrix_benchmark_report(baseline_path = file.path("tools", "baseline.csv"))
```

## Arguments

- baseline_path:

  Path to the baseline CSV. Defaults to
  `file.path(getwd(), "tools", "baseline.csv")` so it works when called
  from the package source tree. Pass `NULL` to skip baseline reading
  entirely and return only calibration data.

## Value

A list with two elements:

- baseline:

  data.frame with columns `op`, `size`, `backend`, `cold_ms`, `warm_ms`,
  `warm_vs_cold_ratio`, `speedup_vs_cpu`. Rows with missing cold OR warm
  data use `NA` for the missing variant. Empty when the baseline file is
  absent.

- calibration:

  data.frame with columns `backend`, `op`, `threshold_elements`,
  `gpu_wins`. Rows come from the cached calibration; empty when no
  calibration is available.

## Details

This is the user-facing honesty surface for Track 4's speed contract:
users can see (a) which backends are calibrated on their machine, (b)
cold-start vs warm-run ratios per op, and (c) where the dispatcher will
currently route.

## See also

[`amatrix_calibrate`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_calibrate.md),
[`amatrix_calibration_info`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_calibration_info.md)

## Examples

``` r
if (FALSE) { # \dontrun{
rep <- amatrix_benchmark_report()
head(rep$baseline)
head(rep$calibration)
} # }
```
