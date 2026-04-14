# Benchmark Harness Compliance Contract

Every rehabilitated benchmark script must satisfy this checklist. This contract exists to prevent the bugs documented in `.omc/research/bench-audit/bugs.md` and enforce patterns that enable CI integration and reliable regression detection.

---

## (a) Entry Point Signature

Every benchmark script must expose exactly one top-level function with signature:

```r
benchmark_<name>_main(opts = list())
```

where `opts` is a parsed option list (typically from `commandArgs(trailingOnly = TRUE)` via a parser like `parse_args()`).

The script's top-level source must guard autorun:

```r
if (!isTRUE(as.logical(Sys.getenv("AMATRIX_BENCHMARK_NO_AUTORUN", "FALSE")))) {
  benchmark_<name>_main()
}
```

This allows scripts to be sourced in non-autorun mode for testing or integration into harnesses without side effects. Reference: `benchmark-svd-backends.R` (lines 624-631).

---

## (b) Honored Environment Variables

Scripts must respect the following env vars. They may not be required, but if set, they must be honored:

- `AMATRIX_BENCHMARK_NO_AUTORUN` — If `"1"` or `"TRUE"`, suppress top-level function execution (allows sourcing without running).
- `AMATRIX_BENCHMARK_SIZES` — Comma-separated list of size labels (e.g., `"small,medium,large"`). Script must select case sizes from this list.
- `AMATRIX_BENCHMARK_BACKENDS` — Comma-separated backend names (e.g., `"cpu,mlx,opencl"`). Script must skip backends not in this list.
- `AMATRIX_BENCHMARK_REPS` — Number of repetitions (integer). Default: 3. Script uses this to set `reps` parameter.
- `AMATRIX_BENCHMARK_REPO_ROOT` — Absolute path to repo root. If set, use it for resolving relative paths; fall back to normalizePath of parent directories.

Document any script-specific env vars (e.g., `AMATRIX_SVD_FACTOR_CASES`, `AMATRIX_MLX_NATIVE_SPECTRAL`) in the script header.

---

## (c) Backend Skip Pattern

Each backend must be gated on availability. Use `requireNamespace("<backend>", quietly = TRUE)` or the shared helper `ensure_optional_backend_namespace()` (from `benchmark-helpers.R`). If a backend is missing:

1. Emit exactly one line: `skipped: <backend> not installed`
2. Continue to the next backend — never `stop()`.

Reference pattern from `benchmark-helpers.R:64-124`:

```r
if (!requireNamespace("amatrix.mlx", quietly = TRUE)) {
  message("skipped: mlx not installed")
  next  # or continue
}
```

This ensures that missing optional backends do not block the entire script. CPU backend is always required.

---

## (d) Timer Discipline

Timing measurements must follow strict discipline to avoid upload/GC noise:

1. **Prime host→device upload outside the timer** using the shared `prime_backend()` helper (once amatrix-fjo lands). This call happens once per backend, before the rep loop, and is not timed.

2. **Explicit garbage collection before each rep**: Use `benchmark_time_ms()` from `benchmark-helpers.R` (lines 354-372) which calls `gc()` before each rep. Alternatively, implement your own timer that calls `gc()` before each iteration.

3. **Cold vs. warm variants** must use distinct primitives:
   - **Cold**: Clear residency / caches before each rep (e.g., call `drop_svd_factor_cache()` before each `benchmark_time_ms()` rep).
   - **Warm**: Execute one unreplicated warmup call before the timer, then run the timer with the warmed state. Reference: `benchmark-regression.R:238-242` (exact SVD warmup).

4. Default `bench::mark` uses `gc = FALSE`. If you use `bench::mark`, either add `gc = TRUE` or accept the lower timing precision.

Reference: `benchmark-svd-backends.R:120-133` and `benchmark-regression.R:354-372`.

---

## (e) Output Format (raw-results.csv schema)

Results must be written as a CSV with these columns in order:

```
op, size, backend, variant, median_ms, mean_ms, sd_ms, p05_ms, p95_ms, n_reps, rel_err, nnz_metadata
```

Column semantics:

- `op`: Operation name (string, e.g., `"matmul"`, `"svd"`, `"chol"`).
- `size`: Size label (string, e.g., `"small"`, `"medium"`, `"large"`). Never use NxP format like `"1024x128"` as the join key (it is not stable).
- `backend`: Backend name (string, e.g., `"cpu"`, `"mlx"`, `"opencl"`, `"arrayfire"`).
- `variant`: Variant label (string, e.g., `"cold"`, `"warm"`, `"resident"`). Always include variant.
- `median_ms`, `mean_ms`, `sd_ms`, `p05_ms`, `p95_ms`: Timing percentiles in milliseconds (numeric). Compute from reps.
- `n_reps`: Number of repetitions (integer).
- `rel_err`: Relative error vs. CPU reference (numeric). Use `NA` if not applicable (e.g., CPU backend itself, or op has no accuracy check).
- `nnz_metadata`: Non-zero count for sparse matrices (integer or NA). This is **metadata only** — never use it as a join key in regression comparisons.

Reference: `benchmark-svd-backends.R:139-155` (new_row pattern).

---

## (f) Registration in benchmark-regression.R Dispatch Table

To make the script runnable via the main harness, register it in `benchmark-regression.R`'s dispatch table at the top level (around line 1065+):

```r
# Near line 1065 in benchmark-regression.R, add to the dense_ops list (or sparse_ops):
dense_ops <- c("matmul", "crossprod", ..., "your_op_name")
```

Then source the script at initialization (e.g., after line 1000) and call its entry point:

```r
# In initialize_regression_benchmark_context():
source(file.path(repo_root, "tools", "benchmark-your-name.R"), local = FALSE)
```

Update `benchmark-regression.R` line 1120+ (the `run_master` dispatch switch) to invoke your script:

```r
if (op %in% dense_ops) {
  case <- benchmark_dense_case(op, size_label = size_label, ...)
  ...
  if (identical(op, "your_op_name")) {
    results <- benchmark_your_name_main(opts = list(...))
  }
}
```

This ensures the harness can find and invoke your script via `Rscript tools/benchmark-regression.R --op=your_op_name`.

---

## (g) CI Runtime Budget

Every benchmark script must complete in under **90 seconds** on the default `macos-latest` GitHub Actions runner with default sizes. Measure this locally:

```bash
time Rscript tools/benchmark-your-name.R
```

If your script exceeds 90 seconds, either optimize the kernel calls or add a larger size override:

```r
# Allow opt-in to larger size via env var
if (Sys.getenv("AMATRIX_BENCHMARK_SIZES", "") == "large") {
  # Use larger matrices
} else {
  # Use default (small/medium)
}
```

This prevents CI timeouts and allows developers to opt into longer runs for detailed analysis.

---

## (h) Accuracy Check

Every operation that runs on a non-CPU backend must call an accuracy validation function (once amatrix-teo lands) to ensure GPU results match CPU reference within tolerance:

```r
assert_backend_accuracy(ref_cpu, gpu_result, op, tol = tol_for_op)
```

Helper signature (amatrix-teo):

```r
assert_backend_accuracy <- function(ref_result, gpu_result, op_name, tol = 1e-4) {
  # Compute relative error: rel_err <- norm(gpu - ref) / norm(ref)
  # If rel_err > tol, throw with op_name and tol in message
}
```

Tolerance defaults by op family (amatrix-teo specifies these; for now document in your script header):

- Factorizations (QR, SVD, Cholesky, LU): `tol = 1e-4` (GPU fast precision)
- Linear solves: `tol = 1e-3`
- Eigenvalues: `tol = 1e-3`
- Matrix products: `tol = 1e-6`

Reference: `benchmark-svd-backends.R:135-137` (relative_sv_error pattern).

---

## Forward Compatibility Notes

This contract will be extended as in-progress helpers land:

- **amatrix-fjo** implements `prime_backend()` — use it once merged (section d).
- **amatrix-teo** implements `assert_backend_accuracy()` — use it once merged (section h).
- **amatrix-5dx** implements per-rep sd/CI calculations — cross-reference in accuracy checks (section h).

When these tickets close, update the cross-references in this document.

---

## Checklist for Script Authors

- [ ] Entry point function named `benchmark_<name>_main(opts = list())`
- [ ] Top-level autorun guard checks `AMATRIX_BENCHMARK_NO_AUTORUN`
- [ ] All backends gated on `requireNamespace()` or `ensure_optional_backend_namespace()`
- [ ] Missing backends emit skip message and continue (no `stop()`)
- [ ] Timer uses `benchmark_time_ms()` or `bench::mark` with explicit `gc = TRUE`
- [ ] Cold/warm variants use distinct warmup primitives
- [ ] Output CSV has all 12 required columns in order
- [ ] Script registered in `benchmark-regression.R` dispatch table
- [ ] Script completes in <90 seconds on default sizes
- [ ] Accuracy check implemented for non-CPU backends (or documented as pending amatrix-teo)
- [ ] Script-specific env vars documented in header
- [ ] Tested with `AMATRIX_BENCHMARK_NO_AUTORUN=1` to verify sourcing without autorun

