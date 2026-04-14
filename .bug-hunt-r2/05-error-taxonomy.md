# amatrix — Error Taxonomy Audit (Round 2, bug hunt)

Companion to **amatrix-cng** (R1). Round 1 flagged that exported functions should raise
typed conditions rather than bare `stop("...")`. This document enumerates every site,
proposes a stable class hierarchy, and maps each site to the class it should adopt.

## 1. Counts

| Metric | Count |
|---|---|
| Total `stop(` calls in `R/` | **121** |
| Typed (`class=` / `rlang::abort(class=)`) | **0** |
| Unclassed bare `stop()` | **121** |
| `warning(` calls | **4** (all unclassed) |
| `rlang::abort(` calls | **0** |
| Files touched | 19 |
| Exported functions | 120 (103 `export()` + 1 `exportMethods()` + 16 `S3method()`) |
| Exported-reachable unclassed sites | **≈115 of 121** |
| Purely-internal (dev-guard) sites | ≈6 (e.g. `R/backend-planning.R:418 stop(e)` rethrow; some `.amatrix_*` helpers never in an exported call chain) |

Every single `stop()` in the package is unclassed. This is a 100% violation rate
against the R1 rule. The few "internal-only" sites are still bugs because they may
propagate through exported wrappers.

### stop() density per file

| File | Count |
|---|---:|
| R/wrappers.R | 15 |
| R/sinkhorn.R | 15 |
| R/backend-registry.R | 12 |
| R/svd-factor.R | 12 |
| R/chol-factor.R | 10 |
| R/qr.R | 10 |
| R/models-lm.R | 9 |
| R/resident-handle.R | 9 |
| R/irlba.R | 8 |
| R/qr-downdate.R | 5 |
| R/bind-resident.R | 4 |
| R/policy.R | 3 |
| R/constructors.R | 2 |
| R/residency.R | 2 |
| R/backend-calibration.R | 1 |
| R/backend-cpu.R | 1 |
| R/backend-planning.R | 1 |
| R/kronmatrix.R | 1 |
| R/product-plan.R | 1 |

## 2. Proposed class hierarchy

All conditions inherit from a common root `amatrix_error` so users can catch "any amatrix
error" with a single handler, and from `error`/`condition` as required by R.

```
amatrix_error                                 (parent; all typed amatrix errors)
├── amatrix_error_invalid_input               (bad user argument)
│   ├── amatrix_error_invalid_type            (wrong R type / class)
│   ├── amatrix_error_invalid_value           (out-of-range, non-finite, negative, etc.)
│   └── amatrix_error_dim_mismatch            (length/nrow/ncol mismatch)
├── amatrix_error_invalid_factor              (wrong factor object passed: amChol/amLU/amSVD/amQR)
├── amatrix_error_backend                     (backend-plane issues)
│   ├── amatrix_error_backend_unavailable     (not registered / probe failed)
│   ├── amatrix_error_backend_unsupported_op  (backend lacks capability/feature)
│   ├── amatrix_error_backend_contract        (registration / contract violations)
│   └── amatrix_error_backend_failed          (backend kernel returned bad result)
├── amatrix_error_residency                   (resident-handle lifecycle)
├── amatrix_error_unsupported_op              (amatrix itself doesn't support the path)
├── amatrix_error_numerical                   (rank-deficient, singular, non-convergence)
├── amatrix_error_missing_dependency          (optional package like irlba not installed)
└── amatrix_error_internal                    (invariant broken; should be unreachable)
```

Suggested constructors in `R/errors.R` (to be created by fix pass):

```r
amatrix_abort <- function(message, class, ..., call = rlang::caller_env()) {
  rlang::abort(message, class = c(class, "amatrix_error"), ..., call = call)
}
```

### Catch examples for users

```r
tryCatch(am_qr(X),
  amatrix_error_dim_mismatch = function(e) { ... },
  amatrix_error_backend_unavailable = function(e) { ... },
  amatrix_error = function(e) { ... }   # fallthrough
)
```

## 3. Recurring failure modes (frequency table)

Sorted by how many unique sites share the same semantic failure. Mode is the proposed class.

| Failure mode | Sites | Class |
|---|---:|---|
| Argument type / class check (`x must be a ...`) | ~32 | `amatrix_error_invalid_type` |
| Argument value check (positive int, non-neg, finite) | ~18 | `amatrix_error_invalid_value` |
| Wrong factor class (amChol/amLU/amSVD/amQR/aMatrix) | 13 | `amatrix_error_invalid_factor` |
| Dimension / length mismatch (`length(w) != nrow(X)`, `nrow(Y) != ncol(u)`, etc.) | 12 | `amatrix_error_dim_mismatch` |
| Backend unavailable / not registered | 5 | `amatrix_error_backend_unavailable` |
| Backend capability / feature missing (residency, broadcast kernels, sparse) | 9 | `amatrix_error_backend_unsupported_op` |
| Backend registration contract violation | 10 | `amatrix_error_backend_contract` |
| Backend kernel failed / returned unusable result | 6 | `amatrix_error_backend_failed` |
| Residency lifecycle (inactive / lost / cannot materialize) | 5 | `amatrix_error_residency` |
| Unsupported amatrix path (LINPACK, sparse sinkhorn, unsupported payload) | 8 | `amatrix_error_unsupported_op` |
| Numerical (singular, rank-deficient, non-convergence) | 6 | `amatrix_error_numerical` |
| Missing R package dependency (`irlba`) | 1 | `amatrix_error_missing_dependency` |
| Internal invariant / rethrow | 2 | `amatrix_error_internal` |

## 4. Site-by-site mapping (file:line → recommended class)

All sites are **⚠️ UNCLASSED** in the current code — the column lists the class to assign.

### R/wrappers.R
| Line | Current call | Recommended class |
|---:|---|---|
| 78 | `stop(sprintf("unsupported bind kind '%s'", kind))` | `amatrix_error_invalid_value` |
| 1098 | `stop("x must be a matrix-like object")` | `amatrix_error_invalid_type` |
| 1187 | `stop("LINPACK is not supported")` | `amatrix_error_unsupported_op` |
| 1271 | `stop("length(w) must equal nrow(X)")` | `amatrix_error_dim_mismatch` |
| 1333 | `stop("length(w) must equal nrow(X)")` | `amatrix_error_dim_mismatch` |
| 1399 | `stop("length(w) must equal nrow(X)")` | `amatrix_error_dim_mismatch` |
| 1404 | `stop("nrow(y) must equal nrow(X)")` | `amatrix_error_dim_mismatch` |
| 1482 | `stop("length(d) must equal nrow(X)")` | `amatrix_error_dim_mismatch` |
| 1494 | `stop("length(d) must equal ncol(X)")` | `amatrix_error_dim_mismatch` |
| 1556 | `stop("lambda must be a scalar or length ncol(X)")` | `amatrix_error_dim_mismatch` |
| 1639 | `stop("n must be supplied when using solve_fn")` | `amatrix_error_invalid_value` |
| 2794 | `stop(arg_name, " must be a list of matrices or a 3-D array ...")` | `amatrix_error_invalid_type` |
| 2830 | `stop("Ls must be a list of amChol objects ...")` | `amatrix_error_invalid_factor` |
| 2839 | `stop("B must be a list or a 3-D array ...")` | `amatrix_error_invalid_type` |
| 2843 | `stop("Ls and B must have the same batch size")` | `amatrix_error_dim_mismatch` |

### R/sinkhorn.R
| Line | Current call | Recommended class |
|---:|---|---|
| 49 | `stop("return_info must be TRUE or FALSE")` | `amatrix_error_invalid_value` |
| 83 | `stop(sprintf("%s must be a single positive integer", name))` | `amatrix_error_invalid_value` |
| 87 | `stop(sprintf("%s must be a single positive integer", name))` | `amatrix_error_invalid_value` |
| 94 | `stop(sprintf("%s must be a single finite numeric value", name))` | `amatrix_error_invalid_value` |
| 98 | `stop(sprintf("%s must be non-negative", name))` | `amatrix_error_invalid_value` |
| 101 | `stop(sprintf("%s must be positive", name))` | `amatrix_error_invalid_value` |
| 108 | `stop("sinkhorn() currently requires a dense matrix or adgeMatrix")` | `amatrix_error_unsupported_op` |
| 126 | `stop("A must be a dense numeric matrix or adgeMatrix")` | `amatrix_error_invalid_type` |
| 158 | `stop("A must be numeric")` | `amatrix_error_invalid_type` |
| 161 | `stop("A must be two-dimensional")` | `amatrix_error_invalid_type` |
| 164 | `stop("sinkhorn() currently requires a square matrix")` | `amatrix_error_dim_mismatch` |
| 167 | `stop("A must contain only finite non-missing values")` | `amatrix_error_invalid_value` |
| 170 | `stop("A must be elementwise non-negative")` | `amatrix_error_invalid_value` |
| 173 | `stop("A must have strictly positive row sums and column sums")` | `amatrix_error_invalid_value` |
| 265 | `stop("resident sinkhorn path requires resident row/col reductions ...")` | `amatrix_error_backend_unsupported_op` |

### R/backend-registry.R
| Line | Current call | Recommended class |
|---:|---|---|
| 122 | `stop("backend must be a named list")` | `amatrix_error_backend_contract` |
| 140 | `stop("backend is missing required fields: ...")` | `amatrix_error_backend_contract` |
| 144 | `stop("backend$capabilities must be a function")` | `amatrix_error_backend_contract` |
| 147 | `stop("backend$features must be a function")` | `amatrix_error_backend_contract` |
| 150 | `stop("backend$precision_modes must be a function")` | `amatrix_error_backend_contract` |
| 157 | `stop("backend$capabilities() must return a character vector")` | `amatrix_error_backend_contract` |
| 160 | `stop("backend$features() must return a character vector")` | `amatrix_error_backend_contract` |
| 163 | `stop("backend$precision_modes() must return a character vector")` | `amatrix_error_backend_contract` |
| 166 | `stop(sprintf(...))` (invalid precision modes) | `amatrix_error_backend_contract` |
| 174 | `stop(sprintf("backend '%s' is already registered", name))` | `amatrix_error_backend_contract` |
| 187 | `stop(sprintf("backend '%s' is not registered", name))` | `amatrix_error_backend_unavailable` |
| 326 | `stop(sprintf("backend '%s' is not registered", name))` | `amatrix_error_backend_unavailable` |

### R/svd-factor.R
| Line | Current call | Recommended class |
|---:|---|---|
| 379 | `stop("subspace SVD did not discover a usable range space")` | `amatrix_error_numerical` |
| 411 | `stop("projected core is numerically rank-deficient")` | `amatrix_error_numerical` |
| 566 | `stop("X must be an aMatrix")` | `amatrix_error_invalid_type` |
| 570 | `stop("k must be a positive integer")` | `amatrix_error_invalid_value` |
| 574 | `stop(sprintf("k (%d) cannot exceed min(dim(X)) (%d)", k, max_k))` | `amatrix_error_invalid_value` |
| 580 | `stop("n_oversamples must be a non-negative integer")` | `amatrix_error_invalid_value` |
| 583 | `stop("n_iter must be a non-negative integer")` | `amatrix_error_invalid_value` |
| 682 | `stop("factor must be an amSVD object")` | `amatrix_error_invalid_factor` |
| 691 | `stop(sprintf("Y has %d rows but factor@u has %d rows", ...))` | `amatrix_error_dim_mismatch` |
| 700 | `stop(sprintf("Y has %d rows but factor@u has %d rows", ...))` | `amatrix_error_dim_mismatch` |
| 734 | `stop("factor must be an amSVD object")` | `amatrix_error_invalid_factor` |
| 739 | `stop(sprintf("Z has %d rows but factor@k is %d", ...))` | `amatrix_error_dim_mismatch` |

### R/chol-factor.R
| Line | Current call | Recommended class |
|---:|---|---|
| 127 | `stop("X must be an adgeMatrix (symmetric positive definite)")` | `amatrix_error_invalid_type` |
| 398 | `stop(arg_name, " must be a list of RHS objects or a 3-D array ...")` | `amatrix_error_invalid_type` |
| 435 | `stop("factor must be an amChol object")` | `amatrix_error_invalid_factor` |
| 499 | `stop("factor must be an amChol object")` | `amatrix_error_invalid_factor` |
| 516 | `stop("all RHS batches must have nrow equal to the factor dimension")` | `amatrix_error_dim_mismatch` |
| 563 | `stop("factor must be an amChol object")` | `amatrix_error_invalid_factor` |
| 589 | `stop("factor must be an amChol object")` | `amatrix_error_invalid_factor` |
| 701 | `stop("factor must be an amChol object")` | `amatrix_error_invalid_factor` |
| 788 | `stop("A must be a square matrix")` | `amatrix_error_dim_mismatch` |
| 820 | `stop("factor must be an amLU object")` | `amatrix_error_invalid_factor` |

### R/qr.R
| Line | Current call | Recommended class |
|---:|---|---|
| 152 | `stop("unsupported QR payload")` | `amatrix_error_unsupported_op` |
| 196 | `stop("unsupported QR payload")` | `amatrix_error_unsupported_op` |
| 337 | `stop("source reconstruction is only supported for explicit_qr payloads")` | `amatrix_error_unsupported_op` |
| 344 | `stop("qr must inherit from amQR")` | `amatrix_error_invalid_factor` |
| 361 | `stop("compact QR factor is unavailable")` | `amatrix_error_unsupported_op` |
| 414 | `stop("qr must inherit from amQR")` | `amatrix_error_invalid_factor` |
| 540 | `stop("explicit QR resident q could not be materialized")` | `amatrix_error_backend_failed` |
| 547 | `stop("explicit QR payload does not contain q")` | `amatrix_error_internal` |
| 722 | `stop("only square matrices can be inverted")` | `amatrix_error_dim_mismatch` |
| 730 | `stop("singular matrix 'a' in solve")` | `amatrix_error_numerical` |

### R/models-lm.R
| Line | Current call | Recommended class |
|---:|---|---|
| 39 | `stop("x must be a dense matrix-like object")` | `amatrix_error_invalid_type` |
| 612 | `stop("weights must be a numeric vector of non-missing non-negative values")` | `amatrix_error_invalid_value` |
| 616 | `stop("weights must have length equal to nrow(X)")` | `amatrix_error_dim_mismatch` |
| 620 | `stop("weights must have positive total weight")` | `amatrix_error_invalid_value` |
| 656 | `stop("block_size must be NULL or a single positive integer")` | `amatrix_error_invalid_value` |
| 661 | `stop("block_size must be NULL or a single positive integer")` | `amatrix_error_invalid_value` |
| 750 | `stop("effective denominator must be positive")` | `amatrix_error_numerical` |
| 769 | `stop("effective denominator must be positive")` | `amatrix_error_numerical` |
| 974 | `stop("lambda must be a single non-negative numeric value")` | `amatrix_error_invalid_value` |

### R/resident-handle.R
| Line | Current call | Recommended class |
|---:|---|---|
| 62 | `stop("x must be an adgeMatrix or matrix")` | `amatrix_error_invalid_type` |
| 67 | `stop("backend '...' does not support residency")` | `amatrix_error_backend_unsupported_op` |
| 109 | `stop("resident_handle is no longer active")` | `amatrix_error_residency` |
| 170 | `stop("backend failed to apply resident vector sweep")` | `amatrix_error_backend_failed` |
| 174 | `stop("backend does not support broadcast_ewise_resident_key")` | `amatrix_error_backend_unsupported_op` |
| 199 | `stop("backend failed to apply resident vector sweep")` | `amatrix_error_backend_failed` |
| 251 | `stop("backend does not support broadcast_ewise_resident")` | `amatrix_error_backend_unsupported_op` |
| 290 | `stop("backend does not support ewise_resident")` | `amatrix_error_backend_unsupported_op` |
| 294 | `stop("rhs must be a scalar or resident_handle")` | `amatrix_error_invalid_type` |

### R/irlba.R
| Line | Current call | Recommended class |
|---:|---|---|
| 76 | `stop("Package 'irlba' is required: install.packages('irlba')")` | `amatrix_error_missing_dependency` |
| 147 | `stop("nv must be >= 1 after clamping to work size")` | `amatrix_error_invalid_value` |
| 334 | `stop("basis_cols must be a non-negative integer")` | `amatrix_error_invalid_value` |
| 774 | `stop("nv must be a positive integer")` | `amatrix_error_invalid_value` |
| 777 | `stop("nu must be a positive integer")` | `amatrix_error_invalid_value` |
| 786 | `stop("block_size must be a positive integer")` | `amatrix_error_invalid_value` |
| 791 | `stop("block_size must not exceed matrix dimensions")` | `amatrix_error_invalid_value` |
| 795 | `stop("n_steps must be a positive integer")` | `amatrix_error_invalid_value` |

### R/qr-downdate.R
| Line | Current call | Recommended class |
|---:|---|---|
| 43 | `stop("row_idx must be a single positive integer")` | `amatrix_error_invalid_value` |
| 48 | `stop(sprintf("row_idx must be between 1 and %d", n_rows))` | `amatrix_error_invalid_value` |
| 57 | `stop("qr_downdate.amQR requires the original matrix X ...")` | `amatrix_error_invalid_value` |
| 67 | `stop("X must be a matrix-like object with rows")` | `amatrix_error_invalid_type` |
| 77 | `stop("qr_downdate.default requires an amQR factor ...")` | `amatrix_error_invalid_factor` |

### R/bind-resident.R
| Line | Current call | Recommended class |
|---:|---|---|
| 42 | `stop("x must be an adgeMatrix, adgCMatrix, matrix, or sparse Matrix")` | `amatrix_error_invalid_type` |
| 62 | `stop(sprintf("backend '%s' is not available", backend_name))` | `amatrix_error_backend_unavailable` |
| 65 | `stop(sprintf("backend '%s' does not support residency", backend_name))` | `amatrix_error_backend_unsupported_op` |
| 92 | `stop(sprintf("backend '%s' does not support sparse residency", backend_name))` | `amatrix_error_backend_unsupported_op` |

### R/policy.R
| Line | Current call | Recommended class |
|---:|---|---|
| 155 | `stop(sprintf("policy must be one of: %s", ...))` | `amatrix_error_invalid_value` |
| 218 | `stop(sprintf(...))` (legacy policy spec) | `amatrix_error_invalid_value` |
| 249 | `stop(sprintf("precision must be one of: %s", ...))` | `amatrix_error_invalid_value` |

### R/constructors.R
| Line | Current call | Recommended class |
|---:|---|---|
| 67 | `stop("x must be a base matrix or dgeMatrix")` | `amatrix_error_invalid_type` |
| 88 | `stop("x must be a base matrix or dgCMatrix")` | `amatrix_error_invalid_type` |

### R/residency.R
| Line | Current call | Recommended class |
|---:|---|---|
| 400 | `stop("deferred adgeMatrix lost its GPU resident data")` | `amatrix_error_residency` |
| 436 | `stop("resident backend returned an unsupported dense materialization type")` | `amatrix_error_backend_failed` |

### R/kronmatrix.R
| Line | Current call | Recommended class |
|---:|---|---|
| 171 | `stop("solve() requires square A and B factor matrices")` | `amatrix_error_dim_mismatch` |

### R/backend-calibration.R
| Line | Current call | Recommended class |
|---:|---|---|
| 379 | `stop(sprintf("unsupported calibration ops: %s", ...))` | `amatrix_error_invalid_value` |

### R/backend-cpu.R
| Line | Current call | Recommended class |
|---:|---|---|
| 125 | `stop("LINPACK is not supported")` | `amatrix_error_unsupported_op` |

### R/backend-planning.R
| Line | Current call | Recommended class |
|---:|---|---|
| 418 | `stop(e)` (rethrow of CPU error) | keep rethrow but wrap as `amatrix_error_internal` with parent in `body` |

### R/product-plan.R
| Line | Current call | Recommended class |
|---:|---|---|
| 190 | `stop("x must be matrix-like")` | `amatrix_error_invalid_type` |

## 5. warning() sites (unclassed)

| File:line | Call | Recommended class |
|---|---|---|
| R/irlba.R:58 | `warning(...)` fallback message | `amatrix_warning_fallback` |
| R/irlba.R:219 | `warning("irlba_native: did not converge after ... restarts")` | `amatrix_warning_nonconvergence` |
| R/wrappers.R:1569 | `warning("matrix has non-positive eigenvalues; result may be complex or NaN", ...)` | `amatrix_warning_numerical` |
| R/policy.R:71 | `warning(...)` | `amatrix_warning_policy` |

Mirror the error tree with an `amatrix_warning` parent.

## 6. Notes on reachability

- `R/backend-registry.R` sites are all reachable through exported `amatrix_register_backend`
  / `amatrix_backend_*` helpers.
- `R/backend-calibration.R:379` is reached via exported `amatrix_calibrate`.
- `R/backend-planning.R:418 stop(e)` sits inside the planner fallback engine used by almost
  every exported dispatch path — this rethrow must preserve / wrap the typed class so
  handlers still match downstream.
- `R/policy.R:218` lives inside `.amatrix_parse_legacy_spec`, called from exported
  `amatrix_set_default_policy`.
- `R/resident-handle.R` sites are reachable through exported `resident_handle`,
  `am_ewise_inplace`, `am_sweep_inplace`, and `amatrix_bind_resident`.

## 7. Recommended fix strategy

1. Add `R/errors.R` with `amatrix_abort()` and `amatrix_warn()` helpers that call
   `rlang::abort()` / `rlang::warn()` with a class vector ending in `"amatrix_error"` /
   `"amatrix_warning"`.
2. Add `R/zzz-errors-doc.R` (or a roxygen block in `R/errors.R`) documenting the class
   tree so users can `?amatrix_errors`.
3. Mechanically rewrite each site per the table above. Favor the most specific class;
   callers can always catch the parent.
4. Replace the `R/backend-planning.R:418 stop(e)` rethrow with
   `amatrix_abort(conditionMessage(e), class = "amatrix_error_internal", parent = e)`.
5. Add `tests/testthat/test-error-classes.R` asserting each exported function raises the
   expected class for its representative failure modes (dim mismatch, invalid factor,
   backend unavailable, etc.).
6. Register the new classes in `NAMESPACE` exports only if helper functions become
   exported; the class names themselves are conventional strings, no export needed.
