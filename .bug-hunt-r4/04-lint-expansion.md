# Round-4 Lint Expansion — Hunter 04

## A. HEAD & Drift Check

- HEAD: `eaf8c43`
- `git describe --always`: `eaf8c43`
- Branch: `main`, no drift from expected HEAD

## B. New Rules Added (P11–P16)

| Rule | Regex / Mechanism | Severity |
|------|-------------------|----------|
| P11_alloc_no_on_exit | `\.amatrix_next_resident_key\(|resident_next_key\(` | HIGH |
| P12_double_drop_same_key | `try\(backend\$resident_drop\((?:scaled_key\|out_key\|y_key\|x_scaled_key\|y_scaled_key)` | HIGH |
| P13_trycatch_blanket_fallback | `error\s*=\s*function\s*\(e\)\s*(?:NULL\|character\(\)\|...)` | MED |
| P14_nan_deferred_sentinel | `rep\s*\(\s*NaN\s*,` | HIGH |
| P15_backend_arg_unused | body-scan: `function(...backend=...)` with no `\bbackend\b` in next 40 lines | MED |
| P16_unclassed_stop_density | `stop\("` — density summary only (top-5 files) | LOW |

## C. Script Runtime

**0.14 seconds** (38 R source files, 16 rules total)

## D. Per-Rule Hit Count

| Pattern | Severity | Hits | Baseline | Net-new |
|---------|----------|------|----------|---------|
| P01_unclassed_stop | HIGH | 89 | 90 | -1 |
| P02_trycatch_swallow_null | HIGH | 65 | 63 | +2 |
| P03_host_cache_valid_set_true | HIGH | 1 | 1 | 0 |
| P04_host_cache_valid_never_false | HIGH | 1 | 1 | 0 |
| P05_double_drop_pattern | HIGH | 30 | 30 | 0 |
| P06_resident_key_alloc_no_on_exit | MED | 45 | 45 | 0 |
| P07_as_matrix_dimnames_drop | LOW | 220 | 212 | +8 |
| P08_missing_call_false_stop | LOW | 86 | 87 | -1 |
| P09_s4_nextmethod_fallthrough | MED | 0 | 0 | 0 |
| P10_trycatch_c_bridge_swallow | MED | 15 | 15 | 0 |
| P11_alloc_no_on_exit | HIGH | 45 | — (new) | — |
| P12_double_drop_same_key | HIGH | 23 | — (new) | — |
| P13_trycatch_blanket_fallback | MED | 46 | — (new) | — |
| P14_nan_deferred_sentinel | HIGH | 1 | — (new) | — |
| P15_backend_arg_unused | MED | 0 | — (new) | — |
| P16_unclassed_stop_density | LOW | 89 | — (new) | — |

**P16 top-5 files by unclassed stop() density:**

| File | stop() count |
|------|-------------|
| `R/wrappers.R` | 17 |
| `R/qr.R` | 10 |
| `R/sinkhorn.R` | 10 |
| `R/models-lm.R` | 9 |
| `R/resident-handle.R` | 9 |

## E. Net-New Hits & Existing Issue Cross-Reference

### P11_alloc_no_on_exit — 45 hits
**Existing issue:** `amatrix-aul` ("GPU key leaks: ~40 alloc sites in wrappers.R use manual try-drop instead of on.exit"). Also `amatrix-4rt` (irlba GPU upload not cleaned up). P11 is a duplicate of P06 with a clearer description; confirms the same 45 sites. No new filing needed.

### P12_double_drop_same_key — 23 hits
**Existing issues:** `amatrix-8kj` (wrappers.R:1283/1345), `amatrix-74d` (wrappers.R:1415-1452), `amatrix-4q9` (chol-factor.R:175-191). The 23 hits capture `try(backend$resident_drop(scaled_key|out_key|y_key|x_scaled_key|y_scaled_key)` — the unconditional-drop half of the double-drop pair. Net-new finding: `backend-planning.R:395` and `models-lm.R:620` are NOT covered by existing issues. **2 net-new sites.**

### P13_trycatch_blanket_fallback — 46 hits
**Existing issue:** `amatrix-e4w` covers the SVD-specific swallow. No existing issue covers the 46-site pattern broadly. Most sites in `backend-registry.R`, `chol-factor.R`, `qr.R`. **Net-new pattern** — not fully covered. Overlaps with P02 (NULL-only) but extends to `character()`, `integer()`, `list()` fallbacks. ~10 non-NULL sites are genuinely new.

### P14_nan_deferred_sentinel — 1 hit
**Existing issue:** `amatrix-lc1` ("NaN-as-deferred-sentinel collides with user NaN data") and `amatrix-ax8` ("amChol stores NaN-sentinel factor_obj"). Hit confirms `constructors.R:199` as the creation site. Fully covered. No new filing.

### P15_backend_arg_unused — 0 hits
Body-scan finds no functions where `backend=` is declared and the body (within 40 lines) never mentions `backend`. The 40-line heuristic may miss long functions; not a confirmed finding.

### P16_unclassed_stop_density — density only
**Existing issue:** `amatrix-6m9` ("Error taxonomy: 121 unclassed stop() calls — full audit and class hierarchy"). Top files `wrappers.R` (17), `qr.R` (10), `sinkhorn.R` (10) are the highest-priority targets. No new filing needed — informs prioritization of amatrix-6m9.

### P02 net-new: +2 hits
Two new NULL-swallow sites added since baseline (likely `backend-planning.R:439,455` — `tryCatch(signalCondition(...), error = function(ee) NULL)`). These are intentional signal-swallows but should be documented. Covered by `amatrix-e4w` scope.

## F. 1-Shot Rscript Promotion Snippets

### P12: backend-planning.R:395 double-drop
```r
# Rscript -e 'source("R/backend-planning.R"); cat(readLines("R/backend-planning.R")[390:400], sep="\n")'
Rscript -e 'lines <- readLines("R/backend-planning.R"); cat(paste(390:400, lines[390:400], sep="\t"), sep="\n")'
```

### P12: models-lm.R:620 double-drop
```r
Rscript -e 'lines <- readLines("R/models-lm.R"); cat(paste(615:625, lines[615:625], sep="\t"), sep="\n")'
```

### P13: non-NULL blanket fallbacks (net-new sites)
```r
Rscript -e '
  files <- list.files("R", pattern="\\.R$", full.names=TRUE)
  rx <- "error\\s*=\\s*function\\s*\\(e\\)\\s*(?:character\\(\\)|integer\\(\\)|numeric\\(\\)|list\\(\\)|logical\\(\\))"
  for (f in files) {
    lines <- readLines(f, warn=FALSE)
    hits <- grep(rx, lines, perl=TRUE)
    if (length(hits)) cat(basename(f), hits, "\n")
  }
'
```

### P14: constructors.R:199 NaN sentinel creation
```r
Rscript -e 'lines <- readLines("R/constructors.R"); cat(paste(195:205, lines[195:205], sep="\t"), sep="\n")'
```

## G. CI Integration Proposal

Add the following step to `.github/workflows/check.yml` after `R CMD check`:

```yaml
- name: Anti-pattern lint (R4 rules)
  run: Rscript tools/lint-anti-patterns.R
  # Exits non-zero only on HIGH-severity net-new hits above baseline.
  # Update baseline via: Rscript tools/lint-anti-patterns.R --baseline
```

To gate on MED-severity regressions as well, add `--strict` support to the script
(flag `high_new_hits + med_new_hits > 0`). Current threshold is HIGH-only to avoid
blocking on known tech-debt patterns (P06/P11 have 45 open sites).
