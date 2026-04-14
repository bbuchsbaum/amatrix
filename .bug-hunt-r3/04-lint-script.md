# Round 3 Bug Hunt — Lint Script (Hunter 04)

Generated: 2026-04-14  
Scope: `R/*.R` (38 source files)  
Script runtime: 0.09 seconds

---

## (a) Script Location and CI Usage

**Script:** `tools/lint-anti-patterns.R`  
**Baseline:** `tools/lint-anti-patterns-baseline.json`

### CI Integration

Add to `.github/workflows/R-CMD-check.yaml` (or a dedicated lint workflow):

```yaml
- name: Anti-pattern lint
  run: Rscript tools/lint-anti-patterns.R
```

The script exits non-zero **only** when HIGH-severity patterns have more hits than the
baseline — i.e., regressions, not pre-existing debt. This means CI won't break on the
~90 pre-existing unclassed `stop()` calls, but will catch any new ones introduced by PRs.

### Usage Modes

```bash
# Interactive report (never exits non-zero):
Rscript tools/lint-anti-patterns.R --no-exit

# CI mode (exits 1 if HIGH net-new hits):
Rscript tools/lint-anti-patterns.R

# Refresh baseline after intentional fixes:
Rscript tools/lint-anti-patterns.R --baseline

# Write JSON for downstream tooling:
Rscript tools/lint-anti-patterns.R --json /tmp/lint-out.json
```

### Dependencies

Base R only. Uses `jsonlite` (prettier JSON) if already installed; falls back to a
hand-rolled serializer otherwise. No new package dependencies introduced.

---

## (b) Hits by Pattern — Total and Net-new vs Round 2

| Pattern | Severity | R2 Count | R3 Count | Net-new vs R2 | Notes |
|---------|----------|----------|----------|---------------|-------|
| P01 Unclassed `stop()` | HIGH | ~90 | 90 | 0 | Matches R2 exactly |
| P02 `tryCatch` → NULL | HIGH | ~80 raw | 63 (excl. 3 intentional files) | — | Intentional probes correctly excluded |
| P03 `host_cache_valid <- TRUE` | HIGH | 1 | 1 | 0 | Same single site |
| P04 `host_cache_valid <- FALSE` | HIGH | **0 (alarm)** | **1** | -1 (improvement) | R2 said 0; now 1 hit at `residency.R:137` — partial fix landed since R2 |
| P05 Double-drop pattern | HIGH | 4 confirmed bug sites | 30 grep hits | +26 (broader pattern) | Grep is wider than R2's manual audit; includes intentional drops mixed in |
| P06 Key alloc without `on.exit` | MED | ~40 (wrappers.R only) | 45 (34 wrappers.R + **11 other files**) | **+11 net-new files** | chol-factor.R, models-lm.R, resident-handle.R, bind-resident.R, backend-planning.R not in R2 scope |
| P07 `as.matrix()` dimnames drop | LOW | ~160 | 212 (excl. backend-cpu.R) | ~52 | More files covered; R2 was R/*.R but may have used narrower globs |
| P08 `stop()` without `call.=FALSE` | LOW | (subset of P01) | 87 | — | New pattern; not in R2 |
| P09 `NextMethod()` in S4 | MED | 0 | 0 | 0 | False positives (callNextMethod hits) eliminated by (?<!call) lookbehind |
| P10 `.Call()` swallow | MED | 3 confirmed (irlba.R) | **15** | **+12** | wrappers.R has 9 unexamined `.Call()` sites not in R2 |

**Total hits (R3): 544 across 10 patterns**

---

## (c) Net-new Bug Candidates

All candidates are **inferred** (lint smell confirmed by grep; underlying bug inferred,
not execution-verified). Confidence labels: `HIGH-inferred` = strong structural evidence;
`MED-inferred` = plausible but needs code-path analysis.

### Candidate 1 — P06: Key allocs without `on.exit` in non-wrappers.R files (HIGH-inferred)

**Files and lines:**
- `R/chol-factor.R:180, 303, 343`
- `R/models-lm.R:579`
- `R/resident-handle.R:81, 136, 181, 258, 319`
- `R/bind-resident.R:88`
- `R/backend-planning.R:378`

**Why inferred:** Round 2 and existing beads issue `amatrix-aul` cover only the ~40
wrappers.R sites. These 11 sites in other files have identical structural risk: key
allocated, then if an error fires between allocation and the manual `try(resident_drop())`
at the end, the key leaks permanently. `resident-handle.R` is particularly high-risk
because it handles in-place GPU operations where error paths are more likely.

**Confidence: HIGH-inferred.** The pattern is identical to the confirmed wrappers.R leak.
Not covered by any open beads issue.

---

### Candidate 2 — P10: `.Call()` bridges in `wrappers.R` not examined in R2 (MED-inferred)

**Lines:** `wrappers.R:1954, 2490, 2494, 2505, 2508, 2528, 2532, 2548, 2555`

**Why inferred:** Round 2 only audited `irlba.R` Lanczos bridges (filed as `amatrix-833`
implicitly via the sweep). The 9 `.Call()` sites in wrappers.R dispatch to ArrayFire/MLX
C bridges (`am_af_dist_sq_bridge`, `am_af_kernel_bridge`, `amatrix_mlx_tcrossprod_bridge`,
`am_sparse_segment_sum_c`). These are called inside dispatch arms that may be inside
`tryCatch(error = function(e) NULL)` blocks — if so, partial C-state mutation before throw
is silently swallowed. Needs code-path trace to confirm tryCatch wrapping.

**Confidence: MED-inferred.** Grep confirms the `.Call()` sites exist; tryCatch wrapping
not verified by grep alone.

---

### Candidate 3 — P04: `host_cache_valid <- FALSE` at `residency.R:137` is new since R2 (LOW-inferred, positive signal)

**Line:** `R/residency.R:137` — `cs$host_cache_valid <- FALSE`

**Why interesting:** Round 2 reported 0 sites for this pattern, declaring the missing
invalidation a confirmed HIGH bug. The script now finds 1 site. This means either:
(a) a fix landed between R2 and now — partial progress on `amatrix-dev`, or  
(b) R2's grep missed this line (different field name: R2 grepped `host_cache_valid <- FALSE`,
   current code uses `cs$host_cache_valid <- FALSE` via a local alias).

Either way, the fix is incomplete: the `broadcast_ewise_resident_inplace_key` path in
`resident-handle.R:155-168` still likely lacks invalidation. The lint confirms the
pattern exists in one place but the coverage may still be insufficient.

**Confidence: LOW-inferred** (positive finding, not a new bug per se — flags incomplete
fix coverage for `amatrix-dev`).

---

## Summary

| | Count |
|--|--|
| Total lint hits | 544 |
| Patterns with net-new vs R2 | 2 (P06 +11 files, P10 +12 sites) |
| Net-new bug candidates | 3 |
| Confirmed HIGH-inferred new candidates | 1 (P06 non-wrappers files) |
| Confirmed MED-inferred new candidates | 1 (P10 wrappers .Call sites) |
| Positive findings (partial fix detected) | 1 (P04 partial invalidation) |

The script and baseline are checked-in to `tools/` and ready for CI.
