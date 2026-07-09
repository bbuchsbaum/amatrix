# Release readiness plan: R-universe → CRAN (2026-07-09)

Goal: `amatrix` convenient and excellent for end users on all major OSes — install core,
optionally add one backend package, get acceleration with zero or one call. Ship on
R-universe first (green dashboard), then CRAN (core package).

## Where we are (evidence-based snapshot, 2026-07-09)

**Working and green:**
- PR-gate CI: R-CMD-check green on macOS / Windows / ubuntu release+devel; pkgdown CI green; site live at bbuchsbaum.github.io/amatrix (HTTP 200).
- Local `main` == `origin/main` at `aa1a39a`; R-universe registry synced to HEAD.
- Monorepo + R-universe `subdir` builds: registry (`bbuchsbaum.r-universe.dev/packages.json`) lists core + 4 backends via subdir entries. Architecture decided and structurally working.
- First-run ergonomics implemented: silent CPU default; installed backends auto-register lazily; MLX default-on with subprocess-contained first probe; `amatrix_use_gpu()` / `amatrix_gpu_status()` with platform-aware install hints; quiet-startup option.
- The four 2026-07-02 P1 release-blocker bugs (cgk, eab, 8dy, 4jue) are fixed in git and closed in mote.

**Broken / blocking:**
1. **R-universe binary builds are red for core `amatrix` on all Linux and macOS builders** (Windows/source/wasm OK). Root-cause signal in the local check log: uncaught C++ `cl::Error: clGetDeviceIDs` aborting the process during "checking dependencies in R code" — the OpenCL probe throws when Suggests are installed on a headless builder. CRAN treats a process abort during check as fatal, so this blocks both channels.
2. **amatrix.mlx and amatrix.metal fail R-universe at the source-build step** (`Version: None`) — the builder can't even package them, so an Apple Silicon user cannot `install.packages("amatrix.mlx")`. This breaks the flagship convenience story.
3. **amatrix.arrayfire**: FAIL on Windows, ERROR on Linux/macOS builders. **amatrix.opencl**: FAIL on macOS builders + wasm, WARNING elsewhere.
4. **nightly-stress has failed 6 consecutive nights** (installed-package context: `callr` subprocess tests can't spawn child R; arrayfire-requiring test doesn't skip). Red nightly blocks our own tier-certification gate (quality-tracking.md §8.1: 7 consecutive green nights).
5. **Sub-package metadata is not releasable**: all five sub-packages are `0.0.0.9000`, self-titled "Scaffold", with placeholder author `person("Ben", "Buchsbaum", email = "bbuchsbaum@example.com")`. `amatrix.models` additionally has an invalid license DCF and 7 undocumented exports, and is absent from the registry.
6. **cran-comments.md is stale** (2026-04-12): claims 0/0/1 which no longer matches; still documents the cross-package `.Call` NOTE that Phase 2.6 was meant to kill.
7. **Tracker split-brain**: active tracker migrated to **mote** (2026-07-02); bd is a stale mirror with 8 already-fixed issues still open (incl. 4 P1s). CLAUDE.md still instructs agents to use bd.

**Real remaining engineering backlog (from mote):**
- `aul` (P1, in progress): GPU alloc sites using manual try-drop instead of `on.exit` — wrappers.R done; remaining: resident-handle.R, bind-resident.R, sinkhorn, backend-planning.R:378.
- Related GPU memory-safety bugs: 36q (double-drop), 3ka (rebind leak), 4rt (irlba upload leak), cth (dangling deferred handle), 7il (blanket tryCatch masking).
- `yux` epic (P1, blocked): benchmark-harness overhaul — gate was a no-op; 23/26 bench scripts orphaned; baseline missing 9 ops (ho0); coverage-table sync (zxn); full rehab tree is large.
- Test-coverage gaps: doj / jbs / lcn (three lists of untested exported ops).
- 0qw epic: GPU-only bug verification deferred until a real-GPU CI runner exists.
- 6b3 epic (docs/CRAN track): "when is amatrix fast?" vignette, generated compat matrix, error taxonomy (6m9), NEWS.md, as-cran checks, revdep, R-universe verify.

**Per-OS acceleration reality:**
| Platform | Backend | Status |
|---|---|---|
| macOS Apple Silicon | amatrix.mlx | The flagship path; auto-on; currently uninstallable from R-universe (blocker 2) |
| macOS (sparse) | amatrix.metal | Experimental, explicit probe |
| Linux + GPU | amatrix.opencl / amatrix.arrayfire | Opt-in via `amatrix_use_gpu()`; experimental/provisional |
| Windows | none verified | CPU-only in practice (opencl/mlx/metal are OS_type: unix; arrayfire Windows path unverified, tw4) |

## Sub-package methodology (proposed standard)

Keep the monorepo + R-universe-subdir model (decided, working). Standardize each sub-package with a release checklist that becomes CI:

1. **Metadata**: real Authors@R (Brad Buchsbaum, real email), real Title/Description (no "Scaffold"), version `0.1.0` at first release, valid LICENSE DCF, NEWS.md, URL/BugReports pointing at the monorepo.
2. **Check gate in CI**: extend R-CMD-check.yaml with a matrix job running `R CMD check` on each `backends/*` subdir on the OS(es) it targets (mlx/metal: macOS runner; opencl: ubuntu+macos; arrayfire: all three). Backends must check cleanly *without* their system library present (configure degrades to mock/unavailable; tests skip with a clear message). This is the same contract R-universe builders exercise.
3. **Graceful-degradation invariant** (the contract): loading a backend namespace on a machine without the device/library must never error, never probe hardware, and register as `available = FALSE` with a reason. Enforced by a shared conformance test sourced from core (`backend-contract.md` §validation) run inside each sub-package's test suite.
4. **Tier honesty**: README/pkgdown claims generated from `backend-certification.md` tiers; a tier is only promoted on green hardware evidence (nightly or manual certified run recorded in the ledger).
5. **CRAN policy**: core goes to CRAN with backends in Suggests + `Additional_repositories` (accepted CRAN pattern). Backends themselves stay R-universe-only for v0.1.x; arrayfire (the only portable one) is the first CRAN candidate later if demand warrants.
6. **Versioning**: sub-packages version independently; core's `amatrix_gpu_status()` reports backend package versions; contract changes bump a `Config/amatrix/contractVersion` field checked at registration.

## Decisions (resolved by Brad, 2026-07-09)

- **D1 — Tracker**: mote is authoritative; beads/`bd` retired (archival only). CLAUDE.md updated.
- **D2 — amatrix.models**: dropped from release scope. The lm/model code lives in core
  (`R/models-lm.R`); the `models/` scaffold stays out of the registry, marked not-for-release.
- **D3 — Benchmark epic scope**: minimal enforcing gate for v0.1.0 (fix gate, baseline the
  9 missing ops = ho0, coverage-table sync = zxn); the 23-script rehab tree (quc) and
  reporting polish (1x2) are post-release.
- **D4 — Windows acceleration**: IN SCOPE — Brad wants a working Windows path. Route:
  `amatrix.opencl` + CLBlast via runtime dynamic loading (see Phase W below). ArrayFire
  Windows (tw4) remains the fallback. Certification requires one manual run on a physical
  Windows machine with GPU drivers (CI runners have no GPU; CI proves build + graceful
  degradation only).
- **Author identity** (all packages): `person("Bradley", "Buchsbaum",
  email = "brad.buchsbaum@gmail.com", role = c("aut", "cre"))`.

### Phase W — Windows acceleration path (folded into Phase B/C)

1. Vendor Khronos OpenCL headers + CLBlast C-API headers (both Apache-2.0) into
   `backends/amatrix.opencl`; remove `OS_type: unix`.
2. Refactor `opencl_bridge.c` to resolve `OpenCL.dll`/`libOpenCL.so`/OpenCL.framework and
   `clblast.dll`/`libclblast` at **runtime** (LoadLibrary/dlopen + GetProcAddress), never at
   link time. Missing library, missing ICD, or zero devices ⇒ `available = FALSE` with a
   reason — never an error, never a process abort. This same refactor is the fix for the
   `clGetDeviceIDs` abort (Phase B1).
3. Add consent-gated `amatrix_install_clblast()` fetching the official CLBlast release
   binary into `tools::R_user_dir("amatrix.opencl")`; `amatrix_use_gpu()` on Windows hints
   at it when OpenCL is present but CLBlast is not.
4. CI: windows-latest job proving build-from-source + graceful degradation (no GPU).
5. Manual certification on real Windows GPU hardware; record in backend-certification.md
   before promoting the tier in README.

## Phased plan

### Phase A — Tracker + hygiene (small)
- A1. Close the 8 stale bd issues (cgk outright; eab/8dy after a local Apple Silicon re-run of the crash repros; 4jue/skj after the Phase B external verification). Resolve D1; update CLAUDE.md.
- A2. `.Rbuildignore`: add `^amatrix\.Rcheck$`, `^amatrix\.models\.Rcheck$`, `^amatrix_.*\.tar\.gz$`; fix the stray `^\.\.Rcheck$` pattern.
- A3. Refresh stale bd memories / record the mote migration.

### Phase B — Green R-universe (release gate #1)
- B1. **Kill the OpenCL probe abort.** Wrap every OpenCL entry point (starting with `clGetDeviceIDs` in `backends/amatrix.opencl/src/opencl_bridge.c`) so `cl::Error`/failures degrade to "unavailable" instead of terminating R. Add a regression test: core `R CMD check` with all backends installed on a headless machine must be abort-free. Verify the local check NOTE disappears.
- B2. **Fix mlx/metal source-build ERROR on R-universe.** Diagnose `Version: None` (likely configure failing hard on the Linux source builder); make `R CMD build` succeed anywhere (configure defers all platform checks to install time / mock bridge). Acceptance: R-universe shows source OK + macOS binaries for mlx; Apple Silicon `install.packages("amatrix.mlx")` works from a fresh library.
- B3. **arrayfire/opencl check hygiene**: tests/examples skip cleanly without the system library; clear the Windows FAIL and macOS FAIL/WARNINGs or document accepted residuals.
- B4. **Sub-package metadata pass** (methodology items 1–2): authors, titles, versions → 0.1.0, NEWS, per-subdir check CI job.
- B5. Execute D2 for amatrix.models (drop/archive or fix).
- B6. Re-verify the flagship UX end-to-end from a fresh library on Apple Silicon and on a Linux box/container: install core (+mlx where apt), one-line startup note, `amatrix_use_gpu()`, `amatrix_gpu_status()` all behave as documented. Then close 4jue/skj.

### Phase C — Green quality gates (release gate #2: "fast, bug free, excellent")
- C1. **Green the nightly**: fix `callr` child-spawn in installed-package context (or skip those blocks with reason there), skip arrayfire-requiring tests when absent. Target: 7 consecutive green nights to unlock tier certification.
- C2. **Finish the GPU memory-safety sweep**: complete `aul` remaining sites; fix 36q, 3ka, 4rt, cth; audit 7il's blanket tryCatch sites. These are correctness bugs regardless of GPU CI.
- C3. **Benchmark gate (per D3 minimal scope)**: make `tools/benchmark-regression.R` an enforcing gate, add the 9 missing baseline ops (ho0), auto-sync the coverage table (zxn); run against `tools/baseline.csv` and record results.
- C4. **Coverage gaps**: add tests for the doj/jbs/lcn op lists (delegate to test-generation agents); fix any bugs they surface.
- C5. Fresh-context adversarial review pass over the changed surfaces (per working style: no self-approval).

### Phase D — Convenience polish
- D1. README + "Get acceleration" vignette with the per-OS matrix (generated from certification tiers, 6b3 task 2); honest Windows story per decision D4.
- D2. "When is amatrix fast?" vignette (6b3 task 1).
- D3. Error-message audit / classed conditions for the top user footguns (6m9) — timebox; the full 121-site taxonomy can trail the release.
- D4. User-facing NEWS.md for 0.1.0.

### Phase E — CRAN submission (core only)
- E1. Run the `/cran-prepare` cycle: rewrite cran-comments.md, `R CMD check --as-cran` locally + win-builder + rhub, spelling/URL checks, confirm the cross-package `.Call` NOTE status (memory says killed; cran-comments still cites it — verify which is true).
- E2. Confirm Suggests-on-Additional_repositories passes incoming checks (backends must be *optional* in every example/test — already the design).
- E3. Submit; monitor; tag release; announce.

### Explicitly deferred (post-release backlog)
- Full benchmark rehab tree (quc + 9 rehab tasks, 1x2 reporting) — unless D3 goes the other way.
- 0qw GPU-CI verification epic (needs a real-GPU runner).
- arrayfire Windows path (tw4); CLBlast portable backend exploration.
- Formal backend-extension API (backends currently reach into `amatrix:::`).
- revdepcheck (no revdeps yet pre-CRAN).

## Sequencing and effort

A (hours) → B (the critical path; B1/B2 are the two unknowns — likely 1–3 days each) →
C (C1 small; C2 a day-plus; C3 per D3 scope; C4 parallelizable) → D (1–2 days, parallel with C) → E (the /cran-prepare cycle, ~a day plus CRAN latency).
B and C are independent enough to interleave; E strictly last.

## Success criteria (unchanged from deployment plan, now enforced)
1. Fresh Apple Silicon machine: `install.packages(c("amatrix","amatrix.mlx"), repos = c("https://bbuchsbaum.r-universe.dev", getOption("repos")))` → `library(amatrix)` → one-line GPU note → accelerated ops, zero env vars.
2. Any other machine: silent, correct CPU; `amatrix_use_gpu()` either enables a backend or explains why not, with install hint.
3. R-universe dashboard: core green on all builders; each backend green on its target OS and source-OK everywhere.
4. Nightly stress: 7 consecutive green; benchmark gate enforcing; zero stop-ship-rule violations.
5. CRAN: core accepted with clean or fully-explained check results.
