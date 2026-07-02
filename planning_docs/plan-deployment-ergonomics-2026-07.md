# amatrix: No-Fuss GPU Matrices for R — Deployment & Ergonomics Plan

## Context

amatrix's goal is GPU access in R with as little fuss as possible — "it just works." The package is architecturally sound (CPU-authoritative core, four optional GPU backend packages, graceful fallback), but today it is neither deployable nor low-fuss:

- **Not published:** the local repo has no git remote; github.com/bbuchsbaum/amatrix is 2.5 months stale (now public, per user); amatrix is not in the bbuchsbaum R-universe registry; no pkgdown deploy CI (docs/ committed by hand, causing the 150-file dirty tree).
- **GPU is off by default:** every backend requires an `AMATRIX_*_PROBE_GPU` env var or internal call before it will even probe the device. A user who installs `amatrix.mlx` gets silent CPU forever, with nothing telling them why.
- **The MLX `Rscript file.R` crash guard is obsolete (verified 2026-07-01):** 6 consecutive `Rscript file.R` runs on this machine (MLX 0.31.1 / mlx-c 0.6.0) exercised Metal init + GPU matmul/crossprod/svd with zero NSException aborts. The crash was documented 2026-04-09 on older MLX; upstream ml-explore/mlx#2691 closed as not-reproducible. Treat as empirically resolved, certify, and retire the guard.
- **Open correctness bugs** on the default CPU path (P0 stale-host-cache amatrix-oe8; confirmed miscompute amatrix-p24) and API gaps (kronecker demotion, missing `[`/`[[`, documented-but-nonexistent `amatrix_release_resident`).
- **Packaging hazards:** backend Suggests unresolvable (no `Additional_repositories`), tracked machine-specific `src/Makevars` with `/opt/homebrew` paths, orphan `vendor/` gitlinks, `.Rbuildignore` leaks, dangling `torch` backend references.

**Confirmed decisions:** repo public (done by user); R-universe now with backend subdir entries + `Additional_repositories`; CRAN-ready but not submitted; CPU-verifiable bug fixes in scope; GPU-only inconclusive beads and benchmark-rehab epics deferred. Defaults taken: pkgdown via CI + delete committed docs/; remove vendor gitlinks; ship `amatrix.models` if it passes check, else fix the migration doc; MLX probing becomes default-on.

**Success criterion:** new user on Apple Silicon: `install.packages(c("amatrix","amatrix.mlx"), repos=c("https://bbuchsbaum.r-universe.dev", ...))` → `library(amatrix)` → GPU active with zero env vars and a visible one-line confirmation; on any other machine: correct, silent CPU with a one-call answer to "why am I not on GPU?"

Track all work as `amatrix-*` beads (no TodoWrite). Copy this plan to `docs/plans/` (or `.planning/`) in the repo as step 0.0. Skip known-stuck beads amatrix-8w6 / amatrix-pms.

---

## Phase 0 — Repo hygiene and push (~1 session)

- **0.1 Reconnect remote.** `git remote add origin git@github.com:bbuchsbaum/amatrix.git && git fetch origin`; confirm `git log main..origin/main` is empty before pushing (local should be strictly ahead of the stale 2026-04-12 remote).
- **0.2 Remove orphan vendor gitlinks.** First `grep -rn "vendor/" backends/*/src/Makevars.in backends/*/configure src/`; if clean, `git rm --cached vendor/{CUDA-QR,MixedPrecisionBlockQR,SVDSolver,gmatrix,irlba}`; gitignore `vendor/`. *Verify:* fresh clone into scratchpad has no submodule warnings.
- **0.3 Fix `.Rbuildignore` leaks.** `^\.bug-hunt-r[0-9]+$`; add `^tmp$`, `^repro-.*\.sessionInfo\.txt$`, `^sessionInfo-probe2\.txt$`, `^tools/benchmark-.*\.R$`, `^tools/benchmark-history\.csv$` (keep `tools/benchmark-helpers.R` if tests source it — check first). Delete stray root probe artifacts and `tmp/`; gitignore them. *Verify:* `R CMD build .` then `tar -tzf` grep for `bug-hunt|benchmark|repro|sessionInfo|tmp/` → empty.
- **0.4 Commit the dirty tree coherently.** Separate commits: bug-hunt doc updates; hygiene from 0.2/0.3. Hold in-flight `R/` changes tied to oe8 for Phase 1 (land with tests). Do NOT recommit regenerated `docs/` (deleted in Phase 4); park on a temp branch/stash. Review untracked files individually; no blanket `git add .`.
- **0.5 Sync and push.** `bd doctor`, `bd sync`, `git push -u origin main`. *Verify:* `git status` → up to date with origin/main.

## Phase 1 — CPU-verifiable correctness (~1–2 sessions; 1.1–1.5 parallelizable across subagents)

Defer GPU-only inconclusive beads (-aul, -3ka, -36q, -7il, -cth, -75h) into one "GPU CI verification" epic bead.

- **1.1 amatrix-oe8 (P0) — verify, audit, pin, close.** The fix may already be in the worktree via `.amatrix_update_resident_aliases()`. Confirm invalidation exists (`grep -n host_cache_valid R/*.R`); enumerate every in-place mutation site (`am_ewise_inplace` / `am_sweep_inplace` / resident in-place surface) and check each routes through the alias update; add pinned tests in `tests/testthat/test-bughunt-residency.R` asserting materialized host data reflects mutations. *Verify:* temporarily reverting the invalidation makes the new tests fail.
- **1.2 amatrix-p24 — pairwise_sqdist_argmin CPU miscompute.** Reproduce the **verbatim original probe** first (this bead was falsely refuted once — bd memory). Fix the `rep(c_norms, each = nrow(...))` recycling in the CPU path (grep near sqdist/argmin in `R/wrappers.R`). New `tests/testthat/test-pairwise-argmin.R`: non-square X/C with distinct dims vs a plain-R double-loop reference, including a case where each-vs-times changes the answer; backend-gated cross-backend variant.
- **1.3 amatrix-lc1 — NaN-sentinel collision.** Implement the empty-vector sentinel per `.bug-hunt-r3/06-fresh-invariant.md`; tests: genuine NaN payloads survive deferred-op round-trips.
- **1.4 API gaps (one commit each, with tests):** amatrix-jnd (`kronecker`/`%x%` methods preserving class); amatrix-x6a (`[` for KronMatrix, materialize-on-subset OK for v1); amatrix-sxs (`[[` scalar extraction, Matrix-consistent); amatrix-vbh (implement + export `amatrix_release_resident()`, no-op-with-message on CPU-only).
- **1.5 amatrix-2f2 — finish copy-on-modify investigation** (resident-backed `Y <- X; Y[1,1] <- 999` must not mutate X). Fresh `Rscript -e` processes. Outcome: filed bug with repro, or documented all-clear.
- **1.6 Regression pinning + special-value conformance (timeboxed).** (a) `tests/testthat/test-regressions-hunt.R`: pin the top ~10 closed wrong-answer hunt bugs. (b) `tests/testthat/test-conformance-special-values.R`: NA/NaN/Inf through arith/matmul/crossprod/apply/sweep/bind vs base R on CPU, written as a backend loop gated by `skip_if_not_installed` + probe env vars for future GPU CI reuse.

**Gate:** `devtools::test()` clean; `Rscript tools/benchmark-regression.R` within tolerance of `tools/baseline.csv`.

## Phase 2 — Packaging and platform (~1 session)

- **2.1 Makevars hygiene (all 4 backends).** `git rm --cached backends/*/src/Makevars`; gitignore the pattern. Harden each `configure` to derive the currently hard-coded paths (opencl: `pkg-config clblast` + brew-prefix/`CLBLAST_HOME` fallbacks; SDK path via `xcrun --show-sdk-path`); `cleanup` must remove generated files. If hardening drags: detect-or-fail-clearly + polish beads — the tracked machine-specific Makevars is the actual bug. *Verify:* per backend, `rm src/Makevars && R CMD INSTALL backends/amatrix.X` regenerates and builds; `git status` stays clean after builds.
- **2.2 Platform declarations.** mlx/metal DESCRIPTIONs: `OS_type: unix` + SystemRequirements (Apple Silicon / Metal + Xcode CLT). arrayfire: `SystemRequirements: ArrayFire` + README note on the arm64-macOS CPU-runtime pin. opencl: `SystemRequirements: OpenCL, CLBlast`; Windows build path unverifiable here — file a bead, document unix-only for now. Core amatrix (portable C) is the Windows story, proven in CI (4.3).
- **2.3 Core DESCRIPTION.** Add `amatrix.metal` to Suggests; add `Additional_repositories: https://bbuchsbaum.r-universe.dev`; `amatrix.models` per decision D7 (ship if check-clean, else fix `planning_docs/package-author-migration.md:108,133`).
- **2.4 Remove dangling torch backend.** `R/policy.R:1,21` + roxygen (~122/142/273); `grep -rn torch R/ man/ tests/ vignettes/`; redocument; fix any tests asserting the policy vector.
- **2.5 Unify backend registration (`backends/*/R/zzz.R`).** Three patterns today (mlx unconditional `overwrite=TRUE`; metal/arrayfire options-gated; opencl try-wrapped). Unify: **register unconditionally, wrapped in `try`** — registration is cheap/safe; probing is what's guarded; the registry already honors per-backend disable options (`R/backend-registry.R:126`). *Verify:* both load orders per backend leave `amatrix_backend_names()` correct; double-load doesn't warn.
- **2.6 Cross-package `.Call` decoupling.** Replace `.Call(..., PACKAGE="amatrix.arrayfire"/"amatrix.mlx")` at `R/wrappers.R:2757–2822` and `R/irlba.R:154,156,184` with exported backend wrapper functions invoked via `getExportedValue()` after `requireNamespace()`; add wrappers/exports in `backends/amatrix.{arrayfire,mlx}/R/`. Kills the cross-package foreign-call check NOTE. (Reverse coupling — backends using `amatrix:::` internals — file a bead for a formal extension API; not this cycle.)

**Gate:** `rcmdcheck --as-cran` on the core tarball with no backends installed: 0E/0W; remaining NOTEs limited to Suggests availability (resolves when 4.1 lands).

## Phase 3 — First-run experience & ergonomics (~1–1.5 sessions; depends on 2.4, 2.5)

The registry already records *why* each backend is unprobed (`.amatrix_backend_health_mark`) and has `amatrix_backend_health_probe()` (10×10 canary) and `amatrix_explain()` exported — the machinery exists; it's never surfaced or made actionable.

- **3.1 Retire the MLX file-entry guard (the bugbear).**
  - (a) **Certify:** `tools/certify-mlx-file-entry.R` + loop: ~20 `Rscript file.R`-entry runs of Metal init + matmul/crossprod/svd (extend the scratchpad probe that already passed 6/6). Record results in `planning_docs/mlx-spectral-benchmark-instability.md` (resolution section); wire into `.github/workflows/nightly-stress.yaml`.
  - (b) **If clean, flip MLX probing to default-ON opt-out** (`AMATRIX_MLX_PROBE_GPU=0` / `options(amatrix.auto_probe=FALSE)` to disable): remove the file-entry guard in `backends/amatrix.mlx/R/backend.R` (`.amatrix_mlx_direct_file_entry`, ~lines 37–72), the env gate in `backends/amatrix.mlx/src/amatrix_mlx_matmul.c` (~2263–2271), and the file-entry probe policy in core `R/backend-registry.R:42–71`.
  - (c) **Belt-and-braces:** first Metal probe of a session runs in a disposable child process (same containment idea as `arrayfire_safe.cpp`); verdict cached per session; a crash kills the child, marks mlx unavailable with an informative reason, and the user's session continues on CPU. This makes default-on safe on machines we can't test.
  - Fallback if certification fails: keep the guard for `Rscript file.R` only and scope default-on to `interactive()`; still ship (c).
  - *Verify:* `Rscript file.R` with no env vars gets GPU (or a contained, explained CPU fallback); `AMATRIX_MLX_PROBE_GPU=0` disables; nightly job green.
- **3.2 `amatrix_use_gpu()` — the documented one-liner.** New `R/backend-enable.R`, exported. Walk `.amatrix_auto_fast_backend_order()` (post-torch) filtered to installed packages; per backend: load namespace, activate probe, `amatrix_backend_health_probe()`; on first healthy one, set fast-path preference and `message()` a one-line report (backend, device, float32 fast-vs-strict caveat). On failure, message each backend's recorded reason; return FALSE invisibly. Prereq refactor: `R/backend-registry.R:92` hardcodes `amatrix_mlx_enable_gpu_probe` — generalize to `spec$enable_probe_fun` and add uniform enable wrappers to metal/arrayfire/opencl (raw `Sys.setenv` sites at `backends/amatrix.{metal,arrayfire,opencl}/R/backend.R:95/125/353`). opencl/arrayfire stay opt-in by default (driver-crash/ICD-hang risk) but become one-call via this function. *Verify:* fresh session with amatrix.mlx → TRUE + mlx dispatch confirmed via `amatrix_backend_status()`; no backends installed → helpful per-platform install pointer, FALSE, no error.
- **3.3 Visible degradation.** (a) Add `.onAttach` to core `R/zzz.R`: one `packageStartupMessage` when a backend package is installed but inactive ("amatrix.mlx installed but GPU inactive — call amatrix_use_gpu()"); silent when pure-CPU or GPU already active. (b) `amatrix_gpu_status()` (new, alongside `amatrix_explain()`): per-backend table of installed/registered/probe-policy/health/reason from the registry's recorded reasons — the one-call "why am I not on GPU?". *Verify:* 4 states (none installed / installed+inactive / active / probe-failed) each produce the right message + status row; check stays clean.
- **3.4 Documented-but-missing surface sweep.** `amatrix_release_resident` lands in 1.4; `amatrix.models` per 2.3; grep every function name in README/vignettes/pkgdown reference against `NAMESPACE`. *Verify:* every advertised function exists and is exported.
- **3.5 "Getting started with GPU" docs.** New `vignettes/gpu.Rmd`: per-platform backend table (backend × OS × SystemRequirements × install command), `amatrix_use_gpu()` walkthrough, fast-vs-strict float32 semantics (1e-4 tolerance, CPU authoritative), fallback behavior, reading `amatrix_gpu_status()`; GPU chunks eval conditionally so it builds everywhere. Rewrite `README.md:118–126` install section with real r-universe commands (finalized in 4.4).

**Gate:** new-user GPU path on this machine = install two packages → `library(amatrix)` → zero calls, zero env vars → visible GPU confirmation; no-GPU path = silent correct CPU.

## Phase 4 — Deployment (~1 session)

- **4.1 R-universe registry** (via `/r-universe-publish`): add to `bbuchsbaum/bbuchsbaum.r-universe.dev` packages.json: `amatrix` (root) + subdir entries for `backends/amatrix.{mlx,opencl,arrayfire,metal}` (+ `models/amatrix.models` per D7). Mac-only backends showing red on linux/windows dashboards is expected (mitigated by 2.2 OS_type). *Verify:* dashboard green for core on all 3 OSes; `install.packages("amatrix", repos=c("https://bbuchsbaum.r-universe.dev","https://cloud.r-project.org"))` into a fresh temp lib works.
- **4.2 pkgdown deploy CI:** `.github/workflows/pkgdown.yaml` (r-lib/actions; `albersdown` via `Config/Needs/website` with the r-universe repo in extra-repositories; deploy to gh-pages). Then `git rm -r docs/`, gitignore `^docs$`, point Pages at gh-pages; audit via `/r-pkgdown-deploy`. *Verify:* Actions green; site serves including the GPU article.
- **4.3 Windows CI for core:** add `windows-latest` (release) to `R-CMD-check.yaml` matrix. *Verify:* green.
- **4.4 README finalization:** run every install command against the live registry in a fresh lib.
- **4.5 cran-comments.md refresh:** post-2.6 NOTE situation; "prepared, not submitted"; `/cran-prepare` is next cycle's entry point.

## Phase 5 — Fresh-eyes validation & close-out (~1 session)

- **5.1 Fresh-eyes walkthrough (the ergonomics gate).** Fresh-context subagent given only the public artifacts (README, pkgdown site, vignettes — not this plan), temp `.libPaths()`: (a) CPU path — install and use per README, confirm nothing GPU-related leaks in; (b) MLX path — get GPU active from docs alone, conformance ops at 1e-4, exercise `amatrix_gpu_status()`, `amatrix_release_resident()`, and degradation messaging (including a non-interactive `Rscript file.R` script). Every stumble → bead; wording/message fixes applied immediately.
- **5.2 Full gates:** `/r-cmd-check` (as-cran + lintr) clean; `devtools::test()` green; `Rscript tools/benchmark-regression.R` vs baseline within tolerance (MLX benchmarks via `Rscript -e` per existing policy until 3.1 certification retires the rule).
- **5.3 Bookkeeping + mandatory push:** close/update all touched beads; file deferred set (GPU CI epic, opencl-Windows, backend extension API, remaining unpinned hunt bugs, benchmark rehab); `bd sync`; `git pull --rebase; bd dolt push; git push; git status` → up to date.

## Sequencing & risk

- Strict order 0→1→2→3→4→5; within Phase 1, 1.1–1.5 parallelize.
- Biggest risks: 2.1 configure hardening (only verifiable on this mac for mlx/metal — fallback: detect-or-fail-clearly + beads) and 3.1's guard removal (bounded by the certification gate, the subprocess-isolated first probe, and the opt-out; can fall back to interactive()-only default-on without blocking the rest).
- Deferred explicitly: GPU-only inconclusive beads, benchmark-harness rehab epic (amatrix-yux), opencl Windows build path, formal backend extension API, CRAN submission.
