# Hunter 05 — parallel workers (callr / mclapply)

## (a) Drift check

- `bd list --status=open | rg -i "parallel|fork|callr|mclapply|future|worker"`
  found no existing parallel-workers issue.
- Existing nearby open issues:
  `amatrix-1ha` (bare `library(amatrix)` generic gap),
  `amatrix-1i1` (deferred `adgeMatrix` unrecoverable after `saveRDS/readRDS`).
- Initial worker probe was confounded by `amatrix-1ha`: `rowSums(X)` failed
  before any worker transfer logic was exercised.

## (b) Scenario

- Started with `tmp/bug-hunt/parallel/probe.R`:
  plain dense `adgeMatrix`, `callr::r()`, `parallel::mclapply()`,
  `saveRDS/readRDS`, and worker `%*%`.
- Because the bare-attach generic gap dominated the result, I pivoted to a
  tighter resident-object surface relevant to worker transfer:
  `resident_handle` materialization in a fresh subprocess after
  `library(amatrix)`.
- Minimal fresh-process regression test:
  `tests/testthat/test-regression-resident-handle-as-matrix.R`.

## (c) Findings

- Reproduced a new bug on the resident-handle surface:
  documented generic `as.matrix(h)` does not dispatch to
  `as.matrix.resident_handle()` in a fresh subprocess after
  `library(amatrix)`.
- Direct namespace call works:
  `amatrix:::as.matrix.resident_handle(h)` returns the expected matrix.
- Generic calls fail:
  `as.matrix(h)` and `base::as.matrix(h)` return the error
  `cannot coerce type 'environment' to vector of type 'any'`.
- This is filing-worthy and now captured as a regression repro. It is not
  actually a worker-specific registry corruption bug; the parallel-worker hunt
  surfaced it because resident handles are the realistic object that users would
  try to move or materialize across process boundaries.

## (d) Proposed bd create

- Title: `resident_handle cannot be materialized with generic as.matrix() after library(amatrix)`
- Priority: `P2`
- Type: `bug`
- Why P2:
  documented public API for `resident_handle` materialization is broken on the
  ordinary installed-package path, but the failure is loud and a workaround
  exists through direct/internal conversion paths.

## (e) Limitations

- I did not reach a true fork/worker residency-corruption repro; the hunt
  pivoted once a narrower, real resident-handle bug reproduced first.
- `callr` was not reliable in this environment for the test harness, so the
  permanent repro uses a base `Rscript` subprocess instead.
