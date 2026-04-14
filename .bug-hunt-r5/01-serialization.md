# Hunter 01 — Serialization round-trip scenario

## (a) Drift check

- Installed version: 0.1.0
- DESCRIPTION last touched: 2026-04-12 15:12:43 -0400 (commit before HEAD)
- HEAD commit: 2695f6e 2026-04-14 18:30:43 -0400 "Add amatrix_benchmark_report() and hardware-aware calibration cache"
- Discrepancy: installed namespace is from a prior build; HEAD is 2 commits ahead of DESCRIPTION touch.
  Proceeded anyway — discrepancy does not affect serialization surface.

## (b) Scenario

Probed saveRDS/readRDS round-trips across fresh R sessions (separate Rscript -e calls):

1. **Resident handle dangling**: normal adgeMatrix (CPU backend, no GPU residency) — save/reload
2. **@env / finalizer_env slot**: does the environment survive across serialization?
3. **Deferred NaN sentinel**: new_adgeMatrix_deferred (host_deferred=TRUE, @x=NaN) — save/reload
4. **Post-reload arithmetic without library(amatrix)**: fresh R with only library(Matrix)
5. **Empty namespace**: fresh R with no libraries

## (c) Findings

### H1 — Resident handle dangling (CPU backend)
No bug. CPU backend stores data in @x (host). After reload: finalizer_env contents = {cache_state, object_id},
host_deferred=FALSE. as.matrix() returns correct values. The finalizer_env IS serialized correctly
(environment serializes its bindings). The residency registry (.amatrix_state$residency) is per-session
and empty in the new process, but CPU-backend objects fall through to reading @x, which is intact.

### H2 — @env / finalizer_env across save/load
No bug for normal objects. The environment slot serializes and deserializes correctly including its
bindings (object_id, cache_state child env). The finalizer is NOT preserved (expected — R does not
serialize finalizers), but that only means the GC cleanup hook is lost, which is acceptable for a
reloaded object.

### H3 — Deferred NaN sentinel (CONFIRMED BUG — two related issues)

**Bug A (P1, amatrix-90k):** `as.matrix()` on a reloaded deferred adgeMatrix silently returns NaN.

Root cause: `as.matrix` is not an S4 generic. `setMethod("as.matrix","adgeMatrix")` registers a
method but S3 dispatch from `base::as.matrix()` bypasses it, falling through to the Matrix-package
inherited coercion that reads `@x` directly (the NaN sentinel vector). Result: caller gets a
matrix of NaN with no error, no warning.

Evidence:
```r
# Session 1
library(amatrix)
X <- amatrix:::new_adgeMatrix_deferred(dim=c(2L,3L), preferred_backend="cpu")
saveRDS(X, "/tmp/deferred.rds")

# Session 2
library(amatrix)
X2 <- readRDS("/tmp/deferred.rds")
as.matrix(X2)        # => NaN NaN NaN  (SILENT — no error)
# but:
amatrix_materialize_host(X2)  # => error: deferred adgeMatrix lost its GPU resident data
show(X2)             # => error after printing header
# and:
f <- selectMethod("as.matrix", "adgeMatrix"); f(X2)  # => error (correct)
```

**Bug B (P2, amatrix-1i1):** Deferred adgeMatrix is completely unrecoverable after round-trip.
No serialization hook materializes host data before saving, and no initialize/load hook detects the
dead state on reload. The object arrives in the new session with host_deferred=TRUE, host_x=NULL,
empty residency registry — a permanent dead end.

### H4 — Post-reload without library(amatrix)
No bug. R auto-loads the amatrix namespace via `requireNamespace` / lazy loading when the S4 class
is encountered. Output shows "Loading required namespace: amatrix". Methods dispatch correctly,
arithmetic works, as.matrix works.

### H5 — Empty namespace (no libraries at all)
Same behavior as H4: amatrix namespace is auto-loaded on first access. No crash, no silent failure.

## (d) Proposed bd create

Both filed:
- **amatrix-90k** (P1): silent NaN from as.matrix() on deferred adgeMatrix after reload
- **amatrix-1i1** (P2): no serialization hook for deferred adgeMatrix (unrecoverable after round-trip)

## (e) Limitations

- Only CPU backend tested (no MLX/ArrayFire available in this environment). Resident-handle
  dangling for GPU backends (H1 GPU path) was not runtime-verified — only the CPU path was probed.
- The `as.matrix` S3/S4 dispatch discrepancy (Bug A) likely also affects `as.numeric`, `as.vector`,
  `as.array` on deferred objects — not probed separately.
- model-cache-backed objects (amChol, amQR) were not tested for serialization; out of scope given
  CPU-only environment.
- H3 Bug A is the more dangerous because it produces incorrect results silently; H3 Bug B produces
  errors which at least alert the user.
