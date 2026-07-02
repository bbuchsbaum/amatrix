# Hunter 01 — copy-on-modify semantics
## (a) Drift check
- `bd search "copy-on-modify"` returned no matching open issue.
- `bd list --status=open` showed no existing bug explicitly covering alias mutation of resident-backed `adgeMatrix`/`adgCMatrix`.

## (b) Scenario
- Fresh-process probes via `Rscript -e` with `pkgload::load_all('.')`.
- Registered a mock residency-capable backend (`make_recording_backend`) to force live resident bindings.
- Probed three cases:
- `X <- amatrix_bind_resident(adgeMatrix(...)); Y <- X; Y[1,1] <- 999`
- `Z <- A %*% B` on resident-capable backend, then `W <- Z; W[1,1] <- 777`
- `S <- amatrix_bind_resident(adgCMatrix(...)); T <- S; T[1,1] <- 111`

## (c) Findings
- No bug found.
- In every case, the original object retained its original value (`X[1,1] == 1`, `Z[1,1] == 1`, `S[1,1] == 1`) after alias mutation.
- The mutated alias lost its resident binding (`.amatrix_resident_key(...) == NULL`) while the original retained its resident key.
- Dense alias replacement produced a distinct `object_id` for the mutated object, consistent with copy-on-modify rather than device-side aliasing.

## (d) Proposed bd create
- None.

## (e) Limitations
- Probes used the in-repo mock backend rather than a real GPU backend; this validates amatrix's alias/residency bookkeeping, not backend-specific device kernels.
- I did not probe alias mutation through every replacement path (`dimnames<-`, `diag<-`, logical subassignment) because the core resident-key alias risk did not reproduce.
