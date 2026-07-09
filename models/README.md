# models/ — not part of the v0.1.0 release

Decision 2026-07-09 (see `planning_docs/plan-release-readiness-2026-07.md`, D2):
`amatrix.models` is **out of release scope**. The supported modeling API lives in the
core package (`R/models-lm.R`: `many_lm`, `lm_fit`, `ridge_fit`, `wls_fit`, `covariance`,
`correlation`, `array_lm` wrappers). This scaffold is retained for a possible future
split-out; it is not on the R-universe registry and must not be added without a fresh
`R CMD check` pass and a reversal of D2.
