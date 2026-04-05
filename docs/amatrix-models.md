# amatrix.models

`amatrix.models` is the intended package boundary for model-core workflows built on top of the `amatrix` substrate.

Current center of gravity:

- `many_lm(..., method = "qr")`

If you are evaluating the package for practical adoption, start there. That is the workflow surface where the current QR/residency work is paying off most clearly.

Current exported functions:

- `lm_fit(...)`
- `ridge_fit(...)`
- `wls_fit(...)`
- `many_lm(...)`
- `array_lm(...)`
- `covariance(...)`
- `correlation(...)`

These are thin wrappers over the current `amatrix` model-core functions:

- `amatrix::am_lm_fit()`
- `amatrix::am_ridge_fit()`
- `amatrix::am_wls_fit()`
- `amatrix::am_many_lm()`
- `amatrix::am_array_lm()`
- `amatrix::am_covariance()`
- `amatrix::am_correlation()`

The purpose of the package is not to duplicate the substrate. It is to give downstream users a model-oriented namespace with ordinary function names while `amatrix` remains centered on objects, kernels, dispatch, and backend policy.

Current migration rule:

- `amatrix.models` is the preferred user-facing namespace for model-core workflows.
- `amatrix` still contains the implementation today to avoid a circular dependency and to keep kernel/model integration simple while the API is still settling.
- Further migration should happen only when the model-core surface is stable enough that moving implementation out of `amatrix` reduces complexity rather than increasing it.

This is the current package split direction:

- `amatrix`
  - substrate
  - kernels
  - backend policy
- `amatrix.models`
  - repeated-response least squares, with `many_lm(..., method = "qr")` as the flagship workflow
  - ridge
  - array-response regression
  - covariance and correlation helpers
