# Draft Upstream Issue: `af_qr` crashes on OpenCL when called from an R extension, but not from a standalone executable

## Summary

I can reproduce a crash in `af_qr` on macOS when ArrayFire is used from a minimal R shared library extension on the OpenCL backend.

The same QR sequence succeeds:

- in a standalone C executable
- in the same minimal R extension when the ArrayFire backend is forced to `CPU`

This suggests the defect is not in my package logic, and not a generic `af_qr` failure, but rather in the interaction between:

- ArrayFire OpenCL QR
- and an R-loaded in-process extension

## Environment

- ArrayFire: `3.10.0`
- macOS: `14.3`
- R: `4.5.1`
- Homebrew ArrayFire install
- `af_is_lapack_available()` returns `TRUE`

## Reproducer bundle

Included files:

- `arrayfire_qr_repro.c`
  - standalone executable
  - succeeds on `default`, `opencl`, and `cpu`
- `af_r_repro.c`
  - minimal R extension
  - crashes on `opencl`
  - succeeds on `cpu`

## Observed behavior

### Minimal R extension, OpenCL backend

The extension reaches:

```text
[af_r_qr] create_input
[af_r_qr] transpose_input
[af_r_qr] af_qr
```

and then the R process segfaults.

### Minimal R extension, CPU backend

The same extension succeeds:

- `af_qr`
- transpose `Q`
- transpose `R`
- materialize `Q`
- materialize `R`

### Standalone executable

The standalone executable succeeds on:

- `default`
- `opencl`
- `cpu`

for at least:

- `80x80`
- `96x96`
- `128x128`

## Expected behavior

`af_qr` should either:

- succeed, or
- return a handled `af_err`

It should not segfault the R process.

## Additional notes

- This is reproducible outside my package in a tiny R extension.
- In my package, I have mitigated it by treating ArrayFire QR as safe only on the ArrayFire CPU backend by default.

## Suggested next question

Is there a known issue with `af_qr` on the OpenCL backend when called from an embedded/shared-library host like R, or with OpenCL backend initialization/ownership in that setting?
