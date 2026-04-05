# ArrayFire OpenCL QR Crash Under R Extension

This directory contains an upstream-ready reproducer bundle for an ArrayFire QR crash observed from an R-loaded shared library on macOS.

## Summary

- `af_qr` succeeds in a standalone C executable.
- The same `af_qr` path crashes when called from a minimal R extension.
- The crash is backend-specific:
  - `OpenCL` backend: crashes at `96x96`
  - `CPU` backend: succeeds

This isolates the defect boundary to:

`ArrayFire QR + OpenCL backend + in-process R extension`

not to the `amatrix` package itself.

## Files

- [ISSUE_DRAFT.md](/Users/bbuchsbaum/code/amatrix/docs/repros/arrayfire-opencl-qr-r-extension/ISSUE_DRAFT.md)
- [af_r_repro.c](/Users/bbuchsbaum/code/amatrix/docs/repros/arrayfire-opencl-qr-r-extension/af_r_repro.c)
- [arrayfire_qr_repro.c](/Users/bbuchsbaum/code/amatrix/docs/repros/arrayfire-opencl-qr-r-extension/arrayfire_qr_repro.c)
- [run.sh](/Users/bbuchsbaum/code/amatrix/docs/repros/arrayfire-opencl-qr-r-extension/run.sh)

## Environment

- ArrayFire: `3.10.0`
- R: `4.5.1`
- macOS: `14.3`
- Active crashing backend in ArrayFire diagnostics: `4` (`OpenCL`)

## Expected Result

Both repros should either:

- succeed, or
- return a handled ArrayFire error

They should not segfault the R process.

## Observed Result

The standalone executable succeeds on `default`, `opencl`, and `cpu`.

The minimal R extension:

- crashes on `opencl`
- succeeds on `cpu`

## Current `amatrix` mitigation

`amatrix.arrayfire` now treats ArrayFire QR as safe only on the ArrayFire CPU backend by default, so the OpenCL crash path is not used accidentally.
