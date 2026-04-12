# cran-comments.md

## Test environments

* macOS 14 (arm64), R 4.5.1
* GitHub Actions: ubuntu-latest, macOS-latest, windows-latest (R-release, R-devel)

## R CMD check results

0 errors, 0 warnings, 1 note.

### Note on foreign function calls

```
checking foreign function calls ... NOTE
Foreign function calls to a different package:
  .Call("...", ..., PACKAGE = "amatrix.arrayfire")
  .Call("...", ..., PACKAGE = "amatrix.mlx")
  .Call("...", ..., PACKAGE = "amatrix.opencl")
```

These are intentional cross-package C bridge calls to the optional
backend packages (`amatrix.mlx`, `amatrix.opencl`, `amatrix.arrayfire`).
Each backend lives in a separate package so it can carry its own
system dependencies (MLX framework, OpenCL/CLBlast, ArrayFire) without
forcing them on users who only need the CPU path.

The bridge functions are exported from the backend packages and called
explicitly with `PACKAGE =` to ensure the right `.so` is used. The calls
are guarded by `requireNamespace()` so the core `amatrix` package
loads cleanly even when no backend is installed.

This is the same pattern used by `parallel`, `Rmpi`, `gbm` and other
packages that interface to optional companion packages.

## Downstream dependencies

This is a new package (initial CRAN submission). No reverse dependencies.
