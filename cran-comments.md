# cran-comments.md

## Test environments

* macOS 14 (arm64), R 4.5.x (local)
* GitHub Actions: ubuntu-latest, macOS-latest, windows-latest (R-release, R-devel)

## R CMD check results

0 errors, 0 warnings, 1 note (`R CMD check --as-cran`, 2026-07-24).

### Note on CRAN incoming feasibility

```
New submission

Suggests or Enhances not in mainstream repositories:
  amatrix.mlx, amatrix.arrayfire, amatrix.opencl, amatrix.metal
Availability using Additional_repositories specification:
  amatrix.mlx         yes   https://bbuchsbaum.r-universe.dev
  amatrix.arrayfire   yes   https://bbuchsbaum.r-universe.dev
  amatrix.opencl      yes   https://bbuchsbaum.r-universe.dev
  amatrix.metal       yes   https://bbuchsbaum.r-universe.dev
```

The four suggested backend packages provide optional GPU acceleration and
live in separate packages so each can carry its own system dependencies
(MLX framework, OpenCL/CLBlast, ArrayFire, Metal) without forcing them on
users who only need the CPU path. They are available from the
`Additional_repositories` R-universe listed in DESCRIPTION. The core
package loads, runs, and passes its full test suite with none of them
installed; all GPU-dependent tests skip with a reason.

## Downstream dependencies

This is a new package (initial CRAN submission). No reverse dependencies.
