# CRAN Release Plan for `amatrix` With Optional Backends

## Bottom Line

`amatrix` can go to CRAN without pulling the backend packages out of the
monorepo, but the nested backend packages cannot function as bundled runtime
dependencies of the CRAN release.

The CRAN package should ship as the CPU-first core. Optional backends should
remain separate packages and be discovered conditionally at runtime.

## Current Repo Shape

- Main CRAN-facing package: `amatrix`
- Optional backend packages:
  - `backends/amatrix.mlx`
  - `backends/amatrix.opencl`
  - `backends/amatrix.arrayfire`
  - `backends/amatrix.metal`
- `.Rbuildignore` excludes `^backends$`, so the backend package directories are
  not included in the `amatrix` source tarball.
- `DESCRIPTION` currently keeps backend packages in `Suggests`, not `Imports`.

This is the correct structural direction for CRAN.

## CRAN Constraints That Matter

1. Packages in `Depends`, `Imports`, or `LinkingTo` should be available from
   CRAN or Bioconductor. Nested packages inside the repo do not satisfy this.
2. Non-CRAN optional packages may appear in `Suggests` only if they are used
   conditionally and do not block install, load, or checks.
3. If non-CRAN packages remain in `Suggests`, provide an
   `Additional_repositories` entry or another clear installation route.
4. The CRAN package must remain fully functional without any GPU backend
   installed.

## Recommended Release Architecture

### 1. Release `amatrix` to CRAN as the CPU core

- Treat CPU as the authoritative backend.
- Ensure all examples, checks, and documentation succeed with no optional
  backends installed.
- Make correctness and API guarantees derive from the CPU path.

### 2. Keep backends as separate optional packages

- Do not move backend packages into `Imports`.
- Keep backend loading conditional with `requireNamespace()` or
  `loadNamespace()`.
- Allow `amatrix` to auto-register optional backends only when they are already
  installed.

### 3. Distribute backends independently

- Keep them in the monorepo if desired.
- Build and release each backend from its own subdirectory.
- Prefer hosting them outside CRAN first, e.g. via `r-universe` or a private
  package repository.

### 4. Add repository metadata only if needed

- If backend packages remain in `Suggests` for the CRAN package, add an
  `Additional_repositories` field pointing to the repository that serves them.
- If that creates CRAN friction, remove the backend packages from `Suggests`
  entirely and document them as optional post-install extensions.

## What Not To Do

- Do not expect CRAN to install `backends/amatrix.*` from inside the `amatrix`
  tarball.
- Do not make `amatrix` installation depend on MLX, OpenCL, ArrayFire, or Metal
  being present.
- Do not move non-CRAN backend packages into strong dependency fields.
- Do not let optional backends leak into `.onLoad()` or examples in a way that
  makes CRAN checks environment-dependent.

## Concrete Pre-CRAN Checklist

1. Audit `DESCRIPTION`.
   - Keep backend packages out of `Imports`.
   - Decide whether they stay in `Suggests` or move out entirely for the CRAN
     release.
2. Audit namespace loading.
   - Ensure backend discovery is conditional and failure-safe.
   - Ensure `amatrix` loads cleanly when none of the optional backend packages
     are installed.
3. Audit examples, vignettes, and tests.
   - Skip backend-specific checks unless the backend package is installed and
     healthy.
   - Ensure CRAN-facing examples are CPU-safe and deterministic.
4. Decide distribution channel for backends.
   - `r-universe` is the simplest realistic choice.
5. If keeping backend packages in `Suggests`, add `Additional_repositories`.
6. Write installation guidance.
   - `install.packages("amatrix")`
   - optional backend install instructions separately

## Recommended Decision

Ship `amatrix` on CRAN as a standalone CPU-first package and keep all GPU
backends as separately installable optional packages. Keep the monorepo layout,
but do not treat nested backend packages as part of the CRAN dependency graph.

