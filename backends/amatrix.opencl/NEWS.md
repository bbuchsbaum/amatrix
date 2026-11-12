# amatrix.opencl 0.1.0

* First tagged release of the OpenCL/CLBlast backend for `amatrix`.
* OpenCL and CLBlast are now located at run time (`dlopen`/`LoadLibrary`)
  instead of being linked, so the package builds and installs on any
  platform, including Windows, and a missing runtime degrades to
  "unavailable" instead of an error or a process abort.
* New `amatrix_install_clblast()` fetches the official CLBlast binary where
  a system copy is not available (Homebrew hint on Apple Silicon).
* Explicit probe required (`amatrix_use_gpu()` or
  `AMATRIX_OPENCL_PROBE_GPU=1`).

## Linux and CI

* The package now links R's LAPACK/BLAS on Unix. It previously failed to load
  on Linux with `undefined symbol: dgebrd_` (macOS tolerated the missing link).
* `AMATRIX_OPENCL_ALLOW_CPU_DEVICE=1` lets the backend initialise on a CPU
  OpenCL device such as PoCL when no GPU is present. This gives no speed
  benefit; it exists so the backend's correctness can be tested on machines
  and CI runners without a GPU. GPU devices are still preferred.
