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
