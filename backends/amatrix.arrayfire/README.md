# amatrix.arrayfire

ArrayFire backend for [amatrix](https://github.com/bbuchsbaum/amatrix). It
registers a backend-contract implementation and dispatches dense matmul to the
ArrayFire C API when ArrayFire is available at build time.

## Building

The `configure` script locates ArrayFire in this order:

1. `ARRAYFIRE_PREFIX` (user override)
2. `pkg-config --exists arrayfire`
3. Homebrew (`brew --prefix arrayfire`)
4. Common system paths (`/usr/local`, `/usr`, `/opt/arrayfire`, `/opt/local`)

If ArrayFire is not found the package still builds, but as a stub whose native
path is unavailable at runtime. A Windows path is provided via `configure.win`.

## Runtime backend selection (Apple Silicon caveat)

ArrayFire ships several compute runtimes (CPU, OpenCL, CUDA, oneAPI) selected at
runtime. On arm64 macOS, ArrayFire defaults to its OpenCL runtime, which aborts
inside `clGetDeviceIDs` on Apple Silicon hosts. To keep the backend usable,
`amatrix.arrayfire` pins the runtime to ArrayFire's **CPU** backend on Apple
Silicon unless you explicitly request another runtime (see
`R/backend.R`, `.amatrix_arrayfire_configured_runtime_backend()`).

Override the runtime with either the environment variable or the R option
(the environment variable takes precedence):

```r
Sys.setenv(AMATRIX_ARRAYFIRE_BACKEND = "opencl")   # or "cpu", "cuda", "oneapi"
options(amatrix.arrayfire.backend = "cpu")
```

Unrecognized values are ignored with a warning. On non-Apple-Silicon platforms
no runtime is pinned and ArrayFire's own default applies.

GPU probing is gated behind `AMATRIX_ARRAYFIRE_PROBE_GPU=1` so that a default
session does not touch the GPU.
