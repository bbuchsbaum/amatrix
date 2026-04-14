# amatrix.opencl — Implementation Plan (V1)

## Goal

Add a focused, reliable portable GPU backend for non-Apple platforms using
CLBlast (OpenCL BLAS) plus a small set of custom OpenCL kernels.

## Scope

Dense resident BLAS + elementwise/reduction kernels, plus backend-contract
parity for dense `solve`, `chol`, `qr`, and `covariance`. Those factor-style
ops are currently host-backed in V1, but they dispatch through the OpenCL
backend cleanly so higher-level model code can stay on the same execution path.

Sparse ops and exact dense spectral decompositions still fall back to CPU.
The intended OpenCL spectral story is iterative/subspace methods built from
dense products, not LAPACK-complete exact SVD on device. f32 only in V1.

## Why CLBlast

ArrayFire is a high-level broad library with integration risk (probe crashes,
sparse gaps, shape-blind dispatch — partly our integration bugs, partly upstream
fragility).  CLBlast is a lower-level focused BLAS building block: fewer moving
parts, well-tested GEMM, auto-tuning per device.  It is not a drop-in AF
replacement — it covers dense BLAS, not a tensor runtime.  The rest is our code.

Both projects are maintained (AF: Sep 2025 release, CLBlast: Mar 2026 release).

## Interface Contract

The backend object must be a third instance of the existing MLX/AF pattern.
Same function signatures, same resident lifecycle, same registration mechanism.
Model the interface by reading the MLX and AF backends line-by-line; the CLBlast
API maps to the *inside* of the bridges — everything R sees must be identical.

### Capabilities and Features

```r
capabilities = function() c("matmul", "crossprod", "tcrossprod", "ewise",
                             "broadcast_ewise", "rowSums", "colSums",
                             "solve", "chol", "qr", "covariance")
features     = function() c("dense_f32", "resident_dense", "custom_ops",
                             "solve", "chol", "qr", "covariance")
precision_modes = function() "fast"
```

## Op Mapping

### CLBlast-native ops

| amatrix op | CLBlast call | Notes |
|------------|-------------|-------|
| `matmul(A, B)` | `CLBlastSgemm(ColMajor, NoTrans, NoTrans, ...)` | Column-major throughout |
| `crossprod(x)` | `CLBlastSsyrk(ColMajor, Upper, Trans, ...)` | **Must fill lower triangle** before returning — SYRK only writes one half |
| `crossprod(x, y)` | `CLBlastSgemm(ColMajor, Trans, NoTrans, ...)` | Standard GEMM with transA |
| `tcrossprod(x)` | `CLBlastSsyrk(ColMajor, Lower, NoTrans, ...)` | **Must fill upper triangle** |
| `tcrossprod(x, y)` | `CLBlastSgemm(ColMajor, NoTrans, Trans, ...)` | Standard GEMM with transB |
| `ewise *` (matrix) | `CLBlastShad(n, alpha, x_buf, ..., y_buf, ..., z_buf, ...)` | Hadamard product |
| `ewise *` (scalar) | `CLBlastSscal(n, alpha, x_buf, ...)` on a copy | Or custom scalar_mul kernel |

### Custom OpenCL kernels (V1 — required for resident_handle)

These must ship in V1.  Without them, `broadcast_ewise_resident`,
`rowSums_resident`, and `colSums_resident` fall back to CPU, breaking the
resident chain that `resident_handle` depends on (see `R/resident-handle.R:96`
and `inst/examples/doubly-stochastic.R:150`).

| Kernel | Lines | Description |
|--------|-------|-------------|
| `ewise_add` | ~15 | `c[i] = a[i] + b[i]` or `c[i] = a[i] + scalar` |
| `ewise_sub` | ~15 | `c[i] = a[i] - b[i]` or `c[i] = a[i] - scalar` |
| `ewise_div` | ~15 | `c[i] = a[i] / b[i]` or `c[i] = a[i] / scalar` |
| `scalar_mul` | ~10 | `c[i] = a[i] * scalar` (when RHS is scalar, not matrix) |
| `broadcast_sweep` | ~25 | `c[i,j] = a[i,j] OP v[i]` (margin=1) or `v[j]` (margin=2) |
| `row_sum` | ~20 | Per-row reduction: `out[i] = sum_j a[i,j]` |
| `col_sum` | ~20 | Per-column reduction: `out[j] = sum_i a[i,j]` |
| `sym_fill` | ~10 | Copy upper→lower or lower→upper triangle (for SYRK results) |

Kernels are compiled once at first use via `clCreateProgramWithSource` +
`clBuildProgram`, cached as static `cl_program`/`cl_kernel` handles.

### CPU fallback (V1)

Exact dense `svd`, backend-native `rsvd`, `eigen`, all sparse ops, `argmax`,
`scatter_mean`, `segment_sum`, `segment_mean`, and batched GEMM. Dense
`solve`/`chol`/`qr`/`covariance` are exposed through the OpenCL backend
contract but currently reuse host math internally.

### Spectral posture (V1)

OpenCL should be treated as an iterative spectral backend, not an exact-SVD
backend. The intended path is:

- `block_lanczos()` on fast OpenCL matrices
- `svd_factor(method = "auto" | "subspace")` choosing OpenCL when the matrix is
  large enough and the requested rank is moderate

Exact `am_svd()` remains CPU-backed unless a future backend-native dense SVD is
added.

## Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Column order | `CLBlastLayoutColMajor` | Match R, avoid AF's historical row-major bug |
| Sync model | `clFinish(queue)` after every op | Simple for V1; async in V2 |
| Precision | f32 only | Avoids runtime `cl_khr_fp64` probing; matches AF/MLX posture |
| Probe safety | Env gate `AMATRIX_OPENCL_PROBE_GPU=1` | Follow MLX pattern; default to safe (unavailable) |
| Device selection | First available GPU | Multi-device in V2 |
| Self-crossprod | SYRK + sym_fill kernel | Must fill untouched triangle — amatrix expects full matrix |
| Scalar ewise | Custom kernel or CLBlastSscal on copy | Conformance suite expects `A * 3.0` on all backends |
| Registration | Coordinated core + backend | Add to `.amatrix_optional_backend_specs()`, gated with `amatrix.enable_opencl` |

## Package Structure

```
backends/amatrix.opencl/
  DESCRIPTION
  NAMESPACE
  configure                 # pkg-config --cflags --libs clblast; detect OpenCL
  src/
    Makevars.in             # @OPENCL_CFLAGS@ @CLBLAST_CFLAGS@ @LIBS@
    opencl_bridge.c         # ~900 LOC: context/queue, registry, CLBlast wrappers
    opencl_kernels.c        # ~350 LOC: custom kernel compilation + dispatch
    init.c                  # ~80 LOC: R .Call registration
  R/
    backend.R               # ~450 LOC: backend object, capabilities, R wrappers
    zzz.R                   # ~15 LOC: .onLoad — conditional, NOT eager
  inst/
    PROBE_SAFETY.md         # doc: why env gate exists, how to enable
```

## C Bridge Design

### Global state

```c
static cl_platform_id   g_platform = NULL;
static cl_device_id     g_device   = NULL;
static cl_context       g_context  = NULL;
static cl_command_queue  g_queue    = NULL;
static int              g_initialized = 0;
```

### Init (gated by env var)

```c
SEXP amatrix_opencl_native_available_bridge(void) {
#ifdef HAVE_CLBLAST
  const char *probe = getenv("AMATRIX_OPENCL_PROBE_GPU");
  if (probe == NULL || strcmp(probe, "1") != 0)
    return ScalarLogical(0);  // safe default
  // clGetPlatformIDs → clGetDeviceIDs → clCreateContext → clCreateCommandQueue
  // All error-checked; return ScalarLogical(0) on any failure
#else
  return ScalarLogical(0);
#endif
}
```

### Resident registry (~256 slots)

```c
typedef struct {
  char key[64];
  cl_mem buffer;
  int nrow, ncol;
  int in_use;
} ocl_resident_entry;

// store:  clCreateBuffer(CL_MEM_READ_WRITE) + clEnqueueWriteBuffer (f64→f32 + col-major)
// has:    scan registry for key
// drop:   clReleaseMemObject + clear slot
// materialize: clEnqueueReadBuffer + f32→f64 conversion
```

Reserve slot AFTER successful `clCreateBuffer` + `clEnqueueWriteBuffer`
(unlike AF which reserves first — OOM leaves ghost entries).

### GEMM wrapper

```c
// CLBlastSgemm(CLBlastLayoutColMajor, transA, transB,
//              m, n, k, alpha, a_buf, lda, b_buf, ldb,
//              beta, c_buf, ldc, &g_queue, &event);
// clFinish(g_queue);  // sync V1
```

### SYRK + sym_fill

```c
// CLBlastSsyrk(CLBlastLayoutColMajor, CLBlastTriangleUpper, CLBlastTransposeYes,
//              n, k, alpha, a_buf, lda, beta, c_buf, ldc, &g_queue, &event);
// clFinish(g_queue);
// enqueue sym_fill kernel to copy upper→lower triangle
// clFinish(g_queue);
```

## R Backend Object

Matches MLX/AF pattern exactly.  Every `*_resident` function includes
`defer = FALSE` parameter.

```r
amatrix_opencl_backend <- function() {
  cpu <- .amatrix_cpu_backend()
  list(
    capabilities    = function() c("matmul", "crossprod", "tcrossprod",
                                   "ewise", "broadcast_ewise",
                                   "rowSums", "colSums"),
    features        = function() c("dense_f32", "resident_dense", "custom_ops"),
    precision_modes = function() "fast",
    available       = function() amatrix_opencl_is_available(),
    supports = function(op, x, y = NULL) {
      if (!inherits(x, "adgeMatrix")) return(FALSE)
      if (!x@precision %in% "fast") return(FALSE)
      if (!op %in% c("matmul","crossprod","tcrossprod","ewise",
                      "broadcast_ewise","rowSums","colSums")) return(FALSE)
      dims <- dim(x)
      min(dims) >= getOption("amatrix.opencl.matmul_min_dim", 128L)
    },
    # Cold ops
    matmul     = function(x, y)    amatrix_opencl_matmul(x, y),
    crossprod  = function(x, y)    amatrix_opencl_crossprod(x, y),
    tcrossprod = function(x, y)    amatrix_opencl_tcrossprod(x, y),
    ewise      = function(x, lhs, rhs, op, ...) amatrix_opencl_ewise(lhs, rhs, op),
    broadcast_ewise = function(x, lhs, v, margin, op, ...)
      base::sweep(as.matrix(lhs), MARGIN = margin, STATS = v, FUN = op),
    rowSums = function(x, na.rm = FALSE, dims = 1L) amatrix_opencl_axis_sums(x, 1L),
    colSums = function(x, na.rm = FALSE, dims = 1L) amatrix_opencl_axis_sums(x, 0L),
    # Residency
    resident_store       = function(key, x) amatrix_opencl_resident_store(key, x),
    resident_has         = function(key)    amatrix_opencl_resident_has(key),
    resident_drop        = function(key)    amatrix_opencl_resident_drop(key),
    resident_materialize = function(key)    amatrix_opencl_resident_materialize(key),
    # Resident ops (all with defer)
    matmul_resident   = function(x_key, y_key, out_key, defer = FALSE) ...,
    crossprod_resident  = function(x_key, y_key = NULL, out_key, defer = FALSE) ...,
    tcrossprod_resident = function(x_key, y_key = NULL, out_key, defer = FALSE) ...,
    ewise_resident      = function(lhs_key, rhs, op, out_key, defer = FALSE) ...,
    broadcast_ewise_resident = function(lhs_key, v, margin, op, out_key, defer = FALSE) ...,
    rowSums_resident = function(x_key, na.rm = FALSE, dims = 1L) ...,
    colSums_resident = function(x_key, na.rm = FALSE, dims = 1L) ...
  )
}
```

## Registration

Coordinated core + backend change:

**In core `amatrix` package** (`R/backend-registry.R`):
```r
# Add to .amatrix_optional_backend_specs():
opencl = list(
  package        = "amatrix.opencl",
  register_fun   = "amatrix_opencl_register",
  enable_option  = "amatrix.enable_opencl",
  available_fun  = "amatrix_opencl_is_available"
)
```

**In `amatrix.opencl/R/zzz.R`:**
```r
.onLoad <- function(libname, pkgname) {
  # Conditional — only register if explicitly enabled
  if (isTRUE(getOption("amatrix.enable_opencl", FALSE))) {
    amatrix_opencl_register(overwrite = FALSE)
  }
}
```

**Probe function:**
```r
amatrix_opencl_enable_probe <- function(register = TRUE) {
  Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = "1")
  if (register) amatrix_opencl_register(overwrite = TRUE)
}
```

## Build Requirements

User needs: OpenCL ICD loader + CLBlast.

```sh
# Detection in configure:
pkg-config --cflags --libs clblast    # CLBlast headers + lib
# Fallback: -I/usr/local/include -L/usr/local/lib -lclblast -lOpenCL
```

`Makevars.in`:
```
PKG_CPPFLAGS = -DHAVE_CLBLAST @CLBLAST_CFLAGS@
PKG_LIBS = @CLBLAST_LIBS@ -lOpenCL
```

When CLBlast is not found: all bridges compile as stubs (return `ScalarLogical(0)`
or `error("requires CLBlast")`).  Package installs and loads cleanly.

## Testing

### Must-have V1 tests

1. **Cross-backend conformance** — reuse `test-cross-backend-conformance.R` with
   `.GPU_TOL = 1e-4` for all claimed ops.

2. **Non-square crossprod/tcrossprod suite** — dedicated tests for rectangular
   inputs (m != n), verifying SYRK + sym_fill produces full symmetric output
   matching `base::crossprod()`.

3. **Scalar ewise** — `A * 3.0`, `A + 1.0`, `A / 2.0` matching conformance
   expectations (`test-cross-backend-conformance.R:78`).

4. **Probe safety** — `library(amatrix.opencl)` + `amatrix_backend_names()` do
   not crash when OpenCL is absent or misconfigured.  Test on a system without
   OpenCL installed.

5. **Resident round-trip** — store → materialize matches original within f32 tol.

6. **Resident chain** — store A, store B, `matmul_resident` → materialize →
   matches `A %*% B`.

7. **resident_handle Sinkhorn** — 50 iterations of doubly-stochastic via handle
   on opencl backend produces correct result (row/col sums within 1e-4).

### Primary test platform

Linux + NVIDIA (most common HPC GPU, best OpenCL runtime quality).
Performance baselines are machine-specific (not CI) per `planning_docs/quality-tracking.md`.

## Estimated Scope

| Component | LOC | Complexity |
|-----------|-----|------------|
| `opencl_bridge.c` (context, registry, CLBlast) | ~900 | Medium |
| `opencl_kernels.c` (7 custom kernels) | ~350 | Medium |
| `init.c` | ~80 | Low |
| `configure` | ~70 | Medium |
| `R/backend.R` | ~450 | Low (pattern from AF/MLX) |
| `R/zzz.R` | ~15 | Low |
| Core `backend-registry.R` change | ~5 | Low |
| Core `helper-conformance.R` change | ~10 | Low |
| Tests | ~150 | Low |
| **Total** | **~2030** | |

## V2 Roadmap (not V1)

- f64 support (runtime `cl_khr_fp64` detection, `precision_modes = c("strict", "fast")`)
- Async execution (`clWaitForEvents` instead of `clFinish`)
- Batched GEMM (`CLBlastSgemmBatched`)
- Multi-device selection
- solve/chol on GPU via custom kernels or MAGMA integration
- Linux + AMD ROCm testing
