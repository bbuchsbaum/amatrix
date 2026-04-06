# amatrix Backend Contract

This document defines the v1 backend plugin contract for `amatrix`.

The contract is intentionally narrow. Backends are optional plugins that
accelerate a small dense-first operation surface. Public materialized matrix
values remain Matrix-compatible S4 objects. Specialized S4 view and factor
intermediates are allowed when they preserve ordinary matrix semantics without
forcing host materialization. Sparse semantics are CPU-native first and
`Matrix` remains the official sparse backend in v1.

## v1 Rules

- User-visible materialized matrix values are always valid Matrix-compatible S4
  objects first.
- Specialized S4 view or factor intermediates are allowed when they preserve
  ordinary matrix semantics.
- Backends operate on host-side inputs per cold call.
- Backends do not own persistent device state in ordinary materialized matrix
  objects.
- Sparse is `Matrix`-first in v1.
- GPU sparse is opportunistic and limited to explicitly supported kernels later.
- Unsupported operations must fall back predictably to CPU semantics.

## Cold vs Resident Execution

The current core now distinguishes between two dense execution regimes:

- `cold`
  - A fresh host-side object is evaluated for one operation.
  - Backend selection is based on ordinary `supports(op, x, y)`.
- `resident`
  - A dense object already has live backend-resident state.
  - The backend may be selected even when the cold path would not be chosen, as long as it exposes the corresponding `*_resident` implementation.

This distinction matters because the endgame is not "one fast kernel." It is chained execution where backend residency avoids repeated upload and conversion costs.

The current public helpers surface this directly:

- `amatrix_backend_plan(x, op, y = NULL)`
  - returns `chosen_path`, which is currently `"cold"` or `"resident"`
- `amatrix_backend_matrix(x, ops, y_map = list())`
  - returns `chosen_path` and `resident_reuse` for each operation
- `amatrix_execution_info(x, ops, y_map = list())`
  - returns object metadata, residency state, and per-op planning in one object

## Structural Views And Explicit BLAS Control

The product syntax contract is intentionally asymmetric:

- Common structural cases should stay elegant under ordinary R syntax.
  - `%*%` remains the ordinary product entry point.
  - `crossprod()` and `tcrossprod()` remain first-class transpose-product APIs.
  - The current code uses a `src_id`-linked transpose shortcut as a stepping
    stone. The cleaner target is a dedicated zero-copy transpose view.
- Full DGEMM control is explicit.
  - `alpha`, `beta`, an accumulator `C`, and explicit transpose flags belong in
    an `am_gemm()`-style kernel, not in a general lazy expression system behind
    `+` and `*`.

This keeps the object surface boring while still leaving room for aggressive
kernel implementations.

## Required Backend Shape

A backend is a named list registered with `amatrix_register_backend(name, backend, overwrite = FALSE)`.

The list must contain these required fields:

```r
list(
  capabilities = function() character(),
  features = function() character(),
  precision_modes = function() c("strict", "fast"),
  available = function() TRUE,
  supports = function(op, x, y = NULL) FALSE,
  matmul = function(x, y) NULL,
  crossprod = function(x, y = NULL, ...) NULL,
  tcrossprod = function(x, y = NULL, ...) NULL,
  ewise = function(x, lhs, rhs = NULL, op, ...) NULL,
  rowSums = function(x, na.rm = FALSE, dims = 1L) NULL,
  colSums = function(x, na.rm = FALSE, dims = 1L) NULL
)
```

## Planned Optional Product Extension

The current required v1 contract does not include GEMM as a separate backend
entry point. The intended explicit power-user extension is:

```r
gemm = function(x, y, C = NULL,
                alpha = 1, beta = 0,
                transA = FALSE, transB = FALSE) NULL
```

and, where residency is supported:

```r
gemm_resident = function(x_key, y_key, c_key = NULL, out_key,
                         alpha = 1, beta = 0,
                         transA = FALSE, transB = FALSE) NULL
```

This is planned optional surface, not part of the current minimal backend
registration contract.

### Required Fields

- `capabilities()`
  - Returns a character vector naming the operations the backend is designed to implement.
  - The order is preserved and should be deliberate. `amatrix` does not sort it away.
  - In v1 this is descriptive metadata plus a consistency target for `supports()`.
- `features()`
  - Returns a character vector describing backend traits beyond individual operation names.
  - This is the stable capability vocabulary that higher-level planning and model-core code should use when backend-name checks would otherwise leak into the codebase.
  - Current vocabulary in core is:
    - `dense_f64`
    - `dense_f32`
    - `resident_dense`
    - `unified_memory`
    - `solve`
    - `chol`
    - `svd`
    - `sparse_spmm`
    - `custom_ops`
- `precision_modes()`
  - Returns the precision policies the backend can honestly support.
  - Current required values are a subset of:
    - `strict`
    - `fast`
  - A backend that only has float32-oriented accelerator paths should advertise `fast`, not `strict`.
- `available()`
  - Returns a single `TRUE` or `FALSE`.
  - `FALSE` means the backend is installed or registered but not usable on the current machine or configuration.
  - Example causes: missing native library, unsupported platform, disabled build flag.
- `supports(op, x, y = NULL)`
  - Returns whether the backend can handle this specific call.
  - This is where class restrictions, dtype restrictions, shape restrictions, and operand restrictions belong.
  - `supports()` should be stricter than `capabilities()` when necessary.
- `matmul()`, `crossprod()`, `tcrossprod()`, `ewise()`, `rowSums()`, `colSums()`
  - These are the callable implementations for the dense-first v1 accelerated surface.
  - They receive host-side inputs, not device handles stored on the object.

## Current Operation Vocabulary

The core package currently reasons about these operation names:

```r
c(
  "matmul", "crossprod", "tcrossprod",
  "ewise", "rowSums", "colSums",
  "solve", "chol", "qr", "svd", "eigen", "diag"
)
```

The CPU backend advertises all of them. Optional backends currently only need the dense-first subset:

```r
c("matmul", "crossprod", "tcrossprod", "ewise", "rowSums", "colSums")
```

Do not advertise `solve`, `chol`, `qr`, `svd`, `eigen`, or `diag` unless the backend actually implements them and the semantics are covered by parity tests.

## Dense-First, Sparse-First

The architectural split for v1 is:

- Dense acceleration is the purpose of optional backends.
- Sparse execution defaults to `Matrix` and CPU paths.
- Public sparse objects stay CSC-oriented and `dgCMatrix`-compatible.

That means:

- A dense backend may support `adgeMatrix` for `matmul` and `ewise`.
- The same backend should usually reject `adgCMatrix` in `supports()`.
- Sparse `%*%`, sparse factorization, sparse indexing, sparse coercion, and structural sparse operations should continue to route through `Matrix` unless the project explicitly adds a sparse kernel later.

This is the rule to preserve:

> Sparse is CPU-native first; GPU sparse is opportunistic.

## Selection Precedence

Backend selection is driven by `amatrix_backend_plan(x, op, y = NULL)`.

Candidate order is:

1. `x@preferred_backend`
2. `x@policy`
3. `amatrix_default_policy()`
4. `"cpu"`

The first candidate that is:

- registered,
- available,
- precision-compatible,
- and reports `supports(op, x, y) == TRUE`

is chosen.

If no candidate qualifies, `cpu` is chosen as the fallback.

This means:

- `preferred_backend` wins over policy fallback.
- `policy` only matters if the preferred backend is unavailable or unsupported.
- CPU is always the semantic safety net.

## Introspection Helpers

The core package exposes helpers that backend authors should use while developing:

- `amatrix_backend_names()`
  - Lists registered backends.
- `amatrix_backend_capabilities(name)`
  - Returns the backend's declared capability vector in declared order.
- `amatrix_backend_features(name)`
  - Returns the backend's declared feature vocabulary in declared order.
- `amatrix_backend_status(names = amatrix_backend_names())`
  - Returns a data frame of availability, precision modes, backend features, residency capability, and operation summaries.
- `amatrix_backend_plan(x, op, y = NULL)`
  - Explains which backend would be chosen for one call, and whether that choice is a cold path or resident reuse.
- `amatrix_backend_matrix(x, ops, y_map = list())`
  - Summarizes chosen backends across several operations, including resident reuse.
- `amatrix_execution_info(x, ops, y_map = list())`
  - Returns a compact execution summary for one object, including residency state and per-op plans.

These helpers are the fastest way to debug selection behavior before adding native code.

## Minimal Dense Backend Skeleton

This is the intended v1 shape for a backend package:

```r
mybackend_capabilities <- function() {
  c("matmul", "crossprod", "tcrossprod", "ewise", "rowSums", "colSums")
}

mybackend_features <- function() {
  c("dense_f64")
}

mybackend_precision_modes <- function() {
  c("strict", "fast")
}

mybackend_is_available <- function() {
  TRUE
}

mybackend_backend <- function() {
  cpu <- amatrix:::.amatrix_cpu_backend()
  capabilities <- mybackend_capabilities()

  list(
    capabilities = function() capabilities,
    features = function() mybackend_features(),
    precision_modes = function() mybackend_precision_modes(),
    available = function() mybackend_is_available(),
    supports = function(op, x, y = NULL) {
      methods::is(x, "adgeMatrix") &&
        x@precision %in% mybackend_precision_modes() &&
        op %in% capabilities
    },
    matmul = function(x, y) {
      cpu$matmul(x, y)
    },
    crossprod = function(x, y = NULL, ...) {
      cpu$crossprod(x, y = y, ...)
    },
    tcrossprod = function(x, y = NULL, ...) {
      cpu$tcrossprod(x, y = y, ...)
    },
    ewise = function(x, lhs, rhs = NULL, op, ...) {
      cpu$ewise(x, lhs = lhs, rhs = rhs, op = op, ...)
    },
    rowSums = function(x, na.rm = FALSE, dims = 1L) {
      cpu$rowSums(x, na.rm = na.rm, dims = dims)
    },
    colSums = function(x, na.rm = FALSE, dims = 1L) {
      cpu$colSums(x, na.rm = na.rm, dims = dims)
    }
  )
}

mybackend_register <- function(overwrite = TRUE) {
  amatrix::amatrix_register_backend("mybackend", mybackend_backend(), overwrite = overwrite)
  invisible("mybackend")
}
```

This is intentionally close to the existing `amatrix.mlx` and `amatrix.arrayfire` stubs.

## Correct Fallback Behavior

A backend should reject unsupported calls early in `supports()`. The core package will then choose CPU and preserve object semantics.

Example:

```r
backend$supports <- function(op, x, y = NULL) {
  methods::is(x, "adgeMatrix") && op %in% c("matmul", "ewise")
}
```

Then:

```r
x <- adgeMatrix(matrix(c(4, 1, 1, 3), nrow = 2), preferred_backend = "mybackend")

amatrix_backend_plan(x, "matmul", y = diag(2))$chosen
# "mybackend"

amatrix_backend_plan(x, "solve")$chosen
# "cpu"
```

That fallback is correct and expected. Backends should not fake support for operations they cannot implement semantically.

## Return Values and Semantics

In v1, backend methods should behave as if they are specialized host-side implementations.

- Accept host-side `Matrix`-like inputs.
- Return results compatible with the CPU path.
- Leave serialization, copying, and subassignment semantics unchanged from ordinary host-backed objects.

The core wrappers rewrap numeric matrix-like results back into `amatrix` classes when appropriate. Backends do not need to construct `adgeMatrix` or `adgCMatrix` directly unless there is a strong reason.

## Residency in Current Backends

Current backend posture is intentionally asymmetric:

- `amatrix.mlx`
  - is the first backend with dense resident chaining
  - is the current path toward a transparent GPU-backed user experience
- `amatrix.arrayfire`
  - now exposes resident hooks for dense objects as well
  - remains most compelling today for portable dense kernels such as `matmul`
  - is still behind MLX in product validation for chained execution and flagship workflow documentation

That is still an intentional product posture. The MLX path is the one being pushed toward transparency first. ArrayFire remains a correct portable backend, but MLX is the more fully validated transparent resident backend today.

## What Not To Do in v1

Do not do these yet:

- store external device handles in user-visible `amatrix` objects
- require ALTREP for correctness
- advertise sparse support unless the sparse kernel is explicitly implemented and tested
- silently densify large sparse inputs for `svd`, `eigen`, or similar operations
- widen the public contract before the parity harness covers the new surface

## Current Reference Implementations

Use these packages as reference shapes:

- [R/backend.R](/Users/bbuchsbaum/code/amatrix/R/backend.R)
- [backends/amatrix.mlx/R/backend.R](/Users/bbuchsbaum/code/amatrix/backends/amatrix.mlx/R/backend.R)
- [backends/amatrix.arrayfire/R/backend.R](/Users/bbuchsbaum/code/amatrix/backends/amatrix.arrayfire/R/backend.R)

They are stubs, but they define the current contract accurately.
