# Virtual base class for backend-aware matrices

`aMatrix` is the abstract base from which all concrete amatrix classes
inherit. It carries backend-dispatch metadata that controls which
compute backend (CPU, GPU, etc.) is used for operations on the matrix.

## Slots

- `preferred_backend`:

  Single string naming the preferred compute backend; one of `"cpu"`,
  `"mlx"`, `"metal"`, `"arrayfire"`, or `"torch"`.

- `policy`:

  Single string controlling dispatch policy; one of `"auto"`, `"cpu"`,
  `"mlx"`, `"metal"`, `"arrayfire"`, or `"torch"`.

- `precision`:

  Single string; either `"strict"` (double precision, exact results) or
  `"fast"` (backend may use lower precision).

- `object_id`:

  Non-empty string uniquely identifying this object within the session;
  used for caching and residency tracking.

- `src_id`:

  String recording the `object_id` of the object this was derived from,
  or `""` for originals.

- `finalizer_env`:

  Environment used to manage GPU-resident memory and deferred host-copy
  state.
