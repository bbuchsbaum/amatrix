# amatrix.mlx 0.1.0

* First tagged release of the MLX backend for `amatrix` (Apple Silicon GPU).
* Dense float32 products, Cholesky/QR/SVD paths, and resident (GPU-memory)
  execution; activates automatically on first use via `amatrix`'s contained
  probe.
* Source builds succeed on any platform: without `mlx-c` a mock bridge is
  compiled and the backend reports itself unavailable.
