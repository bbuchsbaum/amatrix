# amatrix.metal 0.1.0

* First tagged release of the experimental Metal sparse backend for `amatrix`.
* Sparse matrix times dense matrix products on macOS via a direct
  Objective-C++ bridge; explicit probe required (`amatrix_use_gpu()`).
* Builds a plain-C++ mock bridge on platforms without the Metal frameworks.
