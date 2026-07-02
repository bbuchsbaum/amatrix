# amatrix package

Package-level roxygen directives for namespace generation.

## Startup message

On attach, amatrix prints a one-line note listing any installed GPU
backend packages (see
[`amatrix_use_gpu`](https://bbuchsbaum.github.io/amatrix/reference/amatrix_use_gpu.md)).
Pure-CPU sessions with no backend packages installed are silent. To
suppress the note explicitly, set
`options(amatrix.quiet_startup = TRUE)` or the environment variable
`AMATRIX_QUIET=1` (either `"1"` or `"true"`) before loading the package.
The note is also skipped when
`options(amatrix.optional_backends = FALSE)` disables all optional
backends, and any backend disabled via
`options(amatrix.disable_mlx = TRUE)` (or the analogous
`amatrix.disable_metal`, `amatrix.disable_opencl`,
`amatrix.disable_arrayfire`) is omitted from it.

## See also

Useful links:

- <https://bbuchsbaum.github.io/amatrix/>

- <https://github.com/bbuchsbaum/amatrix>

- Report bugs at <https://github.com/bbuchsbaum/amatrix/issues>

## Author

**Maintainer**: Brad Buchsbaum <brad.buchsbaum@gmail.com>
