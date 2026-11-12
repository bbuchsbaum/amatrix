# Tests run via testthat / devtools::test() / R CMD check use the -e launch
# mode, which is safe for Metal device initialisation.  Activate the GPU probe
# so that skip_if_not(amatrix_mlx_is_available()) reflects actual hardware.
Sys.setenv(AMATRIX_MLX_PROBE_GPU = "1")
Sys.setenv(AMATRIX_OPENCL_PROBE_GPU = "1")

# The suite exercises dispatch mechanics on small matrices. Switch off the
# default "keep small / transfer-bound work on the CPU when uncalibrated" rule
# so those tests still reach the backend under test; the rule itself is
# covered by its own tests.
options(amatrix.uncalibrated_min_elements = 0)
