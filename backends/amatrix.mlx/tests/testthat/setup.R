# Tests run via testthat / devtools::test() / R CMD check use the -e launch
# mode, which is safe for Metal device initialisation.  Activate the GPU probe
# so that amatrix_mlx_is_available() reflects actual hardware.
Sys.setenv(AMATRIX_MLX_PROBE_GPU = "1")
