Sys.setenv(AMATRIX_METAL_PROBE_GPU = "1")

# The suite exercises dispatch mechanics on small matrices. Switch off the
# default "keep small / transfer-bound work on the CPU when uncalibrated" rule
# so those tests still reach the backend under test; the rule itself is
# covered by its own tests.
options(amatrix.uncalibrated_min_elements = 0)
