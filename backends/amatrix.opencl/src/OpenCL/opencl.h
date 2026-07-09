/*
 * Minimal <OpenCL/opencl.h> compatibility shim provided by amatrix.opencl.
 *
 * On Apple platforms the vendored CLBlast C header includes
 * <OpenCL/opencl.h> (the macOS framework umbrella header). amatrix.opencl
 * resolves every OpenCL and CLBlast symbol at runtime via dlopen, so it does
 * not link against or compile against the system OpenCL.framework. This shim
 * redirects that include to the vendored Khronos core C API headers so the
 * build is identical and self-contained on every platform.
 *
 * This file is authored by the amatrix project (MIT), not by Apple or The
 * Khronos Group. See inst/COPYRIGHTS.
 */
#ifndef AMATRIX_OPENCL_OPENCL_SHIM_H
#define AMATRIX_OPENCL_OPENCL_SHIM_H

#ifndef CL_TARGET_OPENCL_VERSION
#define CL_TARGET_OPENCL_VERSION 120
#endif

#include <CL/cl.h>

#endif /* AMATRIX_OPENCL_OPENCL_SHIM_H */
