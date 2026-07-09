/*
 * Minimal <CL/opencl.h> compatibility shim provided by amatrix.opencl.
 *
 * The vendored CLBlast C header (clblast_c.h) includes <CL/opencl.h> to obtain
 * the OpenCL C types (cl_mem, cl_command_queue, cl_event, ...). amatrix.opencl
 * resolves every OpenCL and CLBlast symbol at runtime via dlopen/LoadLibrary,
 * so it never uses the OpenCL/GL or extension APIs that the upstream umbrella
 * header pulls in. This shim provides only the core C API declarations, which
 * is all clblast_c.h and the amatrix bridge require, and avoids vendoring the
 * cl_gl.h / cl_ext.h extension headers.
 *
 * This file is authored by the amatrix project (MIT), not by The Khronos Group.
 * The core OpenCL headers it includes (CL/cl.h, CL/cl_platform.h,
 * CL/cl_version.h) are the unmodified Khronos headers (Apache-2.0); see
 * inst/COPYRIGHTS.
 */
#ifndef AMATRIX_CL_OPENCL_SHIM_H
#define AMATRIX_CL_OPENCL_SHIM_H

#ifndef CL_TARGET_OPENCL_VERSION
#define CL_TARGET_OPENCL_VERSION 120
#endif

#include <CL/cl.h>

#endif /* AMATRIX_CL_OPENCL_SHIM_H */
