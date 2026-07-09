/*
 * amatrix.opencl runtime loader for OpenCL and CLBlast.
 *
 * amatrix.opencl links against NEITHER OpenCL nor CLBlast. Every symbol from
 * both libraries is resolved at run time via dlopen()/dlsym() (POSIX) or
 * LoadLibrary()/GetProcAddress() (Windows) into the function-pointer tables
 * declared below. This makes the compiled shared object completely
 * self-contained: it can be built and loaded on any platform, with or without
 * a GPU, without the OpenCL ICD loader or CLBlast being present. Nothing
 * enumerates devices, opens a context, or maps a vendor ICD until the gated,
 * lazy probe (amatrix_opencl_try_init) is explicitly invoked.
 *
 * The spelling macros at the bottom map the ordinary OpenCL/CLBlast call
 * spellings used throughout opencl_bridge.c (e.g. clGetDeviceIDs, CLBlastSgemm)
 * onto the resolved function pointers, so the bridge source needs no call-site
 * edits. Define AMATRIX_CL_LOADER_NO_MACROS before including this header (the
 * loader implementation does) to suppress those macros.
 */
#ifndef AMATRIX_CL_LOADER_H
#define AMATRIX_CL_LOADER_H

#ifndef CL_TARGET_OPENCL_VERSION
#define CL_TARGET_OPENCL_VERSION 120
#endif

#include <stddef.h>
#include <CL/cl.h>
#include "clblast_c.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ---------------------------------------------------------------------------
 * OpenCL C API function-pointer typedefs (OpenCL 1.2 core).
 * ------------------------------------------------------------------------- */
typedef cl_int (CL_API_CALL *amatrix_pfn_clGetPlatformIDs)(
  cl_uint, cl_platform_id *, cl_uint *);
typedef cl_int (CL_API_CALL *amatrix_pfn_clGetDeviceIDs)(
  cl_platform_id, cl_device_type, cl_uint, cl_device_id *, cl_uint *);
typedef cl_int (CL_API_CALL *amatrix_pfn_clGetDeviceInfo)(
  cl_device_id, cl_device_info, size_t, void *, size_t *);
typedef cl_context (CL_API_CALL *amatrix_pfn_clCreateContext)(
  const cl_context_properties *, cl_uint, const cl_device_id *,
  void (CL_CALLBACK *)(const char *, const void *, size_t, void *),
  void *, cl_int *);
typedef cl_command_queue (CL_API_CALL *amatrix_pfn_clCreateCommandQueue)(
  cl_context, cl_device_id, cl_command_queue_properties, cl_int *);
typedef cl_mem (CL_API_CALL *amatrix_pfn_clCreateBuffer)(
  cl_context, cl_mem_flags, size_t, void *, cl_int *);
typedef cl_int (CL_API_CALL *amatrix_pfn_clEnqueueReadBuffer)(
  cl_command_queue, cl_mem, cl_bool, size_t, size_t, void *,
  cl_uint, const cl_event *, cl_event *);
typedef cl_int (CL_API_CALL *amatrix_pfn_clEnqueueWriteBuffer)(
  cl_command_queue, cl_mem, cl_bool, size_t, size_t, const void *,
  cl_uint, const cl_event *, cl_event *);
typedef cl_int (CL_API_CALL *amatrix_pfn_clEnqueueCopyBuffer)(
  cl_command_queue, cl_mem, cl_mem, size_t, size_t, size_t,
  cl_uint, const cl_event *, cl_event *);
typedef cl_int (CL_API_CALL *amatrix_pfn_clEnqueueReadBufferRect)(
  cl_command_queue, cl_mem, cl_bool, const size_t *, const size_t *,
  const size_t *, size_t, size_t, size_t, size_t, void *,
  cl_uint, const cl_event *, cl_event *);
typedef cl_int (CL_API_CALL *amatrix_pfn_clEnqueueWriteBufferRect)(
  cl_command_queue, cl_mem, cl_bool, const size_t *, const size_t *,
  const size_t *, size_t, size_t, size_t, size_t, const void *,
  cl_uint, const cl_event *, cl_event *);
typedef cl_program (CL_API_CALL *amatrix_pfn_clCreateProgramWithSource)(
  cl_context, cl_uint, const char **, const size_t *, cl_int *);
typedef cl_int (CL_API_CALL *amatrix_pfn_clBuildProgram)(
  cl_program, cl_uint, const cl_device_id *, const char *,
  void (CL_CALLBACK *)(cl_program, void *), void *);
typedef cl_kernel (CL_API_CALL *amatrix_pfn_clCreateKernel)(
  cl_program, const char *, cl_int *);
typedef cl_int (CL_API_CALL *amatrix_pfn_clSetKernelArg)(
  cl_kernel, cl_uint, size_t, const void *);
typedef cl_int (CL_API_CALL *amatrix_pfn_clEnqueueNDRangeKernel)(
  cl_command_queue, cl_kernel, cl_uint, const size_t *, const size_t *,
  const size_t *, cl_uint, const cl_event *, cl_event *);
typedef cl_int (CL_API_CALL *amatrix_pfn_clWaitForEvents)(
  cl_uint, const cl_event *);
typedef cl_int (CL_API_CALL *amatrix_pfn_clReleaseMemObject)(cl_mem);
typedef cl_int (CL_API_CALL *amatrix_pfn_clReleaseKernel)(cl_kernel);
typedef cl_int (CL_API_CALL *amatrix_pfn_clReleaseProgram)(cl_program);
typedef cl_int (CL_API_CALL *amatrix_pfn_clReleaseCommandQueue)(cl_command_queue);
typedef cl_int (CL_API_CALL *amatrix_pfn_clReleaseContext)(cl_context);
typedef cl_int (CL_API_CALL *amatrix_pfn_clReleaseEvent)(cl_event);

/* ---------------------------------------------------------------------------
 * CLBlast function-pointer typedefs (single-precision routines we use).
 * ------------------------------------------------------------------------- */
typedef CLBlastStatusCode (*amatrix_pfn_CLBlastSgemm)(
  CLBlastLayout, CLBlastTranspose, CLBlastTranspose, size_t, size_t, size_t,
  float, cl_mem, size_t, size_t, cl_mem, size_t, size_t, float, cl_mem,
  size_t, size_t, cl_command_queue *, cl_event *);
typedef CLBlastStatusCode (*amatrix_pfn_CLBlastSsyrk)(
  CLBlastLayout, CLBlastTriangle, CLBlastTranspose, size_t, size_t,
  float, cl_mem, size_t, size_t, float, cl_mem, size_t, size_t,
  cl_command_queue *, cl_event *);
typedef CLBlastStatusCode (*amatrix_pfn_CLBlastStrsm)(
  CLBlastLayout, CLBlastSide, CLBlastTriangle, CLBlastTranspose,
  CLBlastDiagonal, size_t, size_t, float, cl_mem, size_t, size_t, cl_mem,
  size_t, size_t, cl_command_queue *, cl_event *);
typedef CLBlastStatusCode (*amatrix_pfn_CLBlastShad)(
  size_t, float, cl_mem, size_t, size_t, cl_mem, size_t, size_t, float,
  cl_mem, size_t, size_t, cl_command_queue *, cl_event *);

/* ---------------------------------------------------------------------------
 * Resolved function pointers (defined in amatrix_cl_loader.c). They are NULL
 * until the corresponding library has been loaded successfully.
 * ------------------------------------------------------------------------- */
extern amatrix_pfn_clGetPlatformIDs         amatrix_p_clGetPlatformIDs;
extern amatrix_pfn_clGetDeviceIDs           amatrix_p_clGetDeviceIDs;
extern amatrix_pfn_clGetDeviceInfo          amatrix_p_clGetDeviceInfo;
extern amatrix_pfn_clCreateContext          amatrix_p_clCreateContext;
extern amatrix_pfn_clCreateCommandQueue     amatrix_p_clCreateCommandQueue;
extern amatrix_pfn_clCreateBuffer           amatrix_p_clCreateBuffer;
extern amatrix_pfn_clEnqueueReadBuffer      amatrix_p_clEnqueueReadBuffer;
extern amatrix_pfn_clEnqueueWriteBuffer     amatrix_p_clEnqueueWriteBuffer;
extern amatrix_pfn_clEnqueueCopyBuffer      amatrix_p_clEnqueueCopyBuffer;
extern amatrix_pfn_clEnqueueReadBufferRect  amatrix_p_clEnqueueReadBufferRect;
extern amatrix_pfn_clEnqueueWriteBufferRect amatrix_p_clEnqueueWriteBufferRect;
extern amatrix_pfn_clCreateProgramWithSource amatrix_p_clCreateProgramWithSource;
extern amatrix_pfn_clBuildProgram           amatrix_p_clBuildProgram;
extern amatrix_pfn_clCreateKernel           amatrix_p_clCreateKernel;
extern amatrix_pfn_clSetKernelArg           amatrix_p_clSetKernelArg;
extern amatrix_pfn_clEnqueueNDRangeKernel   amatrix_p_clEnqueueNDRangeKernel;
extern amatrix_pfn_clWaitForEvents          amatrix_p_clWaitForEvents;
extern amatrix_pfn_clReleaseMemObject       amatrix_p_clReleaseMemObject;
extern amatrix_pfn_clReleaseKernel          amatrix_p_clReleaseKernel;
extern amatrix_pfn_clReleaseProgram         amatrix_p_clReleaseProgram;
extern amatrix_pfn_clReleaseCommandQueue    amatrix_p_clReleaseCommandQueue;
extern amatrix_pfn_clReleaseContext         amatrix_p_clReleaseContext;
extern amatrix_pfn_clReleaseEvent           amatrix_p_clReleaseEvent;

extern amatrix_pfn_CLBlastSgemm             amatrix_p_CLBlastSgemm;
extern amatrix_pfn_CLBlastSsyrk             amatrix_p_CLBlastSsyrk;
extern amatrix_pfn_CLBlastStrsm             amatrix_p_CLBlastStrsm;
extern amatrix_pfn_CLBlastShad              amatrix_p_CLBlastShad;

/* ---------------------------------------------------------------------------
 * Loader API. All functions are idempotent and never throw, error, or abort.
 * ------------------------------------------------------------------------- */

/* Load the OpenCL ICD loader and resolve every CL_* pointer above. Returns 1
 * on success (all symbols resolved), 0 on any failure. On failure a reason is
 * available via amatrix_cl_opencl_reason(). */
int amatrix_cl_load_opencl(void);

/* Load CLBlast and resolve every CLBlast* pointer above. Returns 1 on success,
 * 0 on failure; failure reason via amatrix_cl_clblast_reason(). */
int amatrix_cl_load_clblast(void);

int amatrix_cl_opencl_loaded(void);
int amatrix_cl_clblast_loaded(void);
const char *amatrix_cl_opencl_reason(void);
const char *amatrix_cl_clblast_reason(void);

/* Register an additional directory (e.g. tools::R_user_dir("amatrix.opencl"))
 * to search for the CLBlast shared library. Pass NULL or "" to clear. */
void amatrix_cl_set_clblast_dir(const char *dir);

#ifdef __cplusplus
}
#endif

/* ---------------------------------------------------------------------------
 * Spelling macros: rewrite ordinary OpenCL/CLBlast call spellings to the
 * resolved pointers. Included AFTER the real prototypes from <CL/cl.h> and
 * "clblast_c.h" above, so those declarations are unaffected.
 * ------------------------------------------------------------------------- */
#ifndef AMATRIX_CL_LOADER_NO_MACROS
#define clGetPlatformIDs          amatrix_p_clGetPlatformIDs
#define clGetDeviceIDs            amatrix_p_clGetDeviceIDs
#define clGetDeviceInfo           amatrix_p_clGetDeviceInfo
#define clCreateContext           amatrix_p_clCreateContext
#define clCreateCommandQueue      amatrix_p_clCreateCommandQueue
#define clCreateBuffer            amatrix_p_clCreateBuffer
#define clEnqueueReadBuffer       amatrix_p_clEnqueueReadBuffer
#define clEnqueueWriteBuffer      amatrix_p_clEnqueueWriteBuffer
#define clEnqueueCopyBuffer       amatrix_p_clEnqueueCopyBuffer
#define clEnqueueReadBufferRect   amatrix_p_clEnqueueReadBufferRect
#define clEnqueueWriteBufferRect  amatrix_p_clEnqueueWriteBufferRect
#define clCreateProgramWithSource amatrix_p_clCreateProgramWithSource
#define clBuildProgram            amatrix_p_clBuildProgram
#define clCreateKernel            amatrix_p_clCreateKernel
#define clSetKernelArg            amatrix_p_clSetKernelArg
#define clEnqueueNDRangeKernel    amatrix_p_clEnqueueNDRangeKernel
#define clWaitForEvents           amatrix_p_clWaitForEvents
#define clReleaseMemObject        amatrix_p_clReleaseMemObject
#define clReleaseKernel           amatrix_p_clReleaseKernel
#define clReleaseProgram          amatrix_p_clReleaseProgram
#define clReleaseCommandQueue     amatrix_p_clReleaseCommandQueue
#define clReleaseContext          amatrix_p_clReleaseContext
#define clReleaseEvent            amatrix_p_clReleaseEvent
#define CLBlastSgemm              amatrix_p_CLBlastSgemm
#define CLBlastSsyrk              amatrix_p_CLBlastSsyrk
#define CLBlastStrsm              amatrix_p_CLBlastStrsm
#define CLBlastShad               amatrix_p_CLBlastShad
#endif /* AMATRIX_CL_LOADER_NO_MACROS */

#endif /* AMATRIX_CL_LOADER_H */
