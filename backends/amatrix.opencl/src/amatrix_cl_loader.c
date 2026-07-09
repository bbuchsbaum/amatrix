/*
 * amatrix.opencl runtime loader implementation.
 *
 * Resolves the OpenCL ICD loader and CLBlast entirely at run time. No OpenCL
 * or CLBlast symbol is referenced at link time, so the compiled object has no
 * external GPU dependencies. Every entry point here is written in plain C,
 * returns a status (never throws or aborts), and is safe to call on a machine
 * with no GPU, no ICD loader, and no CLBlast installed.
 */
#define AMATRIX_CL_LOADER_NO_MACROS
#include "amatrix_cl_loader.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
typedef HMODULE amatrix_dl_handle;
#else
#  include <dlfcn.h>
typedef void *amatrix_dl_handle;
#endif

/* ---- pointer storage ---------------------------------------------------- */
amatrix_pfn_clGetPlatformIDs         amatrix_p_clGetPlatformIDs = NULL;
amatrix_pfn_clGetDeviceIDs           amatrix_p_clGetDeviceIDs = NULL;
amatrix_pfn_clGetDeviceInfo          amatrix_p_clGetDeviceInfo = NULL;
amatrix_pfn_clCreateContext          amatrix_p_clCreateContext = NULL;
amatrix_pfn_clCreateCommandQueue     amatrix_p_clCreateCommandQueue = NULL;
amatrix_pfn_clCreateBuffer           amatrix_p_clCreateBuffer = NULL;
amatrix_pfn_clEnqueueReadBuffer      amatrix_p_clEnqueueReadBuffer = NULL;
amatrix_pfn_clEnqueueWriteBuffer     amatrix_p_clEnqueueWriteBuffer = NULL;
amatrix_pfn_clEnqueueCopyBuffer      amatrix_p_clEnqueueCopyBuffer = NULL;
amatrix_pfn_clEnqueueReadBufferRect  amatrix_p_clEnqueueReadBufferRect = NULL;
amatrix_pfn_clEnqueueWriteBufferRect amatrix_p_clEnqueueWriteBufferRect = NULL;
amatrix_pfn_clCreateProgramWithSource amatrix_p_clCreateProgramWithSource = NULL;
amatrix_pfn_clBuildProgram           amatrix_p_clBuildProgram = NULL;
amatrix_pfn_clCreateKernel           amatrix_p_clCreateKernel = NULL;
amatrix_pfn_clSetKernelArg           amatrix_p_clSetKernelArg = NULL;
amatrix_pfn_clEnqueueNDRangeKernel   amatrix_p_clEnqueueNDRangeKernel = NULL;
amatrix_pfn_clWaitForEvents          amatrix_p_clWaitForEvents = NULL;
amatrix_pfn_clReleaseMemObject       amatrix_p_clReleaseMemObject = NULL;
amatrix_pfn_clReleaseKernel          amatrix_p_clReleaseKernel = NULL;
amatrix_pfn_clReleaseProgram         amatrix_p_clReleaseProgram = NULL;
amatrix_pfn_clReleaseCommandQueue    amatrix_p_clReleaseCommandQueue = NULL;
amatrix_pfn_clReleaseContext         amatrix_p_clReleaseContext = NULL;
amatrix_pfn_clReleaseEvent           amatrix_p_clReleaseEvent = NULL;

amatrix_pfn_CLBlastSgemm             amatrix_p_CLBlastSgemm = NULL;
amatrix_pfn_CLBlastSsyrk             amatrix_p_CLBlastSsyrk = NULL;
amatrix_pfn_CLBlastStrsm             amatrix_p_CLBlastStrsm = NULL;
amatrix_pfn_CLBlastShad              amatrix_p_CLBlastShad = NULL;

/* ---- loader state ------------------------------------------------------- */
static int g_opencl_loaded = 0;
static int g_clblast_loaded = 0;
static amatrix_dl_handle g_opencl_handle = NULL;
static amatrix_dl_handle g_clblast_handle = NULL;
static char g_opencl_reason[256] = "OpenCL runtime has not been probed";
static char g_clblast_reason[256] = "CLBlast has not been probed";
static char g_clblast_dir[1024] = "";

static void amatrix_cl__set_reason(char *dst, const char *msg) {
  strncpy(dst, msg, 255);
  dst[255] = '\0';
}

/* ---- thin dynamic-linker shims ----------------------------------------- */
static amatrix_dl_handle amatrix_dl_open(const char *path) {
#ifdef _WIN32
  return LoadLibraryA(path);
#else
  return dlopen(path, RTLD_LAZY | RTLD_LOCAL);
#endif
}

static void *amatrix_dl_sym(amatrix_dl_handle handle, const char *name) {
#ifdef _WIN32
  /* GetProcAddress returns FARPROC; route through a void* for a portable
   * object-to-function-pointer conversion. */
  FARPROC p = GetProcAddress(handle, name);
  void *out = NULL;
  memcpy(&out, &p, sizeof(out) < sizeof(p) ? sizeof(out) : sizeof(p));
  return out;
#else
  return dlsym(handle, name);
#endif
}

static amatrix_dl_handle amatrix_dl_open_first(const char *const *paths, size_t n) {
  size_t i;
  for (i = 0; i < n; ++i) {
    amatrix_dl_handle h;
    if (paths[i] == NULL || paths[i][0] == '\0') {
      continue;
    }
    h = amatrix_dl_open(paths[i]);
    if (h != NULL) {
      return h;
    }
  }
  return NULL;
}

/* ---- OpenCL --------------------------------------------------------------
 * Resolve the ICD loader and every core symbol we use. All-or-nothing: if any
 * symbol is missing, the load fails and no partial pointer table is exposed.
 */
int amatrix_cl_load_opencl(void) {
  amatrix_dl_handle h = NULL;

  if (g_opencl_loaded) {
    return 1;
  }

#if defined(_WIN32)
  {
    static const char *const paths[] = {"OpenCL.dll"};
    h = amatrix_dl_open_first(paths, sizeof(paths) / sizeof(paths[0]));
  }
#elif defined(__APPLE__)
  {
    static const char *const paths[] = {
      "/System/Library/Frameworks/OpenCL.framework/OpenCL",
      "/System/Library/Frameworks/OpenCL.framework/Versions/A/OpenCL",
      "libOpenCL.dylib"
    };
    h = amatrix_dl_open_first(paths, sizeof(paths) / sizeof(paths[0]));
  }
#else
  {
    static const char *const paths[] = {"libOpenCL.so.1", "libOpenCL.so"};
    h = amatrix_dl_open_first(paths, sizeof(paths) / sizeof(paths[0]));
  }
#endif

  if (h == NULL) {
    amatrix_cl__set_reason(
      g_opencl_reason,
      "OpenCL runtime library not found (install GPU drivers that provide the "
      "OpenCL ICD loader)");
    return 0;
  }

#define AMATRIX_RESOLVE_CL(name)                                             \
  do {                                                                       \
    void *sym_ = amatrix_dl_sym(h, #name);                                   \
    if (sym_ == NULL) {                                                      \
      amatrix_cl__set_reason(                                               \
        g_opencl_reason,                                                     \
        "OpenCL runtime is missing required symbol " #name);                \
      goto opencl_fail;                                                      \
    }                                                                        \
    amatrix_p_##name = (amatrix_pfn_##name)sym_;                            \
  } while (0)

  AMATRIX_RESOLVE_CL(clGetPlatformIDs);
  AMATRIX_RESOLVE_CL(clGetDeviceIDs);
  AMATRIX_RESOLVE_CL(clGetDeviceInfo);
  AMATRIX_RESOLVE_CL(clCreateContext);
  AMATRIX_RESOLVE_CL(clCreateCommandQueue);
  AMATRIX_RESOLVE_CL(clCreateBuffer);
  AMATRIX_RESOLVE_CL(clEnqueueReadBuffer);
  AMATRIX_RESOLVE_CL(clEnqueueWriteBuffer);
  AMATRIX_RESOLVE_CL(clEnqueueCopyBuffer);
  AMATRIX_RESOLVE_CL(clEnqueueReadBufferRect);
  AMATRIX_RESOLVE_CL(clEnqueueWriteBufferRect);
  AMATRIX_RESOLVE_CL(clCreateProgramWithSource);
  AMATRIX_RESOLVE_CL(clBuildProgram);
  AMATRIX_RESOLVE_CL(clCreateKernel);
  AMATRIX_RESOLVE_CL(clSetKernelArg);
  AMATRIX_RESOLVE_CL(clEnqueueNDRangeKernel);
  AMATRIX_RESOLVE_CL(clWaitForEvents);
  AMATRIX_RESOLVE_CL(clReleaseMemObject);
  AMATRIX_RESOLVE_CL(clReleaseKernel);
  AMATRIX_RESOLVE_CL(clReleaseProgram);
  AMATRIX_RESOLVE_CL(clReleaseCommandQueue);
  AMATRIX_RESOLVE_CL(clReleaseContext);
  AMATRIX_RESOLVE_CL(clReleaseEvent);

#undef AMATRIX_RESOLVE_CL

  g_opencl_handle = h;
  g_opencl_loaded = 1;
  amatrix_cl__set_reason(g_opencl_reason, "OpenCL runtime loaded");
  return 1;

opencl_fail:
  /* Leave the ICD loader mapped but report failure; callers treat OpenCL as
   * unavailable. Pointers are not published because g_opencl_loaded stays 0. */
  return 0;
}

/* ---- CLBlast ------------------------------------------------------------- */
void amatrix_cl_set_clblast_dir(const char *dir) {
  if (dir == NULL) {
    g_clblast_dir[0] = '\0';
    return;
  }
  strncpy(g_clblast_dir, dir, sizeof(g_clblast_dir) - 1);
  g_clblast_dir[sizeof(g_clblast_dir) - 1] = '\0';
}

/* Base filenames to try inside amatrix_cl_set_clblast_dir() and on the loader
 * search path, per platform. */
static const char *const *amatrix_clblast_basenames(size_t *n) {
#if defined(_WIN32)
  static const char *const names[] = {"clblast.dll"};
#elif defined(__APPLE__)
  static const char *const names[] = {"libclblast.dylib", "libclblast.1.dylib"};
#else
  static const char *const names[] = {"libclblast.so", "libclblast.so.1"};
#endif
  *n = sizeof(names) / sizeof(names[0]);
  return names;
}

static amatrix_dl_handle amatrix_clblast_open(void) {
  amatrix_dl_handle h = NULL;
  const char *env = getenv("AMATRIX_CLBLAST_LIB");
  size_t nbase = 0;
  const char *const *base = amatrix_clblast_basenames(&nbase);
  size_t i;

  /* 1. Explicit full path from the environment. */
  if (env != NULL && env[0] != '\0') {
    h = amatrix_dl_open(env);
    if (h != NULL) {
      return h;
    }
  }

  /* 2. A directory registered from R (tools::R_user_dir). */
  if (g_clblast_dir[0] != '\0') {
    for (i = 0; i < nbase; ++i) {
      char path[1200];
#if defined(_WIN32)
      snprintf(path, sizeof(path), "%s\\%s", g_clblast_dir, base[i]);
#else
      snprintf(path, sizeof(path), "%s/%s", g_clblast_dir, base[i]);
#endif
      h = amatrix_dl_open(path);
      if (h != NULL) {
        return h;
      }
    }
  }

  /* 3. Bare names, resolved via the platform loader search path. */
  h = amatrix_dl_open_first(base, nbase);
  if (h != NULL) {
    return h;
  }

  /* 4. Common install locations not always on the default search path. */
#if defined(__APPLE__)
  {
    static const char *const extra[] = {
      "/opt/homebrew/lib/libclblast.dylib",
      "/opt/homebrew/opt/clblast/lib/libclblast.dylib",
      "/usr/local/lib/libclblast.dylib",
      "/usr/local/opt/clblast/lib/libclblast.dylib"
    };
    h = amatrix_dl_open_first(extra, sizeof(extra) / sizeof(extra[0]));
    if (h != NULL) {
      return h;
    }
  }
#elif !defined(_WIN32)
  {
    static const char *const extra[] = {
      "/usr/local/lib/libclblast.so",
      "/usr/lib/libclblast.so",
      "/usr/lib/x86_64-linux-gnu/libclblast.so"
    };
    h = amatrix_dl_open_first(extra, sizeof(extra) / sizeof(extra[0]));
    if (h != NULL) {
      return h;
    }
  }
#endif

  return NULL;
}

int amatrix_cl_load_clblast(void) {
  amatrix_dl_handle h = NULL;

  if (g_clblast_loaded) {
    return 1;
  }

  h = amatrix_clblast_open();
  if (h == NULL) {
    amatrix_cl__set_reason(
      g_clblast_reason,
      "CLBlast library not found; run amatrix_install_clblast() or set "
      "AMATRIX_CLBLAST_LIB to enable GPU BLAS/linear algebra");
    return 0;
  }

#define AMATRIX_RESOLVE_CLBLAST(name)                                        \
  do {                                                                       \
    void *sym_ = amatrix_dl_sym(h, #name);                                   \
    if (sym_ == NULL) {                                                      \
      amatrix_cl__set_reason(                                               \
        g_clblast_reason, "CLBlast library is missing required symbol " #name); \
      goto clblast_fail;                                                     \
    }                                                                        \
    amatrix_p_##name = (amatrix_pfn_##name)sym_;                            \
  } while (0)

  AMATRIX_RESOLVE_CLBLAST(CLBlastSgemm);
  AMATRIX_RESOLVE_CLBLAST(CLBlastSsyrk);
  AMATRIX_RESOLVE_CLBLAST(CLBlastStrsm);
  AMATRIX_RESOLVE_CLBLAST(CLBlastShad);

#undef AMATRIX_RESOLVE_CLBLAST

  g_clblast_handle = h;
  g_clblast_loaded = 1;
  amatrix_cl__set_reason(g_clblast_reason, "CLBlast loaded");
  return 1;

clblast_fail:
  amatrix_p_CLBlastSgemm = NULL;
  amatrix_p_CLBlastSsyrk = NULL;
  amatrix_p_CLBlastStrsm = NULL;
  amatrix_p_CLBlastShad = NULL;
  return 0;
}

int amatrix_cl_opencl_loaded(void) { return g_opencl_loaded; }
int amatrix_cl_clblast_loaded(void) { return g_clblast_loaded; }
const char *amatrix_cl_opencl_reason(void) { return g_opencl_reason; }
const char *amatrix_cl_clblast_reason(void) { return g_clblast_reason; }
