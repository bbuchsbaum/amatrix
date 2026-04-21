.amatrix_state <- new.env(parent = emptyenv())
.amatrix_state$backends <- new.env(parent = emptyenv())
.amatrix_state$default_policy <- "auto"
.amatrix_state$default_precision <- "strict"
.amatrix_state$residency <- new.env(parent = emptyenv())
.amatrix_state$model_cache <- new.env(parent = emptyenv())
.amatrix_state$resident_counter <- 0L
.amatrix_state$object_counter <- 0L
.amatrix_state$session_id <- ""

.onLoad <- function(libname, pkgname) {
  .amatrix_state$session_id <- paste0(
    format(Sys.time(), "%Y%m%d%H%M%OS6"), "-",
    as.hexmode(sample.int(2^31 - 1L, 1L))
  )
  ns <- asNamespace(pkgname)
  registerS3method("as.matrix", "KronMatrix", get("as.matrix.KronMatrix", envir = ns), envir = ns)
  registerS3method("as.matrix", "resident_handle", get("as.matrix.resident_handle", envir = ns), envir = ns)
  .amatrix_cache_init()
  amatrix_register_backend("cpu", .amatrix_cpu_backend(), overwrite = TRUE)
}
