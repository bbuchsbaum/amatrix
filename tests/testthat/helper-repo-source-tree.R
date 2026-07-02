# Locate the amatrix *source* tree (not an installed package directory) so that
# subprocess regression tests can `pkgload::load_all()` it.
#
# The naive "has a DESCRIPTION" check also matches the installed package dir
# returned by getNamespaceInfo(asNamespace("amatrix"), "path") under R CMD
# check. load_all()-ing that directory fails hard in the child process (e.g.
# "getDLLRegisteredRoutines.DLLInfo must specify DLL") because it has no R
# source and no compiled sources to build. A genuine source checkout has all
# three of: a DESCRIPTION, an R/ directory containing .R files, and a src/
# directory (amatrix ships compiled C). Requiring all three means that in an
# installed-package context this returns NULL and the callers' skip_if() guards
# fire cleanly instead of launching a doomed subprocess.
.amatrix_source_tree_dir <- function() {
  candidates <- unique(c(
    tryCatch(getNamespaceInfo(asNamespace("amatrix"), "path"), error = function(e) NULL),
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", "..")
  ))
  candidates <- Filter(Negate(is.null), candidates)

  is_source_tree <- function(dir) {
    if (!file.exists(file.path(dir, "DESCRIPTION"))) {
      return(FALSE)
    }
    r_dir <- file.path(dir, "R")
    has_r_sources <- dir.exists(r_dir) &&
      length(list.files(r_dir, pattern = "\\.[Rr]$")) > 0L
    has_src <- dir.exists(file.path(dir, "src"))
    has_r_sources && has_src
  }

  matches <- Filter(is_source_tree, candidates)
  if (length(matches) == 0L) {
    return(NULL)
  }
  normalizePath(matches[[1L]], winslash = "/", mustWork = TRUE)
}
