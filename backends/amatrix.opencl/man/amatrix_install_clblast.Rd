% Generated manually; keep in sync with the roxygen block in R/install.R.
\name{amatrix_install_clblast}
\alias{amatrix_install_clblast}
\title{Download and install the CLBlast runtime library}
\usage{
amatrix_install_clblast(force = FALSE, quiet = FALSE, version = "1.7.0")
}
\arguments{
\item{force}{Logical. Reinstall even if CLBlast is already present in the
data directory, and consent to downloading without an interactive prompt.
Required in non-interactive sessions. Default \code{FALSE}.}

\item{quiet}{Logical. Suppress progress messages. Default \code{FALSE}.}

\item{version}{Character. CLBlast release version (git tag) to install.
Default \code{"1.7.0"}.}
}
\value{
Invisibly, the path to the installed CLBlast shared library, or
\code{NULL} if the user declined the download.
}
\description{
\code{amatrix.opencl} loads CLBlast (the OpenCL BLAS library that powers
GPU matrix multiply, \code{crossprod}, triangular solves, and the dense
factorizations) at run time; it is never linked at build time. On most
systems CLBlast is not preinstalled. This helper downloads the official
CLBlast release binary for the current platform from GitHub and places its
shared library in the package's user data directory
(\code{tools::R_user_dir("amatrix.opencl", "data")}), where the runtime
loader searches for it. It is most useful on Windows, where CLBlast is
rarely available otherwise.
}
\details{
The download requires consent: in an interactive session you are prompted
before anything is fetched; in a non-interactive session you must pass
\code{force = TRUE} to proceed. Nothing about this function runs at package
load or during \code{R CMD check}; the OpenCL device itself is only touched
later, by the gated probe (see \code{amatrix_opencl_enable_probe}).

No official CLBlast release binary is published for macOS arm64 (Apple
Silicon) or Linux aarch64. On those platforms the function stops with a
pointer to the system package manager (e.g. \code{brew install clblast});
the Homebrew library is discovered automatically. You can always bypass this
helper entirely by setting the \code{AMATRIX_CLBLAST_LIB} environment
variable to the full path of a CLBlast shared library.
}
\seealso{
\code{amatrix_opencl_enable_probe},
\code{amatrix_opencl_diagnostics}
}
\examples{
\dontrun{
# Interactive: prompts before downloading.
amatrix_install_clblast()

# Non-interactive (e.g. a setup script): explicit consent.
amatrix_install_clblast(force = TRUE, quiet = TRUE)
}
}
