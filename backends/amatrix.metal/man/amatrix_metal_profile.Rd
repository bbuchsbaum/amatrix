\name{amatrix_metal_profile}
\alias{amatrix_metal_profile}
\alias{amatrix_metal_profile_enable}
\alias{amatrix_metal_profile_reset}
\title{Profile Metal backend bridge timing}
\usage{
amatrix_metal_profile(reset = FALSE)

amatrix_metal_profile_enable(enabled = TRUE, reset = FALSE)

amatrix_metal_profile_reset()
}
\arguments{
\item{reset}{Logical scalar. For \code{amatrix_metal_profile()}, reset counters after reading them. For \code{amatrix_metal_profile_enable()}, reset counters before changing the enabled state.}
\item{enabled}{Logical scalar. Whether native Metal bridge profiling should collect timing counters.}
}
\value{
\code{amatrix_metal_profile()} returns a named numeric vector of timing counters in milliseconds plus call counts. The enable and reset helpers return invisibly.
}
\description{
Enable, reset, and read lightweight native timing counters for the Metal backend bridge. Counters cover sparse and dense uploads, sparse product command submission and synchronization, deferred pending waits, and host materialization. Profiling is disabled by default to avoid overhead in normal use.
}
\examples{
\dontrun{
amatrix_metal_enable_probe(register = TRUE)
amatrix_metal_profile_enable(TRUE, reset = TRUE)

# Run Metal-backed operations here.

amatrix_metal_profile()
amatrix_metal_profile_enable(FALSE)
}
}
