% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check.many.R
\name{check.many}
\alias{check.many}
\title{Internal-use function for quality checking sets of parameters}
\usage{
check.many(paramlist = NULL, type = c("oneof"))
}
\arguments{
\item{paramlist}{The list of parameters being checked.}

\item{type}{The expected type of parameter: currently only supports "oneof".}
}
\value{
If the check passes, returns nothing. If the check does not pass, throw an error.
}
\description{
Internal-use function for quality checking sets of parameters
}
\examples{
N <- 10
power <- NULL
check.many(list(N, power), "oneof")
}
\keyword{internal}
