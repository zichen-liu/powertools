% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/relrisk.R
\name{relrisk}
\alias{relrisk}
\title{Power calculation for relative risk}
\usage{
relrisk(
  n1 = NULL,
  n.ratio = 1,
  p1 = NULL,
  p2 = NULL,
  RR0 = 1,
  alpha = 0.05,
  power = NULL,
  sides = 2,
  v = FALSE
)
}
\arguments{
\item{n1}{The sample size for group 1.}

\item{n.ratio}{The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).}

\item{p1}{The proportion in group 1.}

\item{p2}{The proportion in group 2.}

\item{RR0}{The relative risk under the null (p2/p1); defaults to 1.}

\item{alpha}{The significance level (type 1 error rate); defaults to 0.05.}

\item{power}{The specified level of power.}

\item{sides}{Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.}

\item{v}{Either TRUE for verbose output or FALSE (default) to output computed argument only.}
}
\value{
A list of the arguments (including the computed one).
}
\description{
Performs power and sample size calculations for a test of relative risk, p2/p1.
Can solve for power, N1, n.ratio or alpha.
}
\examples{
relrisk(n1 = NULL, n.ratio = 1/6, p1 = 0.1, p2 = 0.1 * 2, power = 0.8)
}
