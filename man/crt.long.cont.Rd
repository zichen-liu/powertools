% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crt.long.cont.R
\name{crt.long.cont}
\alias{crt.long.cont}
\title{Power for test of treatment effect in longitudinal cluster randomized trial with baseline measurement}
\usage{
crt.long.cont(
  m = NULL,
  J1 = NULL,
  J.ratio = 1,
  delta = NULL,
  sd = 1,
  icc = 0,
  cac = 0,
  sac = 0,
  alpha = 0.05,
  power = NULL,
  sides = 2,
  v = FALSE
)
}
\arguments{
\item{m}{The number of subjects measured during each cluster-period.}

\item{J1}{The number of clusters in arm 1.}

\item{J.ratio}{The ratio J2/J1 between the number of clusters in the two arms; defaults to 1 (equal clusters per arm).}

\item{delta}{The difference between the intervention and control means under the alternative minus the difference under the null hypothesis.}

\item{sd}{The total standard deviation of the outcome variable; defaults to 1.}

\item{icc}{The within-cluster, within-period intraclass correlation coefficient; defaults to 0.}

\item{cac}{The cluster autocorrelation; defaults to 0.}

\item{sac}{The subject autocorrelation; defaults to 0.}

\item{alpha}{The significance level or type 1 error rate; defaults to 0.05.}

\item{power}{The specified level of power.}

\item{sides}{Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.}

\item{v}{Either TRUE for verbose output or FALSE to output computed argument only.}
}
\value{
A list of the arguments (including the computed one).
}
\description{
Power for test of treatment effect in longitudinal cluster randomized trial with baseline measurement
}
\examples{
crt.long.cont(m = 30, J1 = 8, delta = 0.3, icc = 0.05, cac = 0.4, sac = 0.5)
}
