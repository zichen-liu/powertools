% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prop.2samp.R
\name{prop.2samp}
\alias{prop.2samp}
\title{Power calculation for two-sample proportion test}
\usage{
prop.2samp(
  n1 = NULL,
  n.ratio = 1,
  p1 = NULL,
  p2 = NULL,
  margin = 0,
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

\item{margin}{The margin of noninferiority or superiority; defaults to 0.}

\item{alpha}{The significance level (type 1 error rate); defaults to 0.05.}

\item{power}{The specified level of power.}

\item{sides}{Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.}

\item{v}{Either TRUE for verbose output or FALSE (default) to output computed argument only.}
}
\value{
A list of the arguments (including the computed one).
}
\description{
Performs power and sample size calculations for two-sample tests of proportions
using normal approximation to the binomial. Can solve for power, n1, n.ratio
or alpha.
}
\details{
For a noninferiority or superiority by a margin test, a one-sided test should be used. See Crespi (2025)
for more guidance. For an equivalence test for two proportions, see the prop.test.equiv.
}
\examples{
prop.2samp(n1 = NULL, p1 = 0.6, p2 = 0.8, alpha = 0.025, power = 0.9, sides = 1)
prop.2samp(n1 = NULL, p1 = 0.25, p2 = 0.25, margin = 0.1, alpha = 0.025, power = 0.8, sides = 1)
}
