% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multisite.bin.R
\name{multisite.bin}
\alias{multisite.bin}
\title{Power for a multisite trial with a binary outcome}
\usage{
multisite.bin(
  m = NULL,
  alloc.ratio = 1,
  J = NULL,
  pc = NULL,
  pt = NULL,
  sigma.u = NULL,
  alpha = 0.05,
  power = NULL,
  sides = 2,
  v = FALSE
)
}
\arguments{
\item{m}{The total number of subjects in condition 1 + condition 2.}

\item{alloc.ratio}{The allocation ratio of condition 1/condition 2 within site; defaults to 1.}

\item{J}{The total number of sites.}

\item{pc}{The probability of the outcome in the control condition.}

\item{pt}{The probability of the outcome in the treatment condition.}

\item{sigma.u}{Standard deviation of the treatment effect across sites.}

\item{alpha}{The significance level or type 1 error rate; defaults to 0.05.}

\item{power}{The specified level of power.}

\item{sides}{Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.}

\item{v}{Either TRUE for verbose output or FALSE (default) to output computed argument only.}
}
\value{
A list of the arguments (including the computed one).
}
\description{
Performs power and sample size calculations for a multisite trial with a
binary outcome. Can solve for power, J, m or alpha.
}
\details{
In a multisite trial design, participants are randomized to conditions
within site. Consider using ms.varexplore to select plausible values
for sigma.u.
}
\examples{
multisite.bin(m = 30, J = 25, pc = 0.1, pt = 0.2, sigma.u = 0.4, power = NULL)
multisite.bin(m = 30, J = NULL, pc = 0.1, pt = 0.2, sigma.u = 0.4, power = 0.9)
}
