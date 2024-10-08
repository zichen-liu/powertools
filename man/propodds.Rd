% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/propodds.R
\name{propodds}
\alias{propodds}
\title{Power calculations for ordinal categorical variable under proportional odds assumption}
\usage{
propodds(pC, OR, n1, n.ratio = 1, alpha = 0.05, power = NULL, v = FALSE)
}
\arguments{
\item{pC}{Vector of response probabilities in control group (group 1). Must sum to 1. Categories are ordered from best to worst.}

\item{OR}{Odds ratio when the alternative is true. Must be greater than 1.}

\item{n1}{Sample size for group 1 (control group).}

\item{n.ratio}{The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).}

\item{alpha}{The significance level (type 1 error rate); defaults to 0.05.}

\item{power}{The specified level of power; defaults to 0.8.}

\item{v}{Either TRUE for verbose output or FALSE (default) to output computed argument only.}
}
\value{
A list of the arguments (including the computed one).
}
\description{
Performs power and sample size calculation for a comparison of two groups on
an ordinal categorical response variable.
Assumes that response probabilities follow the proportional
odds assumption. Can solve for power, n1, n.ratio and alpha.
}
\details{
Whitehead J. (1993) Sample size calculations for ordered categorical data.
Statistics in Medicine, 12(24):2257–2271
}
\examples{
library(Hmisc)
pC <- c(0.2, 0.5, 0.2, 0.1)
propodds(pC = pC, OR = 2.5, n1 = 65, n.ratio = 1, alpha = 0.05)
}
