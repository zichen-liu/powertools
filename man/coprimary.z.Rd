% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coprimary.z.R
\name{coprimary.z}
\alias{coprimary.z}
\title{Power calculations for multiple co-primary continuous endpoints assuming known covariance matrix}
\usage{
coprimary.z(
  K,
  n1 = NULL,
  n.ratio = 1,
  delta = NULL,
  Sigma,
  sd,
  rho,
  alpha = 0.025,
  power = NULL,
  tol = .Machine$double.eps^0.25,
  v = FALSE
)
}
\arguments{
\item{K}{The number of endpoints.}

\item{n1}{The sample size for group 1.}

\item{n.ratio}{The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).}

\item{delta}{A vector of length K of the true mean differences mu1k - mu2k; must all be positive.}

\item{Sigma}{The covariance matrix of the K outcomes, of dimension K x K.}

\item{sd}{A vector of length K of the standard deviations of the K outcomes.}

\item{rho}{A vector of length 0.5*K*(K-1) of the correlations among the K outcomes.}

\item{alpha}{The significance level or type 1 error rate; defaults to 0.025. A one-sided test is assumed.}

\item{power}{The specified level of power.}

\item{tol}{The desired accuracy (convergence tolerance) for uniroot.}

\item{v}{Either TRUE for verbose output or FALSE to output computed argument only.}
}
\value{
A list of the arguments (including the computed one).
}
\description{
Power calculations for multiple co-primary continuous endpoints assuming known covariance matrix
}
\examples{
coprimary.z(K = 2, n1 = 100, delta = c(0.4, 0.5), sd = c(1, 1), rho = 0.3,
alpha = 0.025, power = NULL)
}
