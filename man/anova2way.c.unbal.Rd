% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anova2way.c.unbal.R
\name{anova2way.c.unbal}
\alias{anova2way.c.unbal}
\title{Power calculation for two-way unbalanced analysis of variance contrast test}
\usage{
anova2way.c.unbal(
  nmatrix = nmatrix,
  mmatrix = NULL,
  cvec = NULL,
  factor = c("a", "b"),
  sd = NULL,
  Rsq = 0,
  ncov = 0,
  alpha = 0.05,
  v = FALSE
)
}
\arguments{
\item{nmatrix}{A matrix of cell sample sizes (see example).}

\item{mmatrix}{A matrix of cell means (see example).}

\item{cvec}{A vector of contrast coefficients c(c1, c2, ...).}

\item{factor}{Either "a" (rows) or "b" (columns) depending on which factor the contrast test is being made on.}

\item{sd}{The estimated standard deviation within each cell.}

\item{Rsq}{The estimated R^2 for regressing the outcome on the covariates; defaults to 0.}

\item{ncov}{The number of covariates adjusted for in the model; defaults to 0.}

\item{alpha}{The significance level (type 1 error rate); defaults to 0.05.}

\item{v}{Either TRUE for verbose output or FALSE (default) to output computed argument only.}
}
\value{
A list of the arguments (including the computed one).
}
\description{
Calculates power for a test of a contrast between levels of a factor in a two-way
ANOVA with unbalanced data (that is, unequal cell sizes).
This function
only solves for power. For a two-way balanced ANOVA,
(equal cell sizes), anova2way.c.bal can also be used, and will solve for
quantities other than power. For a test of a contrast among cell means, see
anova2way.se.unbal.
}
\examples{
nmatrix <- matrix(c(30, 30, 30, 30, 30, 30), nrow = 2, byrow = TRUE)
mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
anova2way.c.unbal(nmatrix = nmatrix, mmatrix = mmatrix, cvec = c(1, 0, -1),
factor = "b", sd = 2, alpha = 0.05)
}
