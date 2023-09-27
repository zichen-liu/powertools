#' Power for test of average treatment effect

#' @param m The number of subjects per site.
#' @param J The number of sites.
#' @param r The allocation ratio to intervention and control at each site; defaults to 1.
#' @param gamma The average treatment effect under the alternative.
#' @param rho0 The proportion of total variance of the outcome attributable to variation in site-level means.
#' @param rho1 The proportion of total variance of the outcome attributable to variation in the treatment effect across sites.
#' @param ssq.Y The total variance of the outcome variable Y.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.

#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples # Example 14.2
#' pss.multi.ate(m = 20, J = 10, gamma = 3, rho0 = 0.1, rho1 = 0, ssq.Y = 40)
#' pss.multi.ate(m = 20, J = 10, gamma = 3, rho0 = 0.095, rho1 = 0.048, ssq.Y = 48)
#' # Example 14.3
#' pss.multi.ate(m = 20, J = 10, r = 1.5, gamma = 3, rho0 = 0.095, rho1 = 0.048, ssq.Y = 48)

pss.multi.ate <- function(m, J, r = 1, gamma, rho0, rho1, ssq.Y, alpha = 0.05) {
  N <- m * J
  d <- gamma / sqrt(ssq.Y)
  c <- (1 + r)^2 / r
  ncp <- d / sqrt(c * (1 - rho0 + (4 * m / c - 1) * rho1) / N)
  df <- J - 1
  crit <- stats::qt(1 - alpha / 2, df)
  power <- 1 - stats::pt(crit, df, ncp)

  # Generate output text
  METHOD <-"Power for test of average treatment effect"

  # Print output as a power.htest object depending on which inputs were given
  structure(list(m = m, J = J, r = r, gamma = gamma,
                 rho0 = rho0, rho1 = rho1, ssq.Y = ssq.Y,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
