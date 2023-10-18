#' Power for test of average treatment effect
#'
#' @param m The number of subjects per site or the mean cluster size (if unequal number of participants per site).
#' @param m.sd The standard deviation of cluster sizes (provide if unequal number of participants per site); defaults to 0.
#' @param ICC The intraclass correlation between cluster sizes ((provide if unequal number of participants per site); defaults to 0.
#' @param J The number of sites.
#' @param delta The difference between the intervention and control means in the outcome variable.
#' @param var The total variance of the outcome variable; defaults to 1.
#' @param rho0 The proportion of total variance of the outcome attributable to variation in site-level means.
#' @param rho1 The proportion of total variance of the outcome attributable to variation in the treatment effect across sites.
#' @param rho.cov The proportion of variance in the outome attributable to an individual-level covariate; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' pss.multisite.ate(m = 20, J = 10, delta = 3, var = 40, rho0 = 0.1, rho1 = 0)
#' pss.multisite.ate(m = 20, J = 10, delta = 3, var = 48, rho0 = 0.095, rho1 = 0.048)
#' pss.multisite.ate(m = 20, m.sd = 5, ICC = 0.5, J = 10, delta = 3, var = 48, rho0 = 0.095, rho1 = 0.048)
#' pss.multisite.ate(m = 20, J = 10, delta = 3, var = 48, rho0 = 0.095, rho1 = 0.048, rho.cov = 0.5)

pss.multisite.ate <- function (m = NULL, m.mean = NULL, m.sd = 0, ICC = 0,
                               J = NULL, delta = NULL, var = 1,
                               rho0 = NULL, rho1 = NULL, rho.cov = 0,
                               alpha = 0.05, sides = 2) {

  # Calculate power
  N <- m * J
  df <- J - 1
  d <- delta / sqrt(var * (1 - rho.cov^2)) # Adjust for covariate
  # Adjust for relative efficiency
  RE <- pss.multisite.re(m.mean = m, m.sd = m.sd, rho = ICC)$re
  ncp <- d / sqrt(4 * (1 - rho0) / N) / RE
  crit <- stats::qt(1 - alpha / sides, df)
  power <- 1 - stats::pt(crit, df, ncp)

  # Generate output text
  METHOD <-"Power for test of average treatment effect"

  # Print output as a power.htest object depending on which inputs were given
  structure(list(m = m, J = J, delta = delta, var = var,
                 rho0 = rho0, rho1 = rho1, rho.cov = rho.cov,
                 alpha = alpha, power = power, RE = RE,
                 method = METHOD), class = "power.htest")

}

