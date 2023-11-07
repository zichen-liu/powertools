#' Power for test of average treatment effect

#' @param m The number of subjects in the control for each site.
#' @param m.ratio The allocation ratio per site
#' @param J The number of sites.
#' @param gamma The average treatment effect under the alternative.
#' @param rho0 The proportion of total variance of the outcome attributable to variation in site-level means.
#' @param rho1 The proportion of total variance of the outcome attributable to variation in the treatment effect across sites.
#' @param ssq.Y The total variance of the outcome variable Y.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @keywords internal
#'
#' @examples
#' pss.multisite.ate.unbal(m = 20, m.ratio = 1.5, J = 10, gamma = 3, rho0 = 0.0952381, rho1 = 0.04761905, ssq.Y = 48)

pss.multisite.ate.unbal <- function (m = NULL, m.ratio = NULL, J = NULL, gamma = NULL,
                                   rho0 = NULL, rho1 = NULL, ssq.Y = NULL, d = NULL,
                                   alpha = 0.05, sides = 2) {

  N <- m * J
  c <- (1 + m.ratio)^2 / m.ratio
  d <- gamma / sqrt(ssq.Y)
  ncp <- d / sqrt(c * (1 - rho0 + (4 * m / c - 1) * rho1) / N)
  df <- J - 1
  crit <- stats::qt(1 - alpha / sides, df)
  power <- 1 - stats::pt(crit, df, ncp)

  # Generate output text
  METHOD <-"Power for test of average treatment effect"

  # Print output as a power.htest object depending on which inputs were given
  structure(list(m = m, J = J, gamma = gamma,
                 rho0 = rho0, rho1 = rho1, ssq.Y = ssq.Y,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
