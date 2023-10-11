#' Power for test of average treatment effect

#' @param m The number of subjects per site.
#' @param m.ratio The allocation ratio to intervention and control at each site; defaults to 1.
#' @param J The number of sites.
#' @param gamma The average treatment effect under the alternative.
#' @param d The standardized ??. Either d OR gamma and ssq.Y must be specified.
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
#' pss.multi.ate(m = 20, J = 10, m.ratio = 1.5, gamma = 3, rho0 = 0.095, rho1 = 0.048, ssq.Y = 48)
#' # Example 14.5
#' pss.multi.ate(m = 10, J = 19, d = 0.5, rho0 = 0, rho1 = 0.05)

pss.multi.ate <- function(m = NULL, m.ratio = 1, J = NULL, gamma = NULL,
                          rho0 = NULL, rho1 = NULL, ssq.Y = NULL, d = NULL,
                          alpha = 0.05, sides = 2) {

  N <- m * J
  c <- (1 + m.ratio)^2 / m.ratio
  df2 <- J - 1

  if (is.null(d)) {
    d <- gamma / sqrt(ssq.Y)
    ncp <- d / sqrt(c * (1 - rho0 + (4 * m / c - 1) * rho1) / N)
    crit <- stats::qt(1 - alpha / sides, df2)
    power <- 1 - stats::pt(crit, df2, ncp)
  } else {
    DE <- 1 - rho0 + (m - 1) * rho1
    ncp <- (d / sqrt(c * (DE + (4 * m / c - 1) * rho1) / N))^2
    crit <- stats::qf(1 - alpha, 1, df2)
    power <- 1 - stats::pf(crit, 1, df2, ncp)
  }

  # Generate output text
  METHOD <-"Power for test of average treatment effect"

  # Print output as a power.htest object depending on which inputs were given
  structure(list(m = m, m.ratio = m.ratio, J = J, gamma = gamma,
                 rho0 = rho0, rho1 = rho1, ssq.Y = ssq.Y,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
