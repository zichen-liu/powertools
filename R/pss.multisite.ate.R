#' Power for test of average treatment effect
#'
#' @param m The number of subjects per site.
#' @param J The number of sites.
#' @param gamma The average treatment effect under the alternative.
#' @param ssq.Y The total variance of the outcome variable Y.
#' @param d The standardized effect size. Either gamma and ssq.Y OR d must be specified.
#' @param rho0 The proportion of total variance of the outcome attributable to variation in site-level means.
#' @param rho1 The proportion of total variance of the outcome attributable to variation in the treatment effect across sites.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' pss.multisite.ate(m = 20, J = 10, gamma = 3, ssq.Y = 40, rho0 = 0.1, rho1 = 0)
#' pss.multisite.ate(m = 20, J = 10, gamma = 3, ssq.Y = 48, rho0 = 0.095, rho1 = 0.048)
#' pss.multisite.ate(m = 10, J = 19, d = 0.5, rho0 = 0, rho1 = 0.05)

pss.multisite.ate <- function (m = NULL, J = NULL,
                               gamma = NULL, ssq.Y = NULL, d = NULL,
                               rho0 = NULL, rho1 = NULL,
                               alpha = 0.05, sides = 2) {

  # Calculate power
  N <- m * J
  df <- J - 1 # not affected by covariate

  if (is.null(d)) {
    d <- gamma / sqrt(ssq.Y)
    ncp <- d / sqrt(4 * (1 - rho0) / N) # divide by RE
    crit <- stats::qt(1 - alpha / sides, df)
    power <- 1 - stats::pt(crit, df, ncp)
  } else {
    DE <- 1 - rho0 + (m - 1) * rho1
    ncp <- (d / sqrt(4 * DE / N))^2  # divide by RE
    crit <- stats::qf(1 - alpha, 1, df)
    power <- 1 - stats::pf(crit, 1, df, ncp)
  }

  # Generate output text
  METHOD <-"Power for test of average treatment effect"

  # Print output as a power.htest object depending on which inputs were given
  if (is.null(d)) {
    structure(list(m = m, J = J, gamma = gamma, ssq.Y = ssq.Y,
                   rho0 = rho0, rho1 = rho1,
                   alpha = alpha, power = power,
                   method = METHOD), class = "power.htest")
  } else {
    structure(list(m = m, J = J, d = d, rho0 = rho0, rho1 = rho1,
                   alpha = alpha, power = power,
                   method = METHOD), class = "power.htest")
  }

}
