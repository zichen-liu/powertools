#' Power for test of average treatment effect
#'
#' @param m The number of subjects per site or the mean cluster size (if unequal number of participants per site).
#' @param m.sd The standard deviation of cluster sizes (provide if unequal number of participants per site); defaults to 0.
#' @param m.ratio The allocation ratio per site; defaults to 1.
#' @param J The number of sites.
#' @param delta The difference between the intervention and control means in the outcome variable.
#' @param sd The total standard deviation of the outcome variable; defaults to 1.
#' @param rho0 The proportion of total variance of the outcome attributable to variation in site-level means.
#' @param rho1 The proportion of total variance of the outcome attributable to variation in the treatment effect across sites.
#' @param Rsq The estimated R^2 for regressing the outcome on the covariates; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.multisite.ate(m = 20, J = 10, delta = 3, sd = sqrt(40), rho0 = 0.1, rho1 = 0)
#' pss.multisite.ate(m = 20, J = 10, delta = 3, sd = sqrt(48), rho0 = 0.095, rho1 = 0.048)
#' pss.multisite.ate(m = 20, m.ratio = 1.5, J = 10, delta = 0.43, sd = 1, rho0 = 0.095, rho1 = 0.048)
#' pss.multisite.ate(m = 10, J = NULL, delta = 0.5, sd = 1, rho0 = 0, rho1 = 0.05, power = 0.8)
#' pss.multisite.ate(m = 20, m.sd = 5, J = 10, delta = 3, sd = sqrt(48), rho0 = 0.095, rho1 = 0.048)
#' pss.multisite.ate(m = 20, J = 10, delta = 3, sd = sqrt(48), rho0 = 0.095, rho1 = 0.048, Rsq = 0.5^2)

pss.multisite.ate <- function (m = NULL, m.sd = 0, m.ratio = 1, J = NULL,
                               delta = NULL, sd = 1,
                               rho0 = NULL, rho1 = NULL, Rsq = 0,
                               alpha = 0.05, power = NULL, sides = 2) {

  # Calculate power
  p.body <- quote({
    N <- m * J
    df <- J - 1
    d <- delta / (sd * sqrt((1 - Rsq)))
    RE <- pss.multisite.re(m.mean = m, m.sd = m.sd, rho = rho1)$re
    c <- (1 + m.ratio)^2 / m.ratio
    ncp <- d / sqrt(c * (1 - rho0 + (4 * m / c - 1) * rho1) / N) / RE
    crit <- stats::qt(1 - alpha / sides, df)
    1 - stats::pt(crit, df, ncp)
  })

  # Use uniroot to calculate missing argument
  if (is.null(alpha))
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else if (is.null(power))
    power <- eval(p.body)
  else if (is.null(J))
    J <- stats::uniroot(function(J) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
  else if (is.null(m))
    m <- stats::uniroot(function(m) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
  else if (is.null(m.ratio))
    m.ratio <- stats::uniroot(function(m.ratio) eval(p.body) - power, c(2/m, 1e+07), extendInt = "yes", maxiter = 5000)$root
  else if (is.null(delta))
    delta <- stats::uniroot(function(delta) eval(p.body) - power, c(1e-07, 1e+07))$root

  # Generate output text
  METHOD <-"Power for test of average treatment effect"
  rho <- c(rho0, rho1)
  m <- c(m, m * m.ratio)

  # Print output as a power.htest object
  structure(list(m = m, J = J, delta = delta, sd = sd,
                 `rho0, rho1` = rho, Rsq = Rsq,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}

