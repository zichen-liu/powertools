#' Power for test of treatment effect in cluster randomized trials
#'
#' @param m The number of subjects per cluster.
#' @param J The number of clusters.
#' @param delta The difference between the intervention and control means in the outcome variable.
#' @param sd The total standard deviation of the outcome variable; defaults to 1.
#' @param rho The intraclass correlation coefficient; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.crt.parallel(m = 30, J = 16, delta = 0.4, rho = 0.05)

pss.crt.parallel <- function (m = NULL, J = NULL, delta = NULL, sd = 1,
                              rho = 0, alpha = 0.05, power = NULL, sides = 2) {

  # Calculate power
  p.body <- quote({
    N <- m * J
    df <- J - 2
    d <- delta / sd
    ncp <- d / sqrt(4 * (1 + (m - 1) * rho) / N)
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
  else if (is.null(delta))
    delta <- stats::uniroot(function(delta) eval(p.body) - power, c(1e-07, 1e+07))$root

  # Generate output text
  METHOD <- "Power for test of treatment effect in cluster randomized trials"
  NOTE <- "m is the subjects per site"

  # Print output as a power.htest object
  structure(list(m = m, J = J, delta = delta, sd = sd,
                 rho = rho, alpha = alpha, power = power,
                 method = METHOD, note = NOTE), class = "power.htest")

}

