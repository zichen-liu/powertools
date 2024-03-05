#' Power for an individual randomized group treatment trial with continuous outcomes
#'
#' @param m The number of subjects per cluster in the treatment group.
#' @param J The number of clusters in the treatment group.
#' @param n The number of total participants in the control group.
#' @param delta The difference between the intervention and control means in the outcome variable.
#' @param sd The total standard deviation of the outcome variable in the control group; defaults to 1.
#' @param icc The intraclass correlation coefficient in the treatment group; defaults to 0.
#' @param Theta The ratio of the total variance between intervention and control groups; defaults to 1.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.irgtt.cont(m = 10, J = 10, n = 100, delta = 0.4, icc = 0.05, Theta = 1, sides = 1)
#' pss.irgtt.cont(m = 10, J = 12, delta = 0.4, icc = 0.05, Theta = 1, sides = 1, power = 0.8)

pss.irgtt.cont <- function (m = NULL, J = NULL, n = NULL, delta = NULL, sd = 1,
                              icc = 0, Theta = 1, alpha = 0.05, power = NULL, sides = 2) {

  # Calculate power
  p.body <- quote({
    df2 <- J + n
    de <- 1 + (m - 1) * icc
    d <- delta / sd
    lambda <- d / sqrt(Theta * de / (m * J) + 1 / n)
    crit <- stats::qf(1 - alpha / sides, df1 = 1, df2 = df2)
    1 - stats::pf(crit, df1 = 1, df2 = df2, ncp = lambda^2)
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
  else if (is.null(n))
    n <- stats::uniroot(function(n) eval(p.body) - power, c(4 + 1e-10, 1e+07))$root
  else if (is.null(delta))
    delta <- stats::uniroot(function(delta) eval(p.body) - power, c(1e-07, 1e+07))$root

  mjn <- c(m, J, n)

  # Print output as a power.htest object
  METHOD <- "Power for an individual randomized group treatment trial with continuous outcomes"
  structure(list(`m, J, n` = mjn, delta = delta, sd = sd, icc = icc,
                 Theta = Theta, alpha = alpha, power = power, sides = sides,
                 method = METHOD), class = "power.htest")
  # put m, J, n on the same line

}

