#' Power for individual randomized group treatment trials with binary outcomes
#'
#' @param m The number of subjects per cluster in the treatment group.
#' @param J The number of clusters in the treatment group.
#' @param n The number of total participants in the control group.
#' @param p1 The probability of the outcome in the control group.
#' @param p2 The probability of the outcome in the intervention group.
#' @param icc The intraclass correlation coefficient in the treatment group; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.irgtt.bin(m = 20, J = 5, n = 100, p1 = 0.8, p2 = 0.6, icc = 0.04, sides = 2)
#' pss.irgtt.bin(m = 20, J = 6, n = 120, p1 = 0.8, p2 = 0.6, icc = 0.04, sides = 2)

pss.irgtt.bin <- function (m = NULL, J = NULL, n = NULL, p1 = NULL, p2 = NULL,
                            icc = 0, alpha = 0.05, power = NULL, sides = 2) {

  # Calculate power
  p.body <- quote({
    de <- 1 + (m - 1) * icc
    vard <- p2 * (1 - p2) * de / (m * J) + p1 * (1 - p1) / n
    ncp <- abs(p2 - p1) / sqrt(vard)
    crit <- stats::qnorm(alpha / sides, lower.tail = F)
    stats::pnorm(crit, ncp, lower.tail = F)
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

  # Print output as a power.htest object
  METHOD <- "Power for individual randomized group treatment trials with binary outcomes"
  structure(list(m = m, J = J, n = n, p1 = p1, p2 = p2, icc = icc,
                 alpha = alpha, power = power, sides = sides,
                 method = METHOD), class = "power.htest")

}

