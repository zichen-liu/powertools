#' Power for a cluster randomized trial with a binary outcome
#'
#' @param m The number of subjects per cluster
#' @param J The number of clusters
#' @param pc The probability of the outcome in control clusters
#' @param pt The probability of the outcome in treatment clusters
#' @param sigma.u Standard deviation of the cluster random effect.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.crt.parallel.bin(m = 60, J = NULL, pc = 0.25, pt = 0.15, sigma.u = 0.3, power = 0.8)

pss.crt.parallel.bin <- function (m = NULL, J = NULL,
                                  pc = NULL, pt = NULL, sigma.u = NULL,
                                  alpha = 0.05, power = NULL, sides = 2) {

  # Check if the arguments are specified correctly
  if (sides != 1 & sides != 2)
    stop("please specify 1 or 2 sides")

  # Calculate power
  p.body <- quote({
    or <- (pt / (1 - pt)) / (pc / (1 - pc))
    gammaA <- abs(log(or))
    ssq.e <- (1 / 2) * (1 / (pc * (1 - pc)) + 1 / (pt * (1 - pt)))
    var <- 1.2 * (4 * (ssq.e  + m * sigma.u^2)) / (m *J)
    za <- stats::qnorm(1 - alpha / sides)
    stats::pnorm(gammaA / sqrt(var) - za)
  })

  if (is.null(alpha))
    alpha <- stats::uniroot(function(alpha) eval(p.body)  - power, c(1e-10, 1 - 1e-10))$root
  else if (is.null(power))
    power <- eval(p.body)
  else if (is.null(J))
    J <- stats::uniroot(function(J) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
  else if (is.null(m))
    m <- stats::uniroot(function(m) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root

  # Generate output text
  METHOD <-"Power for a cluster randomized trial with a binary outcome"
  p <- c(pc, pt)


  # Print output as a power.htest object depending on which inputs were given
  structure(list(m = m, J = J, `pc, pt` = p, sigma.u = sigma.u,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
