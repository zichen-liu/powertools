#' Number of sites for multisite trials with binary outcomes

#' @param m The number of subjects per site.
#' @param J The number of sites.
#' @param pc The probability of the outcome in the control condition.
#' @param pt The probability of the outcome in the treatment condition.
#' @param prop.t The proportion of subjects allocated to the treatment condition within each site; defaults to 0.5.
#' @param sigma.u Standard deviation of the treatment effect across sites.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' pss.multisite.bin(m = 30, J = 25, pc = 0.1, pt = 0.2, sigma.u = 0.4, power = NULL)

pss.multisite.bin <- function (m = NULL, J = NULL, prop.t = 0.5,
                               pc = NULL, pt = NULL, sigma.u = NULL,
                               alpha = 0.05, power = NULL, sides = 2) {

  # Check if the arguments are specified correctly
  if (sides != 1 & sides != 2)
    stop("please specify 1 or 2 sides")

  # Calculate power
  p.body <- quote({
    or <- (pt / (1 - pt)) / (pc / (1 - pc))
    gammaA <- abs(log(or))
    ssq.e <- (1 / 4) * (1 / ((1 - prop.t) * pc * (1 - pc)) + 1 /
                          (prop.t * pt * (1 - pt)))
    var <- 1.2 * (4 * (ssq.e / m + 0.25 * sigma.u^2)) / J
    za <- stats::qnorm(1 - alpha / sides)
    stats::pnorm(gammaA / sqrt(var) - za)
  })


  # Generate output text
  METHOD <-"Number of sites for multisite trials with binary outcomes"

  if (is.null(alpha))
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  if (is.null(power))
    power <- eval(p.body)
  if (is.null(J))
    J <- stats::uniroot(function(J) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
  if (is.null(m))
    m <- stats::uniroot(function(m) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root

  # Print output as a power.htest object depending on which inputs were given
  structure(list(m = m, prop.t = prop.t, J = J,
                 pc = pc, pt = pt, sigma.u = sigma.u,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
