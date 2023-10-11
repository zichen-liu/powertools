#' Number of sites for multisite trials with binary outcomes

#' @param m The number of subjects per site.
#' @param pc The probability of the outcome in the control condition.
#' @param pt The probability of the outcome in the treatment condition.
#' @param sigma.u Standard deviation of the treatment effect across sites.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' pss.multisite.bin(m = 30, pc = 0.1, pt = 0.2, sigma.u = 0.4, power = 0.9)

pss.multisite.bin <- function (m = NULL, pc = NULL, pt = NULL, sigma.u = NULL,
                               alpha = 0.05, power = NULL, sides = 2) {

  # Check if the arguments are specified correctly
  if (sides != 1 & sides != 2)
    stop("please specify 1 or 2 sides")

  # Calculate J
  or <- (pt / (1 - pt)) / (pc / (1 - pc))
  gamma <- log(or)
  za <- stats::qnorm(1 - alpha / sides)
  zb <- stats::qnorm(power)
  ssq.e <- 0.5 * (1 / (pc * (1 - pc)) + 1 / (pt * (1 - pt)))
  Jvar <- 1.2 * (4 * (ssq.e / m + 0.25 * sigma.u^2))
  J <- (za + zb)^2 * Jvar / gamma^2

  # Generate output text
  METHOD <-"Number of sites for multisite trials with binary outcomes"

  # Print output as a power.htest object depending on which inputs were given
  structure(list(m = m, J = J, pc = pc, pt = pt, sigma.u = sigma.u,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
