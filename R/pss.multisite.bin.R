#' Number of sites for multisite trials with binary outcomes

#' @param m The number of subjects per site.
#' @param alloc.ratio The allocation ratio of intervention/control per site; defaults to 1.
#' @param J The number of sites.
#' @param pc The probability of the outcome in the control condition.
#' @param pt The probability of the outcome in the treatment condition.
#' @param sigma.u Standard deviation of the treatment effect across sites.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.multisite.bin(m = 30, J = 25, pc = 0.1, pt = 0.2, sigma.u = 0.4, power = NULL)
#' pss.multisite.bin(m = 30, J = NULL, pc = 0.1, pt = 0.2, sigma.u = 0.4, power = 0.9)

pss.multisite.bin <- function (m = NULL, alloc.ratio = 1, J = NULL,
                               pc = NULL, pt = NULL, sigma.u = NULL,
                               alpha = 0.05, power = NULL, sides = 2) {

  # Check if the arguments are specified correctly
  if (sides != 1 & sides != 2)
    stop("please specify 1 or 2 sides")

  # Calculate power
  p.body <- quote({
    or <- (pt / (1 - pt)) / (pc / (1 - pc))
    gammaA <- abs(log(or))
    prop.t <- alloc.ratio / (1 + alloc.ratio)
    ssq.e <- (1 / 4) * (1 / ((1 - prop.t) * pc * (1 - pc)) + 1 /
                          (prop.t * pt * (1 - pt)))
    var <- 1.2 * (4 * (ssq.e / m + 0.25 * sigma.u^2)) / J
    za <- stats::qnorm(1 - alpha / sides)
    stats::pnorm(gammaA / sqrt(var) - za)
  })

  if (is.null(alpha))
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else if (is.null(power))
    power <- eval(p.body)
  else if (is.null(J))
    J <- stats::uniroot(function(J) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
  else if (is.null(m))
    m <- stats::uniroot(function(m) eval(p.body) - power, c(2 + 1e-10, 1e+07))$root
  else if (is.null(alloc.ratio))
    alloc.ratio <- uniroot(function(alloc.ratio) eval(p.body) - power, c(2/m, 1e+07))$root
  else if (is.null(prop.t))
    prop.t <- uniroot(function(prop.t) eval(p.body) - power, c(0 + 1e-10, 1 - 1e-10))$root

  # Generate output text
  METHOD <-"Number of sites for multisite trials with binary outcomes"
  NOTE <- "m is the subjects per site split as intervention, control"
  p <- c(pc, pt)
  c <- m / (alloc.ratio + 1)
  t <- alloc.ratio * c
  m <- c(t, c)

  # Print output as a power.htest object depending on which inputs were given
  structure(list(m = m, J = J, `pc, pt` = p, sigma.u = sigma.u,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
