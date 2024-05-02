#' Power calculations for two sample proportion tests
#'
#' @param n1 The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param p1 The proportion in group 1.
#' @param p2 The proportion in group 2.
#' @param margin The margin of noninferiority or superiority; defaults to 0. See margin.sign for guidance on the sign of margin.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.prop.2samp(n1 = NULL, p1 = 0.6, p2 = 0.8, alpha = 0.025, power = 0.9, sides = 1)
#' pss.prop.2samp(n1 = NULL, p1 = 0.25, p2 = 0.25, margin = 0.1, alpha = 0.025, power = 0.8, sides = 1)

pss.prop.2samp <- function (n1 = NULL, n.ratio = 1, p1 = NULL, p2 = NULL, margin = 0,
                            alpha = 0.05, power = NULL, sides = 2) {

  # Check if the arguments are specified correctly
  pss.check.many(list(n1, n.ratio, alpha, power), "oneof")
  pss.check(n1, "int")
  pss.check(n.ratio, "pos")
  pss.check(p1, "req"); pss.check(p1, "unit")
  pss.check(p2, "req"); pss.check(p2, "unit")
  pss.check(margin, "req"); pss.check(margin, "num")
  pss.check(alpha, "unit")
  pss.check(power, "unit")
  pss.check(sides, "req"); pss.check(sides, "vals", valslist = c(1, 2))

  # Calculate test statistic
  p.body <- quote({
    d <- abs(p1 - p2) - margin
    q1 <- 1 - p1
    q2 <- 1 - p2
    ((stats::qnorm(alpha / sides) + stats::qnorm(1 - power))^2 *
    (n.ratio * p1 * q1 + p2 * q2) / (n.ratio * d^2))
  })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(n1))
    n1 <- eval(p.body)
  else if (is.null(power))
    power <- stats::uniroot(function(power) eval(p.body) - n1, c(1e-05, 0.99999))$root
  else if (is.null(n.ratio))
    n.ratio <- stats::uniroot(function(n.ratio) eval(p.body) - n1, c(2/n1, 1e+07))$root
  else if (is.null(alpha))
    alpha <- stats::uniroot(function(alpha) eval(p.body) - n1, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "Two sample comparison of proportions power calculation"
  n <- c(n1, n1 * n.ratio)
  p <- c(p1, p2)

  # Print output as a power.htest object
  structure(list(`n1, n2` = n, `p1, p2` = p, margin = margin,
                 alpha = alpha, power = power, sides = sides,
                 method = METHOD), class = "power.htest")
}

