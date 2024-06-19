#' Power calculations for two sample z tests allowing for unequal sample sizes and/or variances
#'
#' @param n1 The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param delta DeltaA (the true difference mu1 - mu2) - Delta0 (the difference under the null) - margin See margin.sign for guidance on the sign of the margin.
#' @param sd1 The estimated standard deviation for group 1; defaults to 1 (equal standard deviations in the two groups).
#' @param sd.ratio The ratio sd2/sd1 between the standard deviations of the two groups.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' z.test.2samp(n1 = NULL, n.ratio = 1, delta = 0.5, sd1 = 1, power = 0.8, sides = 2)

z.test.2samp <- function (n1 = NULL, n.ratio = 1, delta = NULL,
                          sd1 = 1, sd.ratio = 1,
                          alpha = 0.05, power = NULL, sides = 2,
                          v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(n1, n.ratio, delta, sd1, sd.ratio, alpha, power), "oneof")
  check(n1, "pos")
  check(n.ratio, "pos")
  check(sd1, "pos")
  check(sd.ratio, "pos")
  check(delta, "num")
  check(alpha, "unit")
  check(power, "unit")
  check(sides, "req"); check(sides, "vals", valslist = c(1, 2))
  check(v, "req"); check(v, "bool")

  # Calculate test statistic
  if (sides == 1)
    p.body <- quote({
      d <- abs(delta)
      stats::pnorm(stats::qnorm(alpha) +
                   d / sqrt((sd1 * sd.ratio)^2 / (n1 * n.ratio) + sd1^2 / n1))})
  else if (sides == 2)
    p.body <- quote({
      d <- abs(delta)
      stats::pnorm(stats::qnorm(alpha / 2) +
                   d / sqrt((sd1 * sd.ratio)^2 / (n1 * n.ratio) + sd1^2 / n1)) +
      stats::pnorm(stats::qnorm(alpha / 2) -
                   d / sqrt((sd1 * sd.ratio)^2 / (n1 * n.ratio) + sd1^2 / n1))
    })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(n1)) {
    n1 <- stats::uniroot(function(n1) eval(p.body) - power, c(2, 1e+07))$root
    if (!v) return(n1)
  }
  else if (is.null(n.ratio)) {
    n.ratio <- stats::uniroot(function(n.ratio) eval(p.body) - power,c(2/n1, 1e+07))$root
    if (!v) return(n.ratio)
  }
  else if (is.null(sd1)) {
    sd1 <- stats::uniroot(function(sd1) eval(p.body) - power, delta * c(1e-07, 1e+07))$root
    if (!v) return(sd1)
  }
  else if (is.null(sd.ratio)) {
    sd.ratio <- stats::uniroot(function(sd.ratio) eval(p.body) - power, c(1e-07, 1e+07))$root
    if (!v) return(sd.ratio)
  }
  else if (is.null(delta)) {
    delta <- stats::uniroot(function(delta) eval(p.body) - power,  c(1e-07, 1e+07))$root
    if (!v) return(delta)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Two-sample z test power calculation"
  n <- c(n1, n1 * n.ratio)
  sd <- c(sd1, sd1 * sd.ratio)

  # Print output as a power.htest object
  structure(list(`n1, n2` = n, delta = delta, `sd1, sd2` = sd, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")
}
