#' Power calculations for two sample z tests allowing for unequal sample sizes and/or variances
#'
#' @param n The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param delta DeltaA (the true difference mu1 - mu2) - Delta0 (the difference under the null) - delta. See delta.sign for guidance on the sign of delta.
#' @param sd The estimated standard deviation for group 1; defaults to 1 (equal standard deviations in the two groups).
#' @param sd.ratio The ratio sd2/sd1 between the standard deviations of the two groups.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param strict Use strict interpretation in two-sided case; defaults to TRUE.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.z.test.2samp(n = NULL, n.ratio = 1, delta = 0.5, sd = 1, power = 0.8, sides = 2)

pss.z.test.2samp <- function (n = NULL, n.ratio = 1, delta = NULL,
                              sd = 1, sd.ratio = 1,
                              alpha = 0.05, power = NULL,
                              sides = 2, strict = TRUE) {

  # Check if the arguments are specified correctly
  if (sides != 1 & sides != 2)
    stop("please specify 1 or 2 sides")
  if (sum(sapply(list(n, n.ratio, delta, sd, sd.ratio, power, alpha), is.null)) != 1)
    stop("exactly one of n, n.ratio, delta, sd, sd.ratio, power, and alpha must be NULL")

  # Calculate test statistic
  p.body <- quote({
    d <- abs(delta)
    stats::pnorm(stats::qnorm(alpha / sides) +
                 d / sqrt((sd * sd.ratio)^2 / (n * n.ratio) + sd^2 / n))})
  if (strict & sides == 2)
    p.body <- quote({
      d <- abs(delta)
      stats::pnorm(stats::qnorm(alpha / sides) +
                   d / sqrt((sd * sd.ratio)^2 / (n * n.ratio) + sd^2 / n)) +
      stats::pnorm(stats::qnorm(alpha / sides) -
                   d / sqrt((sd * sd.ratio)^2 / (n * n.ratio) + sd^2 / n))
    })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+07))$root
  else if (is.null(n.ratio))
    n.ratio <- uniroot(function(ratio) eval(p.body) - power,c(2/n, 1e+07))$root
  else if (is.null(sd))
    sd <- uniroot(function(sd) eval(p.body) - power, delta * c(1e-07, 1e+07))$root
  else if (is.null(sd.ratio))
    sd.ratio <- uniroot(function(sd.ratio) eval(p.body) - power, c(1e-07, 1e+07))$root
  else if (is.null(delta))
    delta <- uniroot(function(delta) eval(p.body) - power,  c(1e-07, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "Two-sample z test power calculation"
  n <- c(n, n * n.ratio)
  sd <- c(sd, sd * sd.ratio)

  # Print output as a power.htest object
  structure(list(n = n, delta = delta, sd = sd, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")
}
