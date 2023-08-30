#' Power calculations for one sample z tests
#'
#' @param n The sample size.
#' @param delta muA (the true mean) - mu0 (the mean under the null).
#' @param sd The estimated standard deviation; defaults to 1.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param strict Use strict interpretation in two-sided case; defaults to TRUE.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 2.5
#' pss.z.test.1samp(n = NULL, delta = 6.5 - 5.7, sd = 2, power = 0.8, sides = 2)
#' # Example 2.7
#' pss.z.test.1samp(n = 40, delta = NULL, sd = 1, power = 0.9, sides = 2)
#' # Example 3.6
#' pss.z.test.1samp(n = NULL, delta = 0.6, sd = 1, power = 0.8, sides = 1)

pss.z.test.1samp <- function (n = NULL, delta = NULL, sd = 1,
                              alpha = 0.05, power = NULL,
                              sides = 2, strict = TRUE) {

  # Check if the arguments are specified correctly
  if (sides != 1 & sides != 2)
    stop("please specify 1 or 2 sides")
  if (sum(sapply(list(n, delta, sd, power, alpha), is.null)) != 1)
    stop("exactly one of n, delta, sd, alpha, and power must be NULL")

  # Calculate test statistic
  p.body <- quote({
    d <- abs(delta)
    stats::pnorm(stats::qnorm(alpha / sides) + sqrt(n) * d / sd)
  })
  if (strict && sides == 2)
    p.body <- quote({
      d <- abs(delta)
      stats::pnorm(stats::qnorm(alpha / sides) + sqrt(n) * d / sd) +
      stats::pnorm(stats::qnorm(alpha / sides) - sqrt(n) * d / sd)
    })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+07))$root
  else if (is.null(sd))
    sd <- uniroot(function(sd) eval(p.body) - power, delta * c(1e-07, 1e+07))$root
  else if (is.null(delta))
    delta <- uniroot(function(delta) eval(p.body) - power,  c(1e-07, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "One-sample z test power calculation"

  # Print output as a power.htest object
  structure(list(n = n, delta = delta, sd = sd, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")
}
