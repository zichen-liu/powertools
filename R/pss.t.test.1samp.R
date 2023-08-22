#' Power calculations for one sample t tests
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
#' # Example 3.2
#' pss.t.test.1samp(n = 36, delta = 4.9 - 5.7, sd = 2, sides = 1)
#' # Example 3.3
#' pss.t.test.1samp(n = 36, delta = 6.3 - 5.7, sd = 2, sides = 1)
#' # Example 3.5
#' pss.t.test.1samp(n = 36, delta = 4.9 - 5.7, sd = 2, sides = 2)
#' # Example 3.6
#' pss.t.test.1samp(delta = 0.6, sd = 1, power = 0.8, sides = 1)

pss.t.test.1samp <- function (n = NULL, delta = NULL, sd = 1,
                              alpha = 0.05, power = NULL,
                              sides = c(2, 1), strict = TRUE) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(n, delta, sd, power, alpha), is.null)) != 1)
    stop("exactly one of n, delta, sd, alpha, and power must be NULL")

  # Calculate df and ncp
  p.body <- quote({
    d <- abs(delta)
    stats::pt(stats::qt(alpha / sides, n - 1, lower.tail = FALSE), n - 1,
              sqrt(n) * d / sd, lower.tail = FALSE)
  })
  if (strict && sides == 2)
    p.body <- quote({
      d <- abs(delta)
      ncp <- sqrt(n) * d / sd
      stats::pf(stats::qf(alpha, 1, n - 1, lower.tail = FALSE),
                1, n - 1, ncp^2, lower.tail = FALSE)
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
  METHOD <- "One-sample t test power calculation"

  # Print output as a power.htest object
  structure(list(n = n, delta = delta, sd = sd, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")
}
