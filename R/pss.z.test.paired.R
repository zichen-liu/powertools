#' Power calculations for paired z tests
#'
#' @param n The sample size; the number of pairs.
#' @param delta DeltaA (the true mean difference) - Delta0 (the mean difference under the null).
#' @param sd1 The estimated pre standard deviation; defaults to 1.
#' @param sd2 The estimated post standard deviation; defaults to 1.
#' @param rho The estimated correlation between pre and post measurements on the same individual.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param strict Use strict interpretation in two-sided case; defaults to TRUE.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 3.15 approximation
#' pss.z.test.paired(n = NULL, delta = 4, sd1 = 10, sd2 = 10, rho = 0.4, power = 0.8, sides = 2)

pss.z.test.paired <- function (n = NULL, delta = NULL,
                               sd1 = 1, sd2 = 1, rho = NULL,
                               alpha = 0.05, power = NULL,
                               sides = c(2, 1), strict = TRUE) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(n, delta, power, alpha), is.null)) != 1)
    stop("exactly one of n, delta, alpha, and power must be NULL")
  if (is.null(sd1) | is.null(sd2) | is.null(rho))
    stop("please specify sd1, sd2, and rho")

  # Calculate the standard deviation of differences within pairs
  sigmad <- sqrt(sd1^2 + sd2^2 - 2 * rho * sd1 * sd2)

  # Calculate test statistic
  p.body <- quote({
    d <- abs(delta)
    stats::pnorm(stats::qnorm(alpha / sides) + sqrt(n) * d / sigmad)
  })

  if (strict && sides == 2)
    p.body <- quote({
      d <- abs(delta)
      stats::pnorm(stats::qnorm(alpha / sides) + sqrt(n) * d / sigmad) +
        stats::pnorm(stats::qnorm(alpha / sides) - sqrt(n) * d / sigmad)
    })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+07))$root
  else if (is.null(delta))
    delta <- uniroot(function(delta) eval(p.body) - power,  c(1e-07, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "Paired z test power calculation"
  NOTE <- "n is the number of pairs"
  sd <- c(sd1, sd2)

  # Print output as a power.htest object
  structure(list(n = n, delta = delta, sd = sd, rho = rho,
                 alpha = alpha, power = power, sides = sides,
                 method = METHOD, note = NOTE), class = "power.htest")
}
