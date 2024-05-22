#' Power calculations for paired z tests
#'
#' @param N The sample size; the number of pairs.
#' @param delta DeltaA (the true mean difference) - Delta0 (the mean difference under the null).
#' @param sd1 The estimated pre standard deviation; defaults to 1.
#' @param sd2 The estimated post standard deviation; defaults to 1.
#' @param rho The estimated correlation between pre and post measurements on the same individual; defaults to 0.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.z.test.paired(N = NULL, delta = 4, sd1 = 10, sd2 = 10, rho = 0.4, power = 0.8, sides = 2)

pss.z.test.paired <- function (N = NULL, delta = NULL,
                               sd1 = 1, sd2 = 1, rho = NULL,
                               alpha = 0.05, power = NULL, sides = 2) {

  # Check if the arguments are specified correctly
  pss.check.many(list(N, delta, alpha, power), "oneof")
  pss.check(N, "int")
  pss.check(sd1, "req"); pss.check(sd1, "pos")
  pss.check(sd2, "req"); pss.check(sd2, "pos")
  pss.check(rho, "req"); pss.check(rho, "uniti")
  pss.check(delta, "num")
  pss.check(alpha, "unit")
  pss.check(power, "unit")
  pss.check(sides, "req"); pss.check(sides, "vals", valslist = c(1, 2))

  # Calculate the standard deviation of differences within pairs
  sigmad <- sqrt(sd1^2 + sd2^2 - 2 * rho * sd1 * sd2)

  # Calculate test statistic
  if (sides == 1)
    p.body <- quote({
      d <- abs(delta) / sigmad
      stats::pnorm(stats::qnorm(alpha) + sqrt(N) * d)
    })

  else if (sides == 2)
    p.body <- quote({
      d <- abs(delta) / sigmad
      stats::pnorm(stats::qnorm(alpha / 2) + sqrt(N) * d) +
        stats::pnorm(stats::qnorm(alpha / 2) - sqrt(N) * d)
    })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(N))
    N <- stats::uniroot(function(N) eval(p.body) - power, c(2, 1e+07))$root
  else if (is.null(delta))
    delta <- stats::uniroot(function(delta) eval(p.body) - power,  c(1e-07, 1e+07))$root
  else if (is.null(alpha))
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "Paired z test power calculation"
  NOTE <- "N is the number of pairs"
  sd <- c(sd1, sd2)

  # Print output as a power.htest object
  structure(list(N = N, delta = delta, `sd1, sd2` = sd, rho = rho,
                 alpha = alpha, power = power, sides = sides,
                 method = METHOD, note = NOTE), class = "power.htest")
}
