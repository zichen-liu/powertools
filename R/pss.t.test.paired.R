#' Power calculations for paired t tests
#'
#' @param N The sample size; the number of pairs.
#' @param delta DeltaA (the true mean difference) - Delta0 (the mean difference under the null).
#' @param sd1 The estimated pre standard deviation; defaults to 1.
#' @param sd2 The estimated post standard deviation; defaults to 1.
#' @param rho The estimated correlation between pre and post measurements on the same individual.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.t.test.paired(N = NULL, delta = 4, sd1 = 10, sd2 = 10, rho = 0.4, power = 0.8, sides = 2)

pss.t.test.paired <- function (N = NULL, delta = NULL,
                               sd1 = 1, sd2 = 1, rho = NULL,
                               alpha = 0.05, power = NULL, sides = 2,
                               v = TRUE) {

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
  pss.check(v, "req"); pss.check(v, "bool")

  # Calculate the standard deviation of differences within pairs
  sigmad <- sqrt(sd1^2 + sd2^2 - 2 * rho * sd1 * sd2)

  # Calculate df and ncp
  if (sides == 1)
    p.body <- quote({
      d <- abs(delta) / sigmad
      df <- N - 1
      stats::pt(stats::qt(alpha, df, lower.tail = FALSE), df,
                sqrt(N) * d, lower.tail = FALSE)
    })
  else if (sides == 2)
    p.body <- quote({
      d <- abs(delta)
      ncp <- sqrt(N) * d / sigmad
      df2 <- N - 1
      stats::pf(stats::qf(alpha, 1, df2, lower.tail = FALSE),
                1, df2, ncp^2, lower.tail = FALSE)
    })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(N)) {
    N <- stats::uniroot(function(N) eval(p.body) - power, c(2, 1e+07))$root
    if (!v) return(N)
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
  METHOD <- "Paired t test power calculation"
  NOTE <- "N is the number of pairs"
  sd <- c(sd1, sd2)

  # Print output as a power.htest object
  structure(list(N = N, delta = delta, `sd1, sd2` = sd, rho = rho,
                 alpha = alpha, power = power, sides = sides,
                 method = METHOD, note = NOTE), class = "power.htest")
}
