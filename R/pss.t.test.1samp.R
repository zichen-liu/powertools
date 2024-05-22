#' Power calculations for one sample t tests
#'
#' @param N The sample size.
#' @param delta muA (the true mean) - mu0 (the mean under the null).
#' @param sd The estimated standard deviation; defaults to 1.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.t.test.1samp(N = 36, delta = 4.9 - 5.7, sd = 2, sides = 1)
#' pss.t.test.1samp(N = 36, delta = 6.3 - 5.7, sd = 2, sides = 1)
#' pss.t.test.1samp(N = 36, delta = 4.9 - 5.7, sd = 2, sides = 2)
#' pss.t.test.1samp(delta = 0.6, sd = 1, power = 0.8, sides = 1)

pss.t.test.1samp <- function (N = NULL, delta = NULL, sd = 1,
                              alpha = 0.05, power = NULL, sides = 2) {

  # Check if the arguments are specified correctly
  pss.check.many(list(N, delta, sd, alpha, power), "oneof")
  pss.check(N, "int")
  pss.check(delta, "num")
  pss.check(sd, "pos")
  pss.check(alpha, "unit")
  pss.check(power, "unit")
  pss.check(sides, "req"); pss.check(sides, "vals", valslist = c(1, 2))

  # Calculate df and ncp
  if (sides == 1)
    p.body <- quote({
      d <- abs(delta) / sd
      df <- N - 1
      stats::pt(stats::qt(alpha, df, lower.tail = FALSE), df,
                sqrt(N) * d, lower.tail = FALSE)
  })
  else if (sides == 2)
    p.body <- quote({
      d <- abs(delta) / sd
      ncp <- sqrt(N) * d
      df2 <- N - 1
      stats::pf(stats::qf(alpha, 1, df2, lower.tail = FALSE),
                1, df2, ncp^2, lower.tail = FALSE)
    })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(N))
    N <- stats::uniroot(function(N) eval(p.body) - power, c(2, 1e+07))$root
  else if (is.null(sd))
    sd <- stats::uniroot(function(sd) eval(p.body) - power, delta * c(1e-07, 1e+07))$root
  else if (is.null(delta))
    delta <- stats::uniroot(function(delta) eval(p.body) - power,  c(1e-07, 1e+07))$root
  else if (is.null(alpha))
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "One-sample t test power calculation"

  # Print output as a power.htest object
  structure(list(N = N, delta = delta, sd = sd, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")

}
