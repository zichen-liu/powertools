#' Power calculations for signed-rank test
#'
#' @param N The sample size; number of observations or paired differences.
#' @param ps The probability that the sum of two values exceeds zero when the alternative hypothesis is true.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two-sided hypothesis test.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' signedrank(N = 20, p = 0.87, power = NULL, sides = 2)


signedrank <- function (N = NULL, ps = NULL,
                            alpha = 0.05, power = NULL,
                            sides = 2) {

  # Check if the arguments are specified correctly
  if (sides != 1 & sides != 2)
    stop("please specify 1 or 2 sides")
  if (sum(sapply(list(N, power, alpha), is.null)) != 1)
    stop("exactly one of N, alpha, and power must be NULL")

  # power equation
  p.body <- quote({
    stats::pnorm(stats::qnorm(alpha / sides) + sqrt(3) * sqrt(N) * abs(ps - 0.5))
  })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(N))
    N <- stats::uniroot(function(N) eval(p.body) - power, c(2, 1e+07))$root
  else if (is.null(alpha))
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "Signed-rank test power calculation (normal approximation)"

  # Print output as a power.htest object
  structure(list(N = N, ps = ps, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")
}
