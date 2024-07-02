#' Power calculations for sign test (one sample median)
#'
#' @param N The sample size.
#' @param p The probability of a positive difference when the alternative hypothesis is true.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two-sided hypothesis test.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' signtest(N = 40, p = 0.7, power = NULL, alpha = 0.05, sides = 1)


signtest <- function (N = NULL, p = NULL, alpha = 0.05, power = NULL,
                      sides = 2, v = FALSE) {

  # Check if the arguments are specified correctly
  if (sides != 1 & sides != 2)
    stop("please specify 1 or 2 sides")
  if (sum(sapply(list(N, power, alpha), is.null)) != 1)
    stop("exactly one of N, alpha, and power must be NULL")
  check(v, "req"); check(v, "bool")

  # power equation
  p.body <- quote({
    stats::pnorm(stats::qnorm(alpha / sides) + 2 * sqrt(N) * abs(p - 0.5))
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
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "One-sample sign test power calculation (normal approximation)"

  # Print output as a power.htest object
  structure(list(N = N, p = p, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")
}
