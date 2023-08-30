#' Power calculations for one sample proportion tests
#'
#' @param n The sample size.
#' @param pA The true proportion.
#' @param p0 The proportion under the null.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 6.1
#' pss.prop.1samp(n = NULL, p0 = 0.2, pA = 0.3, power = 0.8, sides = 1)
#' # Example 6.2
#' pss.prop.1samp(n = NULL, p0 = 0.4, pA = 0.5, power = 0.8, sides = 1)
#'

pss.prop.1samp <- function (n = NULL, p0 = NULL, pA = NULL, alpha = 0.05,
                            power = NULL, sides = 2) {

  # Check if the arguments are specified correctly
  if (sides != 1 & sides != 2)
    stop("please specify 1 or 2 sides")
  if (sum(sapply(list(n, power, alpha), is.null)) != 1)
    stop("exactly one of 'n', 'alpha', and 'power' must be NULL")

  # Calculate test statistic
  p.body <- quote({
    d <- abs(pA - p0)
    (stats::pnorm(sqrt(n) * d / sqrt(pA * (1 - pA)) -
                  stats::qnorm(alpha / sides, lower = FALSE)))
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2 + 1e-10, 1e+09))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error", domain = NA)

  # Generate output text
  METHOD <- "One sample proportion test power calculation"

  # Print output as a power.htest object
  structure(list(n = n, p0 = p0, pA = pA, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")
}

