#' Power calculations for a simple linear regression
#'
#' @param n The sample size.
#' @param beta10 The slope regression coefficient under the null hypothesis.
#' @param beta1A The slope regression coefficient under the alternative hypothesis.
#' @param sd.x.sq The sample variance of the covariate X.
#' @param sigma.e The estimated standard deviation of the error terms.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples#' # Example 10.3
#' Yi = beta0 + beta1 * Xi + ei, i = 1,...,n
#' pss.slr(n = 100, beta10 = 1, beta1A = 1.5, sd.x.sq = 25, sigma.e = 10, sides = 1)

pss.slr <- function (n = NULL, beta10 = 0, beta1A = NULL,
                     sd.x.sq = NULL, sigma.e = NULL,
                     alpha = 0.05, power = NULL, sides = 2) {

  # Check if the arguments are specified correctly
  if (sides != 1 & sides != 2)
    stop("please specify 1 or 2 sides")
  if (sum(sapply(list(n, power, alpha), is.null)) != 1)
    stop("exactly one of n, alpha, and power must be NULL")
  if (is.null(beta1A))
    stop("please specify beta1A")
  if (is.null(sd.x.sq) | is.null(sigma.e))
    stop("please specify estimated standard deviations (see documentation)")
  if (!is.null(n) && any(n < 4))
    stop("number of observations must be at least 4")

  # Calculate power
  p.body <- quote({
    sxx <- (n - 1) * sd.x.sq
    ncp <- abs(beta1A - beta10) * sqrt(sxx) / sigma.e
    df <- n - 2
    1 - stats::pt(stats::qt(1 - alpha / sides, df), df, ncp)
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(4, 1e+09))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "Power calculation for a simple linear regression"

  # Print output as a power.htest object
  structure(list(n = n, beta10 = beta10, beta1A = beta1A,
                 sd.x.sq = sd.x.sq, sigma.e = sigma.e,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
