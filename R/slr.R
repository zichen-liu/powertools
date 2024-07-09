#' Power calculations for a simple linear regression
#'
#' @param N The sample size.
#' @param beta10 The slope regression coefficient under the null hypothesis; defaults to 0.
#' @param beta1A The slope regression coefficient under the alternative hypothesis.
#' @param sd.x.sq The sample variance of the covariate X.
#' @param sigma.e The estimated standard deviation of the error terms.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Yi = beta0 + beta1 * Xi + ei, i = 1,...,N
#' slr(N = 100, beta10 = 1, beta1A = 1.5, sd.x.sq = 25, sigma.e = 10, sides = 1)

slr <- function (N = NULL, beta10 = 0, beta1A = NULL,
                 sd.x.sq = NULL, sigma.e = NULL,
                 alpha = 0.05, power = NULL, sides = 2,
                 v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(N, alpha, power), "oneof")
  check(N, "pos"); check(N, "min", min = 4)
  check(alpha, "unit")
  check(power, "unit")
  check(beta10, "req"); check(beta10, "num")
  check(beta1A, "req"); check(beta1A, "num")
  check(sd.x.sq, "req"); check(sd.x.sq, "pos")
  check(sigma.e, "req"); check(sigma.e, "pos")
  check(sides, "req"); check(sides, "vals", valslist = c(1, 2))
  check(v, "req"); check(v, "bool")

  # Calculate power
  p.body <- quote({
    sxx <- (N - 1) * sd.x.sq
    ncp <- abs(beta1A - beta10) * sqrt(sxx) / sigma.e
    df <- N - 2
    1 - stats::pt(stats::qt(1 - alpha / sides, df), df, ncp)
  })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(N)) {
    N <- stats::uniroot(function(N) eval(p.body) - power, c(4, 1e+09))$root
    if (!v) return(N)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Power calculation for a simple linear regression"

  # Print output as a power.htest object
  structure(list(N = N, beta10 = beta10, beta1A = beta1A,
                 sd.x.sq = sd.x.sq, sigma.e = sigma.e,
                 alpha = alpha, power = power, sides = sides,
                 method = METHOD), class = "power.htest")

}
