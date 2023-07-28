#' Power calculations for two sample t tests allowing for unequal sample sizes and/or variances
#'
#' @param n1 The sample size for group 1.
#' @param n2 The sample size for group 2; defaults to the same as n1.
#' @param delta DeltaA (the true difference mu1 - mu2) - Delta0 (the difference under the null).
#' @param sigma1 The estimated standard deviation for group 1; defaults to 1.
#' @param sigma2 The estimated standard deviation for group 2; defaults to the same as sigma2.
#' @param df.method Method for calculating the degrees of freedom: "welch" (default) or "classical".
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param one.or.two.sided Either "one" or "two" (default) to specify a one- or two- sided hypothesis test.
#' @param strict Use strict interpretation in two-sided case; defaults to TRUE.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 3.7
#' pss.t.test.2samp(n1 = 50, delta = 2, sigma1 = 5, one.or.two.sided = "one")
#' # Example 3.9
#' pss.t.test.2samp(n1 = NULL, delta = 0.5, sigma1 = 1, power = 0.8, one.or.two.sided = "two")

pss.t.test.2samp <- function (n1 = NULL, n2 = NULL, delta = NULL,
                              sigma1 = 1, sigma2 = 1,
                              df.method = c("welch", "classical"),
                              alpha = 0.05, power = NULL,
                              one.or.two.sided = c("two", "one"), strict = TRUE) {

  # If n2 or sigma2 aren't specified, default to the same as n1 and sigma1
  if (!is.null(n1) & is.null(n2))
    n2 <- n1
  if (!is.null(sigma1) & is.null(sigma2))
    sigma2 <- sigma1

  # Check if the arguments are specified correctly
  if (sum(sapply(list(n1, delta, sigma1, power, alpha), is.null)) != 1)
    stop("exactly one of n1, delta, sigma1, alpha, and power must be NULL")

  # Assign df method and number of sides
  df.method <- match.arg(df.method)
  one.or.two.sided <- match.arg(one.or.two.sided)
  side <- switch(one.or.two.sided, one = 1, two = 2)

  # Use absolute value of the effect size
  if (!is.null(delta))
    delta <- abs(delta)

  # welch df =
  # classical df = (1 + n1/n2) * n1 - 2
  p.body <- quote({
    stats::pt(stats::qt(alpha / side, n - 1, lower.tail = FALSE), n - 1,
              sqrt(n) * delta / sigma, lower.tail = FALSE)
  })

  if (strict & side == 2)
    p.body <- quote({
      qu <- stats::qt(alpha / side, n - 1, lower.tail = FALSE)
      ncp <- sqrt(n) * delta / sigma
      stats::pt(qu, n - 1, ncp, lower.tail = FALSE) +
        stats::pt(-qu, n - 1, ncp, lower.tail = TRUE)
    })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+07))$root
  else if (is.null(sigma))
    sigma <- uniroot(function(sigma) eval(p.body) - power, delta * c(1e-07, 1e+07))$root
  else if (is.null(delta))
    delta <- uniroot(function(delta) eval(p.body) - power,  c(1e-07, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "One-sample t test power calculation"

  # Print output as a power.htest object
  structure(list(n = n, delta = delta, sigma = sigma, alpha = alpha,
                 power = power, one.or.two.sided = one.or.two.sided,
                 method = METHOD), class = "power.htest")
}
