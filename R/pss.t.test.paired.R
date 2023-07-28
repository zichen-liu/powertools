#' Power calculations for paired t tests
#'
#' @param n The sample size; the number of pairs.
#' @param delta DeltaA (the true mean difference) - Delta0 (the mean difference under the null).
#' @param sigma.pre The estimated pre standard deviation; defaults to 1.
#' @param sigma.post The estimated post standard deviation; defaults to the same as sigma.pre.
#' @param rho The estimated correlation between pre and post measurements on the same individual.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param one.or.two.sided Either "one" or "two" (default) to specify a one- or two- sided hypothesis test.
#' @param strict Use strict interpretation in two-sided case; defaults to TRUE.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 3.15
#' pss.t.test.paired(n = NULL, delta = 4, sigma.pre = 10, sigma.post = 10, rho = 0.4, power = 0.8, one.or.two.sided = "two")

pss.t.test.paired <- function (n = NULL, delta = NULL,
                               sigma.pre = 1, sigma.post = NULL, rho = NULL,
                               alpha = 0.05, power = NULL,
                               one.or.two.sided = c("two", "one"), strict = TRUE) {

  # If sigma.post isnt' specified, default to the same as sigma.pre
  if (!is.null(sigma.pre) && is.null(sigma.post))
    sigma.post <- sigma.pre

  # Check if the arguments are specified correctly
  if (sum(sapply(list(n, delta, sigma.pre, power, alpha), is.null)) != 1)
    stop("exactly one of n, delta, sigma.pre, alpha, and power must be NULL")
  if (is.null(rho))
    stop("please specify a correlation rho")

  # Assign number of sides
  one.or.two.sided <- match.arg(one.or.two.sided)
  side <- switch(one.or.two.sided, one = 1, two = 2)

  # Use absolute value of the effect size
  if (!is.null(delta))
    delta <- abs(delta)

  # Calculate the standard deviation of differences within pairs
  sigma.d <- sqrt(sigma.pre^2 + sigma.post^2 - 2 * rho * sigma.pre * sigma.post)

  # Calculate df and ncp
  p.body <- quote({
    stats::pt(stats::qt(alpha / side, n - 1, lower.tail = FALSE), n - 1,
              sqrt(n) * delta / sigma.d, lower.tail = FALSE)
  })

  # The strict two-sided case uses the F distribution
  if (strict && side == 2)
    p.body <- quote({
      ncp <- sqrt(n) * delta / sigma.d
      stats::pf(stats::qf(alpha, 1, n - 1, lower.tail = FALSE),
                1, n - 1, ncp^2, lower.tail = FALSE)
    })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+07))$root
#  else if (is.null(sigma.pre)){
#    sigma.pre <- uniroot(function(sigma.pre) eval(p.body) - power, delta * c(1e-07, 1e+07))$root
#    sigma.ratio <- uniroot(function(sigma.ratio) eval(p.body) - power, c(1e-07, 1e+07))$root
#    sigma.post <- sigma.pre * sigma.ratio
#  }
  else if (is.null(delta))
    delta <- uniroot(function(delta) eval(p.body) - power,  c(1e-07, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "Paired t test power calculation"
  NOTE <- "n is the number of pairs"
  sigma <- c(sigma.pre, sigma.post)

  # Print output as a power.htest object
  structure(list(n = n, delta = delta, sigma = sigma, rho = rho,
                 alpha = alpha, power = power, one.or.two.sided = one.or.two.sided,
                 method = METHOD, note= NOTE), class = "power.htest")
}
