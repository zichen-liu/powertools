#' Power calculations for two sample t tests allowing for unequal sample sizes and/or variances
#'
#' @param n The sample size for group 1.
#' @param n.ratio The ratio n1/n2 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param delta DeltaA (the true difference mu1 - mu2) - Delta0 (the difference under the null).
#' @param sd The estimated standard deviation for group 1; defaults to 1 (equal standard deviations in the two groups).
#' @param sd.ratio The ratio sd1/sd2 between the standard deviations of the two groups.
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
#' pss.t.test.2samp(n = 50, delta = 2, sd = 5, one.or.two.sided = "one")
#' # Example 3.9a
#' pss.t.test.2samp(n = NULL, n.ratio = 2, delta = 0.5, sd = 1, power = 0.8, one.or.two.sided = "two")
#' # Example 3.9b
#' pss.t.test.2samp(n = 49, n.ratio = 2, delta = 0.5, sd = 1, power = NULL, one.or.two.sided = "two")
#' # Example 3.10
#' pss.t.test.2samp(n = 25, n.ratio = 3, delta = 3, sd = 4, sd.ratio = 1.5, alpha = 0.025, one.or.two.sided = "one")
#' # Example 3.11
#' pss.t.test.2samp(n = NULL, delta = 0.5, sd = 1, power = 0.8, one.or.two.sided = "two")

pss.t.test.2samp <- function (n = NULL, n.ratio = 1, delta = NULL,
                              sd = 1, sd.ratio = 1,
                              df.method = c("welch", "classical"),
                              alpha = 0.05, power = NULL,
                              one.or.two.sided = c("two", "one"), strict = TRUE) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(n, n.ratio, delta, sd, sd.ratio, power, alpha), is.null)) != 1)
    stop("exactly one of n, n.ratio, delta, sd, sd.ratio, power, and alpha must be NULL")

  # Assign df method and number of sides
  df.method <- match.arg(df.method)
  one.or.two.sided <- match.arg(one.or.two.sided)
  side <- switch(one.or.two.sided, one = 1, two = 2)

  # Use absolute value of the effect size
  if (!is.null(delta))
    delta <- abs(delta)

  # Calculate df and ncp
  p.body <- quote({
    nu <- switch(df.method,
                 welch = (sd^2 / n + (sd * sd.ratio)^2 / (n * n.ratio))^2 /
                 ((sd^2 / n)^2 / (n - 1) +
                 ((sd * sd.ratio)^2 / (n.ratio * n))^2 / (n * n.ratio - 1)),
                 classical = (1 + n.ratio) * n - 2)
    stats::pt(stats::qt(alpha / side, nu, lower = FALSE), nu,
              sqrt(n / (1 + sd.ratio^2 / n.ratio)) * delta / sd, lower = FALSE)
  })

  if (strict & side == 2)
    p.body <- quote({
      nu <- switch(df.method,
                   welch = (sd^2 / n + (sd * sd.ratio)^2 / (n * n.ratio))^2 /
                   ((sd^2 / n)^2 / (n - 1) +
                   ((sd * sd.ratio)^2 / (n.ratio * n))^2 / (n * n.ratio - 1)),
                   classical = (1 + n.ratio) * n - 2)
      qu <- stats::qt(alpha / side, nu, lower = FALSE)
      ncp <- sqrt(n / (1 + sd.ratio^2 / n.ratio)) * delta / sd
      stats::pt(qu, nu, ncp, lower = FALSE) + pt(-qu, nu, ncp, lower = TRUE)
    })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+07))$root
  else if (is.null(n.ratio))
    n.ratio <- uniroot(function(ratio) eval(p.body) - power,c(2/n, 1e+07))$root
  else if (is.null(sd))
    sd <- uniroot(function(sd) eval(p.body) - power, delta * c(1e-07, 1e+07))$root
  else if (is.null(sd.ratio))
    sd.ratio <- uniroot(function(sd.ratio) eval(p.body) - power, c(1e-07, 1e+07))$root
  else if (is.null(delta))
    delta <- uniroot(function(delta) eval(p.body) - power,  c(1e-07, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "Two-sample t test power calculation"
  NOTE <- "n is the number in each group"
  n <- c(n, n * n.ratio)
  sd <- c(sd, sd * sd.ratio)

  # Print output as a power.htest object
  structure(list(n = n, delta = delta, sd = sd, alpha = alpha,
                 power = power, one.or.two.sided = one.or.two.sided,
                 method = METHOD, note = NOTE), class = "power.htest")
}
