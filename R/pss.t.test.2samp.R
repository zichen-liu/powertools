#' Power calculations for two sample t tests allowing for unequal sample sizes and/or variances
#'
#' @param n1 The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param delta DeltaA (the true difference mu1 - mu2) - Delta0 (the difference under the null) - d. See delta.sign for guidance on the sign of d.
#' @param sd The estimated standard deviation for group 1; defaults to 1 (equal standard deviations in the two groups).
#' @param sd.ratio The ratio sd2/sd1 between the standard deviations of the two groups.
#' @param df.method Method for calculating the degrees of freedom: "welch" (default) or "classical".
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param strict Use strict interpretation in two-sided case; defaults to TRUE.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.t.test.2samp(n1 = 50, delta = 2, sd = 5, sides = 1)
#' pss.t.test.2samp(n1 = NULL, n.ratio = 2, delta = 0.5, sd = 1, power = 0.8, sides = 2)
#' pss.t.test.2samp(n1 = 49, n.ratio = 2, delta = 0.5, sd = 1, power = NULL, sides = 2)
#' pss.t.test.2samp(n1 = 25, n.ratio = 3, delta = 3, sd = 4, sd.ratio = 1.5, alpha = 0.025, sides = 1)
#' pss.t.test.2samp(n1 = NULL, delta = 0.5, sd = 1, power = 0.8, sides = 2)

pss.t.test.2samp <- function (n1 = NULL, n.ratio = 1, delta = NULL,
                              sd = 1, sd.ratio = 1,
                              df.method = c("welch", "classical"),
                              alpha = 0.05, power = NULL,
                              sides = 2, strict = TRUE) {

  # Check if the arguments are specified correctly
  if (sides != 1 & sides != 2)
    stop("please specify 1 or 2 sides")
  if (sum(sapply(list(n1, n.ratio, delta, sd, sd.ratio, power, alpha), is.null)) != 1)
    stop("exactly one of n1, n.ratio, delta, sd, sd.ratio, power, and alpha must be NULL")

  # Assign df method
  df.method <- match.arg(df.method)

  # Calculate df and ncp
  p.body <- quote({
    d <- abs(delta)
    nu <- switch(df.method,
                 welch = (sd^2 / n1 + (sd * sd.ratio)^2 / (n1 * n.ratio))^2 /
                 ((sd^2 / n1)^2 / (n1 - 1) +
                 ((sd * sd.ratio)^2 / (n.ratio * n1))^2 / (n1 * n.ratio - 1)),
                 classical = (1 + n.ratio) * n1 - 2)
    stats::pt(stats::qt(alpha / sides, nu, lower = FALSE), nu,
              sqrt(n1 / (1 + sd.ratio^2 / n.ratio)) * d / sd, lower = FALSE)
  })

  if (strict & sides == 2)
    p.body <- quote({
      d <- abs(delta)
      nu <- switch(df.method,
                   welch = (sd^2 / n1 + (sd * sd.ratio)^2 / (n1 * n.ratio))^2 /
                   ((sd^2 / n1)^2 / (n1 - 1) +
                   ((sd * sd.ratio)^2 / (n.ratio * n1))^2 / (n1 * n.ratio - 1)),
                   classical = (1 + n.ratio) * n1 - 2)
      qu <- stats::qt(alpha / sides, nu, lower = FALSE)
      ncp <- sqrt(n1 / (1 + sd.ratio^2 / n.ratio)) * d / sd
      stats::pt(qu, nu, ncp, lower = FALSE) + pt(-qu, nu, ncp, lower = TRUE)
    })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n1))
    n1 <- uniroot(function(n1) eval(p.body) - power, c(2, 1e+07))$root
  else if (is.null(n.ratio))
    n.ratio <- uniroot(function(n.ratio) eval(p.body) - power, c(2/n1, 1e+07))$root
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
  n <- c(n1, n1 * n.ratio)
  sd <- c(sd, sd * sd.ratio)

  # Print output as a power.htest object
  structure(list(n = n, delta = delta, sd = sd, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")
}
