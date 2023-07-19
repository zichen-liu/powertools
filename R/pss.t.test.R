#' Power calculations for one and two sample t tests with unequal sample size
#'
#' @param n The sample size (in the smallest group if two groups).
#' @param delta For a two-sample test, the true difference DeltaA-Delta0. For a one-sample test, the true difference meanA-mean0.
#' @param sigma The estimated standard deviation. Defaults to 1. For a paired test, sigma is the standard deviation of *differences* within pairs.
#' @param alpha The significance level or type 1 error rate.
#' @param power The specified level of power.
#' @param n.ratio The ratio n2/n1 between the larger group and the smaller group. Should be a value equal to or greater than 1 since n2 is the larger group. Defaults to 1 (equal group sizes).
#' @param sd.ratio The ratio sd2/sd1 between the standard deviations in the larger group and the smaller group. Defaults to 1 (equal standard deviations in the two groups).
#' @param type Type of t test: "one.sample", "two.sample" (default), or "paired".
#' @param one.or.two.sided Either "one" or "two" (default) to specify a one- or two- sided hypothesis test.
#' @param df.method Method for calculating the degrees of freedom: "welch" (default) or "classical".
#' @param strict Use strict interpretation in two-sided case. Defaults to TRUE.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples pss.t.test(n = 50, delta = 0.5, sigma = 1, n.ratio = 1.4, power = NULL, type = "two.sample", one.or.two.sided = "two")
#' pss.t.test(n = NULL, delta = 4, sigma = 10.95445, alpha = 0.05, power = 0.8, type = "paired", one.or.two.sided = "two")
#'
pss.t.test <- function(n = NULL, delta = NULL, sigma = 1,
                       alpha = 0.05, power = NULL, n.ratio = 1, sd.ratio = 1,
                       type = c("two.sample", "one.sample", "paired"),
                       one.or.two.sided = c("two", "one"),
                       df.method = c("welch", "classical"), strict = TRUE) {

  # Check if the arguments are specified correctly
  type <- match.arg(type)
  if (type == "two.sample") {
    if (sum(sapply(list(n, delta, sigma, power, alpha, n.ratio, sd.ratio), is.null)) != 1)
      stop("exactly one of n, delta, sigma, power, alpha, n.ratio and sd.ratio must be NULL")
    if (!is.null(n.ratio) && n.ratio < 1)
      stop("n.ratio between group sizes cannot be less than 1")
    if (!is.null(sd.ratio) && sd.ratio < 1)
      stop("sd.ratio between group sd's cannot be less than 1")}
  else {
    n.ratio <- 1
    sd.ratio <- 1
    if (sum(sapply(list(n, delta, sigma, power, alpha), is.null)) != 1)
      stop("exactly one of n, delta, sigma, power, and alpha must be NULL")}

  # Assign number of samples and sides
  one.or.two.sided <- match.arg(one.or.two.sided)
  df.method <- match.arg(df.method)
  sample <- switch(type, one.sample = 1, two.sample = 2, paired = 1)
  side <- switch(one.or.two.sided, one = 1, two = 2)

  # Use absolute value of the effect size
  if (!is.null(delta))
    delta <- abs(delta)

  # Copy df calculations from MESS::power_t_test
  p.body <- quote({
    nu <- switch(sample, n - 1, switch(df.method,
          welch = (sigma^2 / n + (sigma * sd.ratio)^2 / (n * n.ratio))^2 /
          ((sigma^2 / n)^2 / (n - 1) + ((sigma * sd.ratio)^2 / (n.ratio * n))^2 /
          (n * n.ratio - 1)),
          classical = (1 + n.ratio) * n - 2))
    stats::pt(stats::qt(alpha / side, nu, lower = FALSE), nu, ncp =
       switch(sample, sqrt(n / sample), sqrt(n / (1 + sd.ratio^2 / n.ratio))) *
       delta / sigma, lower = FALSE)
  })

  if (strict & side == 2)
    p.body <- quote({
      nu <- switch(sample, n - 1, switch(df.method,
            welch = (sigma^2/n + (sigma * sd.ratio)^2 / (n * n.ratio))^2 /
            ((sigma^2 / n)^2/(n - 1) + ((sigma * sd.ratio)^2 / (n.ratio * n))^2 /
            (n * n.ratio - 1)),
            classical = (1 + n.ratio) * n - 2))
      qu <- stats::qt(alpha / side, nu, lower = FALSE)
      ncp <- switch(sample, sqrt(n / sample), sqrt(n / (1 + sd.ratio^2 / n.ratio))) *
             delta / sigma
      stats::pt(qu, nu, ncp = ncp, lower = FALSE) + pt(-qu, nu, ncp = ncp, lower = TRUE)
    })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body)-power, c(2, 1e+07))$root
  else if (is.null(sigma))
    sigma <- uniroot(function(sigma) eval(p.body)-power, delta*c(1e-07, 1e+07))$root
  else if (is.null(delta))
    delta <- uniroot(function(delta) eval(p.body)-power,  c(1e-07, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body)-power, c(1e-10, 1 - 1e-10))$root
  else if (is.null(n.ratio))
    n.ratio <- uniroot(function(n.ratio) eval(p.body)-power, c(2/n, 1e+07))$root
  else if (is.null(sd.ratio))
    sd.ratio <- uniroot(function(sd.ratio) eval(p.body)-power, c(1e-07, 1e+07))$root
  else stop("internal error")

  # For unequal allocation or sd, get both n and/or sd
  if (type == "two.sample" & (n.ratio != 1 | sd.ratio != 1)) {
    n <- c(n, n * n.ratio)
    sigma <- c(sigma, sigma * sd.ratio)
  }

  # Generate output text
  note <- switch(type,
                 paired = "n is the number of *pairs*; sigma is standard deviation of *differences* within pairs",
                 two.sample = "n is the number in *each* group", NULL)

  method <- paste(switch(type, one.sample = "One-sample t test power calculation",
                         two.sample = ifelse(n.ratio == 1, "Two-sample t test power calculation",
                                             "Two-sample t test power calculation with unequal sample sizes"),
                         paired = "Paired t test power calculation"))
  if (type == "two.sample" & sd.ratio != 1) {
    method <- paste0(method, ifelse(n.ratio == 1, " with", " and"), " unequal variances")
  }

  # Print output as a power.htest object
  structure(list(n = n, delta = delta, sigma = sigma, alpha = alpha,
                 power = power, one.or.two.sided = one.or.two.sided,
                 method = method, note = note), class = "power.htest")
}
