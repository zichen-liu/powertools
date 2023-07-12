#' Power calculations for one and two sample z tests with unequal sample size
#'
#' @param n The sample size (in the smallest group if two groups).
#' @param delta For a two-sample test, the true difference DeltaA-Delta0. For a one-sample test, the true difference meanA-mean0.
#' @param sigma The estimated standard deviation. Defaults to 1. For a paired test, sigma is the standard deviation of *differences* within pairs.
#' @param alpha The significance level or type 1 error rate.
#' @param power The specified level of power.
#' @param n.ratio The ratio n2/n1 between the larger group and the smaller group. Should be a value equal to or greater than 1 since n2 is the larger group. Defaults to 1 (equal group sizes).
#' @param sd.ratio The ratio sd2/sd1 between the standard deviations in the larger group and the smaller group. Defaults to 1 (equal standard deviations in the two groups).
#' @param type Type of z test: "one.sample", "two.sample" (default), or "paired".
#' @param one.or.two.sided Either "one" or "two" (default) to specify a one- or two- sided hypothesis test.
#' @param strict Use strict interpretation in two-sided case. Defaults to TRUE.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples pss.z.test(delta = 6.3-5.7, sigma = 2, alpha = 0.05, power = 0.8, type = "one.sample", one.or.two.sided = "one")
#' pss.z.test(n = 40, delta = 2, sigma = 5, sd.ratio = 2, n.ratio = 1.5, alpha = 0.05, type = "two.sample", one.or.two.sided = "two")
#'
pss.z.test <- function(n = NULL, delta = NULL, sigma = 1,
                         alpha = 0.05, power = NULL, n.ratio = 1, sd.ratio = 1,
                         type = c("two.sample", "one.sample", "paired"),
                         one.or.two.sided = c("two", "one"), strict = TRUE) {

  # Check if the arguments are specified correctly
  type <- match.arg(type)
  if (type == "two.sample") {
    if (sum(sapply(list(n, delta, sd, power, alpha, n.ratio, sd.ratio), is.null)) != 1)
      stop("exactly one of n, d, power, alpha, n.ratio and sd.ratio must be NULL")
    if (!is.null(n.ratio) && n.ratio < 1)
      stop("n.ratio between group sizes cannot be less than 1")
    if (!is.null(sd.ratio) && sd.ratio < 1)
      stop("sd.ratio between group sd's cannot be less than 1")}
  else {
    n.ratio <- 1
    sd.ratio <- 1
    if (sum(sapply(list(n, delta, sd, power, alpha), is.null)) != 1)
      stop("exactly one of n, d, power, and alpha must be NULL")}

  # Assign number of samples and sides
  one.or.two.sided <- match.arg(one.or.two.sided)
  sample <- switch(type, one.sample = 1, two.sample = 2, paired = 1)
  side <- switch(one.or.two.sided, one = 1, two = 2)

  # Use absolute value of the effect size
  if (!is.null(delta))
    delta <- abs(delta)

  # For 1 sample, power = z + delta / (sigma/sqrt(n))
  # For 2 sample, power = z + delta / sqrt(s1^2/n1 + s2^2/n2)
  p.body <- quote({
    sd <- switch(sample, sigma/sqrt(n),
                 sqrt((sigma * sd.ratio)^2 / (n * n.ratio) + sigma^2/n))
    stats::pnorm(stats::qnorm(alpha/side) + delta/sd)})
  if (strict & side == 2)
    p.body <- quote({
      sd <- switch(sample, sigma/sqrt(n),
                   sqrt((sigma * sd.ratio)^2 / (n * n.ratio) + sigma^2/n))
      stats::pnorm(stats::qnorm(alpha/side) + delta/sd) +
      stats::pnorm(stats::qnorm(alpha/side) - delta/sd)
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

  method <- paste(switch(type, one.sample = "One-sample z test power calculation",
            two.sample = ifelse(n.ratio == 1, "Two-sample z test power calculation",
            "Two-sample z test power calculation with unequal sample sizes"),
            paired = "Paired z test power calculation"))
  if (type == "two.sample" & sd.ratio != 1) {
    method <- paste0(method, ifelse(n.ratio == 1, " with", " and"), " unequal variances")
  }

  # Print output as a power.htest object
  structure(list(n = n, delta = delta, sigma = sigma, alpha = alpha,
                 power = power, one.or.two.sided = one.or.two.sided,
                 method = method, note = note), class = "power.htest")
}
