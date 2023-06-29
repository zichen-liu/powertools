#' Power calculations for one and two sample z tests with unequal sample size
#'
#' @param n Number of observations (in the smallest group if two groups)
#' @param delta True difference in means
#' @param sd Standard deviation
#' @param sig.level Significance level (Type I error probability)
#' @param power Power of test (1 minus Type II error probability)
#' @param ratio The ratio n2/n1 between the larger group and the smaller group. Should be a value equal to or greater than 1 since n2 is the larger group. Defaults to 1 (equal group sizes). If ratio is set to NULL (i.e., find the ratio) then the ratio might be smaller than 1 depending on the desired power and ratio of the sd's.
#' @param sd.ratio The ratio sd2/sd1 between the standard deviations in the larger group and the smaller group. Defaults to 1 (equal standard deviations in the two groups)
#' @param type Type of t test
#' @param alternative One- or two-sided test
#' @param strict Use strict interpretation in two-sided case. Defaults to TRUE unlike the standard power.t.test function.
#'
#' @return
#' @export
#'
#' @examples power_t_test(delta=300, sd=450, power=.8, sd.ratio=2)
#'
function (n = NULL, delta = NULL, sd = 1, sig.level = 0.05,
          power = NULL, ratio = 1, sd.ratio = 1,
          type = c("two.sample", "one.sample", "paired"),
          alternative = c("two.sided", "one.sided"), strict = TRUE)
{
  type <- match.arg(type)
  if (type == "two.sample") {
    if (sum(sapply(list(n, delta, sd, power, sig.level,
                        ratio, sd.ratio), is.null)) != 1)
      stop("exactly one of n, delta, sd, power, sig.level, ratio and sd.ratio must be NULL")
    if (!is.null(ratio) && ratio < 1)
      stop("ratio between group sizes cannot be less than 1")
    if (!is.null(sd.ratio) && sd.ratio < 1)
      stop("sd.ratio between group sd's cannot be less than 1")
  }
  else {
    ratio <- 1
    sd.ratio <- 1
    if (sum(sapply(list(n, delta, sd, power, sig.level),
                   is.null)) != 1)
      stop("exactly one of n, delta, sd, power, and sig.level must be NULL")
  }
  alternative <- match.arg(alternative)
  tsample <- switch(type, one.sample = 1, two.sample = 2, paired = 1)
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  if (!is.null(delta))
    delta <- abs(delta)
  p.body <- quote({
    sigma <- switch(tsample, sd / sqrt(n), sqrt((sd * sd.ratio)^2/(n * ratio) + sd^2/n))
    pnorm(qnorm(sig.level/tside) + delta/sigma, lower.tail = F)
  })
  if (strict & tside == 2)
    p.body <- quote({
      sigma <- switch(tsample, sd / sqrt(n), sqrt((sd * sd.ratio)^2/(n * ratio) + sd^2/n))
      pnorm(qnorm(sig.level/tside) + delta/sigma, lower.tail = F) +
      pnorm(qnorm(sig.level/tside) - delta/sigma, lower.tail = F)
    })
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+07))$root
  else if (is.null(sd))
    sd <- uniroot(function(sd) eval(p.body) - power, delta * c(1e-07, 1e+07))$root
  else if (is.null(delta))
    delta <- uniroot(function(delta) eval(p.body) - power, sd * c(1e-07, 1e+07))$root
  else if (is.null(sig.level))
    sig.level <- uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else if (is.null(ratio))
    ratio <- uniroot(function(ratio) eval(p.body) - power, c(2/n, 1e+07))$root
  else if (is.null(sd.ratio))
    sd.ratio <- uniroot(function(sd.ratio) eval(p.body) - power, c(1e-07, 1e+07))$root
  else stop("internal error")
  NOTE <- switch(type, paired = "n is number of *pairs*, sd is std.dev. of *differences* within pairs",
                 two.sample = ifelse(ratio == 1, "n is number in *each* group",
                                     "n is vector of number in each group"), NULL)
  if (type == "two.sample" & (ratio != 1 | sd.ratio != 1)) {
    n <- c(n, n * ratio)
    sd <- c(sd, sd * sd.ratio)
  }
  METHOD <- paste(switch(type, one.sample = "One-sample z test power calculation",
                         two.sample = ifelse(ratio == 1, "Two-sample z test power calculation",
                                             "Two-sample z test power calculation with unequal sample sizes"),
                         paired = "Paired z test power calculation"))
  if (type == "two.sample" & sd.ratio != 1) {
    METHOD <- paste0(METHOD, ifelse(ratio == 1, " with",
                                    " and"), " unequal variances")
  }
  structure(list(n = n, delta = delta, sd = sd, sig.level = sig.level,
                 power = power, alternative = alternative, note = NOTE,
                 method = METHOD), class = "power.htest")
}
