#' Power calculations for one and two sample z tests with unequal sample size
#'
#' @param n The sample size (in the smallest group if two groups).
#' @param d The standardized effect size (delta/sigma). Either d OR delta and sigma need specified.
#' @param delta The true difference between means. If delta is specified, sigma needs to be specified too.
#' @param sigma The estimated standard deviation. If sigma is specified, delta needs to be specified too. For a paired test, sigma is the std.dev. of *differences* within pairs.
#' @param alpha The significance level or type 1 error rate.
#' @param power The specified level of power.
#' @param n.ratio The ratio n2/n1 between the larger group and the smaller group. Should be a value equal to or greater than 1 since n2 is the larger group. Defaults to 1 (equal group sizes).
#' @param sd.ratio The ratio sd2/sd1 between the standard deviations in the larger group and the smaller group. Defaults to 1 (equal standard deviations in the two groups).
#' @param type Type of z test ("one.sample", "two.sample", or "paired)
#' @param one.or.two.sided Either "one" or "two" to specify a one or two sided hypothesis test. Default is two-sided.
#' @param strict Use strict interpretation in two-sided case. Defaults to TRUE.
#'
#' @return A list of the arguments (including the computed one)
#' @export
#'
#' @examples power_z_test(d=300, power=.8, sd.ratio=2)
#'
power_z_test <- function(n = NULL, d = NULL, delta = NULL, sigma = 1,
                         alpha = 0.05, power = NULL, n.ratio = 1, sd.ratio = 1,
          type = c("two.sample", "one.sample", "paired"),
          one.or.two.sided = c("two", "one"), strict = TRUE){

  # Calculate d based on delta/sigma if necessary
  if(is.null(d) & is.null(delta)){
    print(paste0("Either d (delta/sigma) OR delta and sigma need to be specified"))
    stop()}
  if(!is.null(d) & !is.null(delta)){
    print(paste0("Either d (delta/sigma) OR delta and sigma need to be specified"))
    stop()}
  if(!is.null(delta) & !is.null(sigma)){
    d <- delta / sigma}

  # Check if the arguments are specified correctly
  type <- match.arg(type)
  if (type == "two.sample") {
    if (sum(sapply(list(n, d, sd, power, alpha, n.ratio, sd.ratio), is.null)) != 1)
      stop("exactly one of n, d, power, alpha, n.ratio and sd.ratio must be NULL")
    if (!is.null(n.ratio) && n.ratio < 1)
      stop("n.ratio between group sizes cannot be less than 1")
    if (!is.null(sd.ratio) && sd.ratio < 1)
      stop("sd.ratio between group sd's cannot be less than 1")}
  else {
    n.ratio <- 1
    sd.ratio <- 1
    if (sum(sapply(list(n, d, sd, power, alpha), is.null)) != 1)
      stop("exactly one of n, d, power, and alpha must be NULL")}

  # Assign number of samples and sides
  one.or.two.sided <- match.arg(one.or.two.sided)
  sample <- switch(type, one.sample = 1, two.sample = 2, paired = 1)
  side <- switch(one.or.two.sided, one = 1, two = 2)

  # Absolute value of effect size
  if (!is.null(d))
    d <- abs(d)

  # For 1 sample, power = z + d / (1/sqrt(n))
  # For 2 sample, power = z + d / sqrt(1/n + s2^2/n2)
  p.body <- quote({
    sd <- switch(sample, 1/sqrt(n), sqrt((sd.ratio)^2/(n * n.ratio) + 1/n))
    stats::pnorm(stats::qnorm(alpha/side) + d/sd, lower.tail = F)})
  if (strict & side == 2)
    p.body <- quote({
      sd <- switch(sample, 1/sqrt(n), sqrt((sd.ratio)^2/(n * n.ratio) + 1/n))
      stats::pnorm(stats::qnorm(alpha/side) + d/sigma, lower.tail = F) +
      stats::pnorm(stats::qnorm(alpha/side) - d/sigma, lower.tail = F)
    })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+07))$root
  else if (is.null(d))
    d <- uniroot(function(d) eval(p.body) - power,  c(1e-07, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else if (is.null(n.ratio))
    n.ratio <- uniroot(function(n.ratio) eval(p.body) - power, c(2/n, 1e+07))$root
  else if (is.null(sd.ratio))
    sd.ratio <- uniroot(function(sd.ratio) eval(p.body) - power, c(1e-07, 1e+07))$root
  else stop("internal error")

  if (type == "two.sample" & (n.ratio != 1 | sd.ratio != 1)) {
    n <- c(n, n * n.ratio)
    sigma <- c(1, sd.ratio)
  }

  return(data.frame(n = n, d = d, sigma = sigma, alpha = alpha,
              one.or.two.sided = one.or.two.sided))
}
