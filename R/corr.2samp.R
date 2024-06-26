#' Power calculations for comparing two correlation coefficients
#'
#' @param n1 The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param rho1 The correlation coefficient in group 1.
#' @param rho2 The correlation coefficient in group 2.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' corr.2samp(n1 = 300, rho1 = 0.3, rho2 = 0.1, sides = 1)

corr.2samp <- function (n1 = NULL, n.ratio = 1, rho1 = NULL, rho2 = NULL,
                        alpha = 0.05, power = NULL, sides = 2,
                        v = FALSE) {

  # Check if the arguments are specified correctly
  if (sides != 1 & sides != 2)
    stop("please specify 1 or 2 sides")
  if (sum(sapply(list(n1, n.ratio, power, alpha), is.null)) != 1)
    stop("exactly one of n1, n.ratio, alpha, and power must be NULL")
  if (is.null(rho1) | is.null(rho2))
    stop("please specify rho1 and rho2")
  if (!is.null(n1) && any(n1 < 4))
    stop("number of observations must be at least 4")
  check(v, "req"); check(v, "bool")

  fisherz <- function(rho){
    0.5 * log((1 + rho)/(1 - rho))
  }

  # Calculate power
  p.body <- quote({
    za <- stats::qnorm(1 - alpha / sides)
    f1 <- fisherz(rho1) + rho1 / (2 * (n1 - 1))
    f2 <- fisherz(rho2) + rho2 / (2 * (n1 * n.ratio - 1))
    DeltaA <- abs(f1 - f2)
    lambda <- DeltaA / sqrt(1 / (n1 - 3) + 1 / (n1 * n.ratio - 3))

    stats::pnorm(lambda - za)
  })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(n1)) {
    n1 <- stats::uniroot(function(n1) eval(p.body) - power, c(4 + 1e-10, 1e+09))$root
    if (!v) return(n1)
  }
  else if (is.null(n.ratio)) {
    n.ratio <- stats::uniroot(function(n.ratio) eval(p.body) - power, c(4/n1, 1e+07))$root
    if (!v) return(n.ratio)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-7, 1 - 1e-7))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Power calculation for comparing two correlation coefficients"
  n <- c(n1, n1 * n.ratio)

  # Print output as a power.htest object
  structure(list(n = n, rho1 = rho1, rho2 = rho2, alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")


}
