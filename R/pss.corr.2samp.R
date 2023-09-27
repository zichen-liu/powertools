#' Power calculations for comparing two correlation coefficients
#'
#' @param n The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param rho1 The correlation coefficient in the first group.
#' @param rho2 The correlation coefficient in the second group.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed one).
#' @import psych
#' @export
#'
#' @examples#' # Example 10.2
#' pss.corr.2samp(n = 300, rho1 = 0.3, rho2 = 0.1)

pss.corr.2samp <- function (n = NULL, n.ratio = 1, rho1 = NULL, rho2 = NULL,
                            alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(n, n.ratio, power, alpha), is.null)) != 1)
    stop("exactly one of n, n.ratio, alpha, and power must be NULL")
  if (is.null(rho1) | is.null(rho2))
    stop("please specify rho1 and rho2")
  if (!is.null(n) && any(n < 4))
    stop("number of observations must be at least 4")

  # Calculate power
  p.body <- quote({
    za <- stats::qnorm(1 - alpha)
    f1 <- psych::fisherz(rho1) + rho1 / (2 * (n - 1))
    f2 <- psych::fisherz(rho2) + rho2 / (2 * (n * n.ratio - 1))
    DeltaA <- abs(f1 - f2)
    lambda <- DeltaA / sqrt(1 / (n - 3) + 1 / (n * n.ratio - 3))

    stats::pnorm(lambda - za)
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(4 + 1e-10, 1e+09))$root
  else if (is.null(n.ratio))
    n.ratio <- uniroot(function(ratio) eval(p.body) - power, c(2/n, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "Power calculation for comparing two correlation coefficients"
  n <- c(n, n * n.ratio)

  # Print output as a power.htest object
  structure(list(n = n, rho1 = rho1, rho2 = rho2, alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")


}
