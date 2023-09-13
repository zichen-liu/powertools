#' Power calculations for comparing two correlation coefficients
#'
#' @param n The sample size in each group.
#' @param rho1 The correlation coefficient in the first group.
#' @param rho2 The correlation coefficient in the second group.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples#' # Example 10.2
#' library(psych)
#' pss.corr.2samp(n = 300, rho1 = 0.3, rho2 = 0.1)

pss.corr.2samp <- function (n = NULL, rho1 = NULL, rho2 = NULL,
                            alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(n, power, alpha), is.null)) != 1)
    stop("exactly one of n, alpha, and power must be NULL")
  if (is.null(rho1) | is.null(rho2))
    stop("please specify rho1 and rho2")
  if (!is.null(n) && any(n < 4))
    stop("number of observations must be at least 4")

  # Calculate power
  p.body <- quote({
    za <- stats::qnorm(1 - alpha)
    f1 <- psych::fisherz(rho1) + rho1 / (2 * (n - 1))
    f2 <- psych::fisherz(rho2) + rho2 / (2 * (n - 1))
    DeltaA <- abs(f1 - f2)
    lambda <- DeltaA / sqrt(1 / (n - 3) + 1 / (n - 3))

    stats::pnorm(lambda - za)
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(4 + 1e-10, 1e+09))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "Power calculation for comparing two correlation coefficients"
  NOTE <- "n is the number in each group"

  # Print output as a power.htest object
  structure(list(n = n, rho1 = rho1, rho2 = rho2, alpha = alpha, power = power,
                 method = METHOD, note = NOTE), class = "power.htest")


}
