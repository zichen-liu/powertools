#' Power calculations for one correlation coefficient
#'
#' @param n The sample size.
#' @param r0 The correlation coefficient under the null hypothesis.
#' @param rA The correlation coefficient under the alternative hypothesis.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples#' # Example 10.1
#' pss.corr.1samp(n = 100, rA = 0.2)
#' pss.corr.1samp(n = 100, r0 = 0.2, rA = 0.4)

pss.corr.1samp <- function (n = NULL, r0 = 0, rA = NULL,
                            alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(n, power, alpha), is.null)) != 1)
    stop("exactly one of n, alpha, and power must be NULL")
  if (!is.null(rA))
    stop("please specify rA")
  if (!is.null(n) && any(n < 4))
    stop("number of observations must be at least 4")

  # Calculate power
  p.body <- quote({
    za <- stats::qnorm(1 - alpha)
    K1 <- (5 + rA^2) / (4 * (n - 1))
    K2 <- (11 + 2 * rA^2 + 3 * rA^4) / (8 * (n - 1)^2)
    K <- 1 + K1 + K2

    om1 <- 0.5 * log((1 + rA) / (1 - rA))
    om2 <- K * rA / (2 * (n - 1))
    om3 <- 0.5 * log((1 + r0) / (1 - r0)) + r0 / (2 * (n - 1))
    omega <- sqrt(n - 3) * (om1 + om2 - om3)

    nu1 <- (4 - rA^2) / (2 * (n - 1))
    nu2 <- (22 - 6 * rA^2 - 3 * rA^4) / (6 * (n - 1)^2)
    nu <- ((n - 3) / (n - 1)) * (1 + nu1 + nu2)

    stats::pnorm((omega - za) / sqrt(nu))
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
  METHOD <- "Single correlation coefficient power calculation"

  # Print output as a power.htest object
  structure(list(n = n, r0 = r0, rA = rA, alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")


}
