#' Power calculations for one correlation coefficient
#'
#' @param n The sample size.
#' @param rho0 The correlation coefficient under the null hypothesis; defaults to 0.
#' @param rhoA The correlation coefficient under the alternative hypothesis.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples#' # Example 10.1
#' pss.corr.1samp(n = 100, rhoA = 0.2)
#' pss.corr.1samp(n = 100, rho0 = 0.2, rhoA = 0.4)

pss.corr.1samp <- function (n = NULL, rho0 = 0, rhoA = NULL,
                            alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(n, power, alpha), is.null)) != 1)
    stop("exactly one of n, alpha, and power must be NULL")
  if (is.null(rhoA))
    stop("please specify rhoA")
  if (!is.null(n) && any(n < 4))
    stop("number of observations must be at least 4")

  # Calculate power
  p.body <- quote({
    za <- stats::qnorm(1 - alpha)
    K1 <- (5 + rhoA^2) / (4 * (n - 1))
    K2 <- (11 + 2 * rhoA^2 + 3 * rhoA^4) / (8 * (n - 1)^2)
    K <- 1 + K1 + K2

    om1 <- 0.5 * log((1 + rhoA) / (1 - rhoA))
    om2 <- K * rhoA / (2 * (n - 1))
    om3 <- 0.5 * log((1 + rho0) / (1 - rho0)) + rho0 / (2 * (n - 1))
    omega <- sqrt(n - 3) * (om1 + om2 - om3)

    nu1 <- (4 - rhoA^2) / (2 * (n - 1))
    nu2 <- (22 - 6 * rhoA^2 - 3 * rhoA^4) / (6 * (n - 1)^2)
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
  structure(list(n = n, rho0 = rho0, rhoA = rhoA, alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")


}
