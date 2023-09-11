#' Power calculations for precision analysis for one mean
#'
#' @param n The sample size.
#' @param h The desired halfwidth.
#' @param sd The estimated standard deviation; defaults to 1.
#' @param d The standardized halfwidth. Either d OR h and sd must be specified.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param sides Either 1 or 2 (default) to specify a one- or two- sided hypothesis test.
#' @param cond Specify using unconditional or conditional probability. Defaults to FALSE.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 9.2
#' pss.ci.mean(n = NULL, d = 0.25, power = 0.8)
#' # Example 9.3
#' library(PowerTOST)
#' pss.ci.mean(n = 73, d = 0.25, cond = TRUE)

pss.ci.mean <- function (n = NULL, d = NULL, h = NULL, sd = 1,
                         alpha = 0.05, power = NULL, sides = 2, cond = FALSE) {

  # Check if the arguments are specified correctly
  if (sides != 1 & sides != 2)
    stop("please specify 1 or 2 sides")
  if (sum(sapply(list(n, power, alpha), is.null)) != 1)
    stop("exactly one of n, alpha, and power must be NULL")
  if (is.null(d) & (is.null(h) | is.null(sd)))
    stop("d OR h and sd must be specified")

  if (is.null(d)) d <- h / sd
  p.body <- quote({
    df <- n - 1
    t <- stats::qt(1 - alpha / sides, df)
    b1 <- d * sqrt(n * df) / t
    if (cond) {
      uQ <- PowerTOST::OwensQ(nu = df, t = t, delta = 0, a = 0, b = b1)
      lQ <- PowerTOST::OwensQ(nu = df, t = 0, delta = 0, a = 0, b = b1)
      (2 / (1 - alpha)) * (uQ - lQ)
    } else {
      stats::pchisq(b1^2, df)
    }
  })

  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- paste0("Precision analysis for one mean\n     using ",
                   ifelse(cond, "", "un"), "conditional probability")

  # Print output as a power.htest object
  structure(list(n = n, d = d, alpha = alpha,
                 power = power, sides = sides,
                 method = METHOD), class = "power.htest")

}
