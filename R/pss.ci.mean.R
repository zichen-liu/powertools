#' Power calculations for precision analysis for one mean
#'
#' @param n The sample size.
#' @param h The desired halfwidth.
#' @param sd The estimated standard deviation; defaults to 1.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param cond Specify using unconditional or conditional probability. Defaults to FALSE.
#'
#' @return A list of the arguments (including the computed one).
#' @import PowerTOST
#' @export
#'
#' @examples
#' # Example 9.2
#' pss.ci.mean(n = NULL, h = 0.25, power = 0.8)
#' pss.ci.mean(n = 62, h = 0.25, power = NULL)
#' # Example 9.3
#' pss.ci.mean(n = 73, h = 0.25, cond = TRUE)

pss.ci.mean <- function (n = NULL, h = NULL, sd = 1,
                         alpha = 0.05, power = NULL, cond = FALSE) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(n, power, alpha), is.null)) != 1)
    stop("exactly one of n, alpha, and power must be NULL")
  if (is.null(h))
    stop("h must be specified")
  d <- h / sd

  p.body <- quote({
    df <- n - 1
    t <- stats::qt(1 - alpha / 2, df)
    b <- d * sqrt(n * df) / t
    if (cond) {
      uQ <- PowerTOST::OwensQ(nu = df, t = t, delta = 0, a = 0, b = b)
      lQ <- PowerTOST::OwensQ(nu = df, t = 0, delta = 0, a = 0, b = b)
      (2 / (1 - alpha)) * (uQ - lQ)
    } else {
      stats::pchisq(b^2, df)
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
  structure(list(n = n, h = h, sd = sd,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
