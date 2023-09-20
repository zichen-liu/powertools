#' Power calculations for precision analysis for a difference between means
#'
#' @param n The sample size for group 1.
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param h The desired halfwidth for the difference in means.
#' @param sd The estimated standard deviation; defaults to 1.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param cond Specify using unconditional or conditional probability. Defaults to FALSE.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 9.4
#' pss.ci.meandiff(n = NULL, h = 0.25, power = 0.8)
#' # Example 9.5
#' library(PowerTOST)
#' pss.ci.meandiff(n = 134, h = 0.25, cond = TRUE)

pss.ci.meandiff <- function (n = NULL, n.ratio = 1, h = NULL, sd = 1,
                              alpha = 0.05, power = NULL, cond = FALSE) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(n, n.ratio, power, alpha), is.null)) != 1)
    stop("exactly one of n, n.ratio, alpha, and power must be NULL")
  if (is.null(h))
    stop("h must be specified")

  d <- h / sd
  p.body <- quote({
    df <- n * (n.ratio + 1) - 2
    t <- stats::qt(1 - alpha / 2, df)
    b <- d * sqrt(n * n.ratio * df) / (t * sqrt(n.ratio + 1))
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
  else if (is.null(n.ratio))
    n.ratio <- uniroot(function(ratio) eval(p.body) - power, c(2/n, 1e+07))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- paste0("Precision analysis for difference between means\n     using ",
                   ifelse(cond, "", "un"), "conditional probability")
  n <- c(n, n * n.ratio)

  # Print output as a power.htest object
  structure(list(n = n, h = h, sd = sd,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
