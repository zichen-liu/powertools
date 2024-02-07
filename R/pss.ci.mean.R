#' Power calculations for precision analysis for one mean
#'
#' @param N The sample size.
#' @param halfwidth The desired halfwidth.
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
#' pss.ci.mean(N = NULL, halfwidth = 0.25, power = 0.8)
#' pss.ci.mean(N = 62, halfwidth = 0.25, power = NULL)
#' pss.ci.mean(N = 73, halfwidth = 0.25, cond = TRUE)

pss.ci.mean <- function (N = NULL, halfwidth = NULL, sd = 1,
                         alpha = 0.05, power = NULL, cond = FALSE) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(N, power, alpha), is.null)) != 1)
    stop("exactly one of N, alpha, and power must be NULL")
  if (is.null(halfwidth))
    stop("halfwidth must be specified")
  d <- halfwidth / sd

  p.body <- quote({
    df <- N - 1
    t <- stats::qt(1 - alpha / 2, df)
    b <- d * sqrt(N * df) / t
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
  else if (is.null(N))
    N <- stats::uniroot(function(N) eval(p.body) - power, c(2, 1e+07))$root
  else if (is.null(alpha))
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- paste0("Precision analysis for one mean\n     using ",
                   ifelse(cond, "", "un"), "conditional probability")

  # Print output as a power.htest object
  structure(list(N = N, halfwidth = halfwidth, sd = sd,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
