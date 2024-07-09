#' Power calculations for precision analysis for one mean
#'
#' @param N The sample size.
#' @param halfwidth The desired halfwidth.
#' @param sd The estimated standard deviation; defaults to 1.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param cond Specify using unconditional or conditional probability. Defaults to FALSE.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @import PowerTOST
#' @export
#'
#' @examples
#' ci.mean(N = NULL, halfwidth = 0.25, power = 0.8)
#' ci.mean(N = 62, halfwidth = 0.25, power = NULL)
#' ci.mean(N = 73, halfwidth = 0.25, cond = TRUE)

ci.mean <- function (N = NULL, halfwidth = NULL, sd = 1,
                     alpha = 0.05, power = NULL, cond = FALSE,
                     v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(N, alpha, power), "oneof")
  check(N, "pos")
  check(alpha, "unit")
  check(power, "unit")
  check(halfwidth, "req"); check(halfwidth, "num")
  check(sd, "req"); check(sd, "pos")
  check(cond, "req"); check(cond, "bool")
  check(v, "req"); check(v, "bool")

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

  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(N)) {
    N <- stats::uniroot(function(N) eval(p.body) - power, c(2, 1e+07))$root
    if (!v) return(N)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- paste0("Precision analysis for one mean\n     using ",
                   ifelse(cond, "", "un"), "conditional probability")

  # Print output as a power.htest object
  structure(list(N = N, halfwidth = halfwidth, sd = sd,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
