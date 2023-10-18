#' Power for test of heterogeneity of treatment effect

#' @param m The number of subjects per site.
#' @param J The number of sites.
#' @param ssq.u1 The variance of the site-level treatment effects (sigma squared u1) under the alternative.
#' @param ssq.e The variance of the observations within sites (sigma squared epsilon) under the alternative.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' pss.multisite.hte(m = 10, J = 30, ssq.u1 = 8, ssq.e = 36)

pss.multisite.hte <- function (m = NULL, J = NULL, ssq.u1 = NULL, ssq.e = NULL,
                               alpha = 0.05) {
  VR <- ssq.u1 / ssq.e
  omega <- 1 + m * VR / 4
  df1 <- J - 1
  df2 <- J * (m - 2)
  crit <- stats::qf(1 - alpha, df1, df2)
  power  <- 1 - stats::pf(crit / omega, df1, df2)

  # Generate output text
  METHOD <-"Power for test of heterogeneity of treatment effect"

  # Print output as a power.htest object depending on which inputs were given
  structure(list(m = m, J = J, ssq.u1 = ssq.u1, ssq.e = ssq.e,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
