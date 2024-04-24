#' Power for test of heterogeneity of treatment effect in multisite trials

#' @param m The total number of subjects per site.
#' @param J The number of sites.
#' @param VR The variance ratio (site-level treatment effects / observations within sites) the under the alternative.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' pss.multisite.hte(m = 10, J = 30, VR = 8 / 36)

pss.multisite.hte <- function (m = NULL, alloc.ratio = 1, J = NULL, VR = NULL, alpha = 0.05) {

  # Check if the arguments are specified correctly
  pss.check(m, "req"); pss.check(m, "int")
  pss.check(alloc.ratio, "req"); pss.check(alloc.ratio, "pos")
  pss.check(J, "req"); pss.check(J, "min", min = 2)
  pss.check(VR, "req"); pss.check(VR, "pos")
  pss.check(alpha, "req"); pss.check(alpha, "unit")

  omega <- 1 + m * VR / 4
  df1 <- J - 1
  df2 <- J * (m - 2)
  crit <- stats::qf(1 - alpha, df1, df2)
  power  <- 1 - stats::pf(crit / omega, df1, df2)

  # Generate output text
  METHOD <-"Power for test of heterogeneity of treatment effect in multisite trials"
  NOTE <- "m1, m2 are the number of subjects within site in condition 1, condition 2\n      (total of m1 + m2 per site)"
  c <- m / (alloc.ratio + 1)
  t <- alloc.ratio * c
  m <- paste0(t, ", ", c)

  # Print output as a power.htest object depending on which inputs were given
  structure(list(`m1, m2` = m, J = J, VR = VR,
                 alpha = alpha, power = power,
                 method = METHOD, note = NOTE), class = "power.htest")

}
