#' Correlation between a cluster mean at baseline and follow up
#'
#' @param m The number of measurements in each cluster at baseline and follow up.
#' @param icc The intraclass correlation coefficient.
#' @param cac The cluster autocorrelation.
#' @param sac The subject autocorrelation.
#' @return The computed correlation.
#' @export
#'
#' @examples
#' pss.crt.r(m = 30, icc = 0.05, cac = 0.4, sac = 0.5)


pss.crt.r <- function (m, icc, cac, sac) {
  pss.check(m, "req"); pss.check(m, "pos")
  pss.check(icc, "req"); pss.check(icc, "uniti")
  pss.check(cac, "req"); pss.check(cac, "uniti")
  pss.check(sac, "req"); pss.check(sac, "uniti")

  denom <- 1 + (m - 1) * icc
  num <- m * icc * cac + (1 - icc) * sac
  r <- num / denom

  return(r)
}
