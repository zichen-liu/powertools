#' Correlation between a cluster mean at baseline and follow up
#'
#' @description
#' For a cluster randomized trial with a continuous outcome, this function calculates the correlation
#' between a cluster's mean at baseline and at follow up. For cross-sectional sampling of subjects,
#' that is, different subjects are measured at baseline and follow up, specify sac = 0.
#'
#'
#' @param m The number of measurements in each cluster at baseline and follow up.
#' @param icc The intraclass correlation coefficient.
#' @param cac The cluster autocorrelation.
#' @param sac The subject autocorrelation.
#' @return The computed correlation.
#' @export
#'
#' @examples
#' crt.means.r(m = 30, icc = 0.05, cac = 0.4, sac = 0.5)


crt.means.r <- function (m, icc, cac, sac) {
  check(m, "req"); check(m, "pos")
  check(icc, "req"); check(icc, "uniti")
  check(cac, "req"); check(cac, "uniti")
  check(sac, "req"); check(sac, "uniti")

  denom <- 1 + (m - 1) * icc
  num <- m * icc * cac + (1 - icc) * sac
  r <- num / denom

  return(r)
}
