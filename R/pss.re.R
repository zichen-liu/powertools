#' Relative efficiency for cluster randomized or multisite trials
#' Relative efficiency of varying to equal cluster sizes
#'
#' @param m The number of subjects per cluster or the mean cluster size (if unequal number of participants per cluster).
#' @param m.sd The standard deviation of cluster sizes (provide if unequal number of participants per cluster).
#' @param icc The intraclass correlation coefficient. For a multisite trial this is icc1. For a CRT this is the average of the 2 icc's.
#' @return A list of the arguments (including the computed RE).
#' @export
#'
#' @examples
#' pss.re(m = 25, m.sd = 15, icc = 0.05)


pss.re <- function (m, m.sd, icc) {
  pss.check(m, "req"); pss.check(m, "pos")
  pss.check(m.sd, "req"); pss.check(m.sd, "min", min = 0)
  pss.check(icc, "req"); pss.check(icc, "uniti")

  cv <- m.sd / m
  K <- (m * icc) / (1 + (m - 1) * icc)
  RE <- 1 - cv^2 * K * (1 - K)
  return(RE)
}
