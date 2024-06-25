#' Relative efficiency for cluster randomized or multisite trials
#' due to varying cluster sizes
#'
#' @param m The number of subjects per cluster or the mean cluster size (if unequal number of participants per cluster).
#' @param m.sd The standard deviation of cluster sizes (provide if unequal number of participants per cluster).
#' @param icc The intraclass correlation coefficient. For a multisite trial this is icc1. For a CRT this is the average of the 2 icc's.
#' @return The computed RE.
#' @export
#'
#' @examples
#' re.clustsize.cont(m = 25, m.sd = 15, icc = 0.05)


re.clustsize.cont <- function (m, m.sd, icc) {
  check(m, "req"); check(m, "pos")
  check(m.sd, "req"); check(m.sd, "min", min = 0)
  check(icc, "req"); check(icc, "uniti")

  cv <- m.sd / m
  K <- (m * icc) / (1 + (m - 1) * icc)
  RE <- 1 - cv^2 * K * (1 - K)
  return(RE)
}
