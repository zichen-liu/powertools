#' Relative efficiency for cluster randomized or multisite trials
#'
#' @param m The number of subjects per cluster or the mean cluster size (if unequal number of participants per cluster).
#' @param m.sd The standard deviation of cluster sizes (provide if unequal number of participants per cluster).
#' @param icc1 The intraclass correlation coefficient.
#' @return A list of the arguments (including the computed RE).
#' @export
#'
#' @examples
#' pss.re(m = 25, m.sd = 15, icc1 = 0.05)


pss.re <- function (m, m.sd, icc1) {
  cv <- m.sd / m
  K <- (m * icc1) / (1 + (m - 1) * icc1)
  RE <- 1 - cv^2 * K * (1 - K)

  METHOD <- "Relative efficiency for cluster randomized or multisite trials"
  structure(list(m = m, m.sd = m.sd, RE = RE, method = METHOD), class = "power.htest")

}
