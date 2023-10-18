#' Relative efficiency due to unequal number of participants per site
#'
#' @param m.mean The mean cluster size.
#' @param m.sd The standard deviation of cluster sizes.
#' @param rho The intraclass correlation among cluster sizes.
#' @return A list of the arguments (including the computed relative efficiency).
#' @export
#'
#' @examples
#' 1 / pss.multisite.re(m.mean = 30, m.sd = 23, rho = 0.05)$re
#' 1 / pss.multisite.re(m.mean = 30, m.sd = 23, rho = 0.1)$re
#'

pss.multisite.re <- function (m.mean = NULL, m.sd = NULL, rho = NULL) {
 cv <- m.sd / m.mean
 K <- (m.mean * rho) / (1 + (m.mean - 1) * rho)
 re <- 1 - cv^2 * K * (1 - K)

 METHOD <- "Relative efficiency due to unequal number of participants per site"

 structure(list(m.mean = m.mean, m.sd = m.sd, rho = rho,
                re = re, method = METHOD), class = "power.htest")
}
