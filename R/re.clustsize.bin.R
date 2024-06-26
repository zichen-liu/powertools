#' Relative efficiency for cluster randomized trials with binary outcomes
#' due to varying cluster sizes
#'
#' @param m The number of subjects per cluster or the mean cluster size (if unequal number of participants per cluster).
#' @param m.sd The standard deviation of cluster sizes (provide if unequal number of participants per cluster).
#' @param pc The probability of the outcome in control clusters.
#' @param pt The probability of the outcome in treatment clusters.
#' @param sigma.u Standard deviation of the cluster random effect.
#' @return The computed RE.
#' @export
#'
#' @examples
#' re.clustsize.bin(m = 60, m.sd = 45, pc = 0.25, pt = 0.15, sigma.u = 0.3)

re.clustsize.bin <- function(m, m.sd, pc, pt, sigma.u){
  check(m, "req"); check(m, "pos")
  check(m.sd, "req"); check(m.sd, "min", min = 0)
  check(pc, "req"); check(pc, "unit")
  check(pt, "req"); check(pt, "unit")
  check(sigma.u, "req"); check(sigma.u, "pos")

  logoddsc <- log(pc / (1 - pc))
  logoddst <- log(pt / (1 - pt))
  gam0 <- (logoddsc + logoddst) / 2
  gam1 <- logoddst - logoddsc

  sigmasqet <- 2 + exp(-gam0 - 2 * gam1) + exp(gam0 + 2 * gam1)
  sigmasqec <- 2 + exp(-gam0 + 2 * gam1) + exp(gam0 - 2 * gam1)

  alpha_t <- sigmasqet / sigma.u^2
  alpha_c <- sigmasqec / sigma.u^2

  cv <- m.sd / m

  lamt <- m / (m + alpha_t)
  lamc <- m / (m + alpha_c)

  re.num <- (1 - cv^2 * lamt * (1 - lamt)) * (1 - cv^2 * lamc * (1 - lamc)) * (lamt + lamc)
  re.denom <- lamt + lamc - cv^2 * (lamt^2 * (1 - lamt) + lamc^2 * (1 - lamt))
  re <- re.num / re.denom

  return(re)
}


