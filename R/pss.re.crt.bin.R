#' Relative efficiency for cluster randomized trials with binary outcomes
#' Relative efficiency of varying to equal cluster sizes
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
#' pss.re.crt.bin(m = 60, m.sd = 45, pc = 0.25, pt = 0.15, sigma.u = 0.3)

pss.re.crt.bin <- function(m, m.sd, pc, pt, sigma.u){
  pss.check(m, "req"); pss.check(m, "pos")
  pss.check(m.sd, "req"); pss.check(m.sd, "min", min = 0)
  pss.check(pc, "req"); pss.check(pc, "unit")
  pss.check(pt, "req"); pss.check(pt, "unit")
  pss.check(sigma.u, "req"); pss.check(sigma.u, "pos")

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


