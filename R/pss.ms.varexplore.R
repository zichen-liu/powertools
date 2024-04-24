pss.ms.varexplore <- function(pc, pt){

  or <- (pt / (1 - pt)) / (pc / (1 - pc))

  logoddsc <- log(pc / (1 - pc))
  logoddst <- log(pt / (1 - pt))
  gam1 <- logoddst - logoddsc

  sigma_u1 <- seq(0.1, 1, by = 0.1)
  or_lo <- exp(gam1 - 1.96 * sigma_u1)
  or_lo <- round(or_lo, 2)
  or_up <- exp(gam1 + 1.96 * sigma_u1)
  or_up <- round(or_up, 2)

  out <- cbind(sigma_u1, or_lo, or_up)

  return(list(or = or, out = out))

}

pss.ms.varexplore(pc = 0.1, pt = 0.2)



## crt sigma u explorer

pss.crt.varexplore <- function(pc, pt){
  logoddsc <- log(pc / (1 - pc))
  logoddst <- log(pt / (1 - pt))
  gam0 <- (logoddsc + logoddst) / 2
  gam1 <- logoddst - logoddsc

  sigma_u <- seq(0.1, 1, by = 0.1)
  pc_lo <- exp(gam0 + gam1 * -0.5 - 1.96 * sigma_u)/(1 + exp(gam0 + gam1 * -0.5 - 1.96 * sigma_u))
  pc_lo <- round(pc_lo, 2)
  pc_up <- exp(gam0 + gam1 * -0.5 + 1.96 * sigma_u)/(1 + exp(gam0 + gam1 * -0.5 + 1.96 * sigma_u))
  pc_up <- round(pc_up, 2)

  pt_lo <- exp(gam0 + gam1 * 0.5 - 1.96 * sigma_u)/(1 + exp(gam0 + gam1 * 0.5 - 1.96 * sigma_u))
  pt_lo <- round(pt_lo, 2)
  pt_up <- exp(gam0 + gam1 * 0.5 + 1.96 * sigma_u)/(1 + exp(gam0 + gam1 * 0.5 + 1.96 * sigma_u))
  pt_up <- round(pt_up, 2)

  pc_pt <- c(pc, pt)
  out <- cbind(sigma_u, pc_lo, pc_up, pt_lo, pt_up)
  return(list(pc_pt = pc_pt, out = out))
}


pss.crt.varexplore(pc = 0.25, pt = 0.15)
