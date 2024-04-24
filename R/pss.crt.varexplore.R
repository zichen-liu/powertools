library(knitr)

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

  out <- data.frame("sigma.u" = sigma_u, "pc.lower" = pc_lo, "pc.upper" = pc_up,
                    "pt.lower" = pt_lo, "pt.upper" = pt_up)
  table <- kable(out, caption = paste("pc:", pc, "; pt:", pt), "simple")
  table <- gsub("^Table:", "", table)
  return(table)
}


pss.crt.varexplore(pc = 0.25, pt = 0.15)
