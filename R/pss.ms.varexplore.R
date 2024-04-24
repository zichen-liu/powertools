library(knitr)

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

  out <- data.frame("sigma.u1" = sigma_u1, "OR.lower" = or_lo, "OR.upper" = or_up)
  table <- kable(out, caption = paste("OR:", or), "simple")
  table <- gsub("^Table:", "", table)
  return(table)
}

pss.ms.varexplore(pc = 0.1, pt = 0.2)
