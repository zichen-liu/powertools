pss.re.crt.bin <- function(m, m.sd, pc, pt, sigma_u){

  logoddsc <- log(pc / (1 - pc))
  logoddst <- log(pt / (1 - pt))
  gam0 <- (logoddsc + logoddst) / 2
  gam1 <- logoddst - logoddsc

  sigmasqet <- 2 + exp(-gam0 - 2 * gam1) + exp(gam0 + 2 * gam1)
  sigmasqec <- 2 + exp(-gam0 + 2 * gam1) + exp(gam0 - 2 * gam1)

  alpha_t <- sigmasqet / sigma_u^2
  alpha_c <- sigmasqec / sigma_u^2

  cv <- m.sd / m

  lamt <- m / (m + alpha_t)
  lamc <- m / (m + alpha_c)

  re.num <- (1 - cv^2 * lamt * (1 - lamt)) * (1 - cv^2 * lamc * (1 - lamc)) * (lamt + lamc)
  re.denom <- lamt + lamc - cv^2 * (lamt^2 * (1 - lamt) + lamc^2 * (1 - lamt))
  re <- re.num / re.denom

  print(re)

}


pss.re.crt.bin(m = 60, m.sd = 45, pc = 0.25, pt = 0.15, sigma_u = 0.3)



