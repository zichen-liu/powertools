#' Variance exploration for multisite trials with binary outcomes
#'
#' @param pc The probability of the outcome in control clusters.
#' @param pt The probability of the outcome in treatment clusters.
#'
#' @return A list of the arguments and a dataframe of outputs.
#' @import knitr
#' @export
#'
#' @examples pss.ms.varexplore(pc = 0.1, pt = 0.2)

pss.ms.varexplore <- function(pc, pt) {

  or <- (pt / (1 - pt)) / (pc / (1 - pc))

  logoddsc <- log(pc / (1 - pc))
  logoddst <- log(pt / (1 - pt))
  gam1 <- logoddst - logoddsc

  sigma_u1 <- seq(0.1, 1, by = 0.1)
  or_lo <- exp(gam1 - 1.96 * sigma_u1)
  or_lo <- round(or_lo, 2)
  or_up <- exp(gam1 + 1.96 * sigma_u1)
  or_up <- round(or_up, 2)

  table <- data.frame("sigma.u1" = sigma_u1, "OR.lower" = or_lo, "OR.upper" = or_up)
  out <- kable(table, caption = paste("OR:", or), "simple")
  print(gsub("^Table:", "", out))
  return(invisible(list(OR = or, table = table)))

}
