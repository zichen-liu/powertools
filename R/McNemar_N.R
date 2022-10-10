#### Calculate N pairs from marginal probabilities and rho ####
# Supply EITHER (p1, p2, and rho) OR (p01 and p10)

# this is an update to Kate's function #
#' Sample Size Calculation for McNemar's Test for paired proportions
#'
#' Calculates number of pairs needed to achieve specified level of power from marginal probabilities
#' and the correlation between probabilities (p1, p2, rho) OR with discordant cell probabilities (p01 and p10).
#'
#' @param p1 Marginal probability of success for outcome 1. Either NULL or supplied alongside p2 and rho.
#' @param p2 Marginal probability of success for outcome 2. Either NULL or supplied alongside p1 and rho.
#' @param rho Correlation between proportions. Either NULL or supplied alongside p1 and p2.
#' @param p01 First discordant cell probability. Either NULL or supplied alongside p10.
#' @param p10 Second discordant cell probability. Either NULL or supplied alongside p01.
#' @param alpha Significance level. Default is 0.05
#' @param power Required power. Default is 0.8.
#' @param one_or_twosides Whether the hypothesis test of interest is one or two-sided. Default is "two".
#'
#' @return Returns a data frame with total number of observations needed (N_obs) and number of pairs needed (N_pairs).
#' @export
#'
#' @examples
#' # supplying marginal probabilities
#' McNemar_N(p1 = 0.8, p2 = 0.9, rho = 0, alpha = 0.05, power = 0.9)
#' # supplying cell probabilities
#' McNemar_N(p10 = 0.18, p01 = 0.08, alpha = 0.05, power = 0.9)


McNemar_N <- function(p1 = NULL, p2 = NULL, rho = NULL, p01 = NULL, p10 = NULL,
                      alpha = 0.05, power = 0.8, one_or_twosides = "two"){

  # if they supply marginal and cell proportions end function
  if((!is.null(p1) | !is.null(p2)) & (!is.null(p01) | !is.null(p10))){
    stop("Please supply either p1 and p2 OR p01 and p10")
  }

  siglevel <- ifelse(one_or_twosides == "two", alpha/2, alpha)

  # calculate p01 and p10 if not supplied
  if(!is.null(p1) & !is.null(p2) & !is.null(rho)){
    p01 <- p2*(1-p1) - rho * sqrt((1-p1)*p1*(1-p2)*p2)
    p10 <- p01 + p1 - p2
  }

  # stop function if no p01 and p10
  if(is.null(p01) | is.null(p10)){
    stop("Either there was an issue calculating p01 and p10 OR they were not supplied")
  }

  # Normal Approximation
  Normal_Npairs <- (stats::qnorm(1-siglevel)*sqrt(p01+p10) + stats::qnorm(power)*sqrt(p10+p01-(p10-p01)^2))^2/(p10-p01)^2
  Normal_Npairs <- ceiling(Normal_Npairs)
  Normal_N <- Normal_Npairs * 2

  df_out <- data.frame("Normal_Approx" = c(Normal_N, Normal_Npairs))
  rownames(df_out) <- c("N_obs", "N_pairs")

  return(df_out)
}
