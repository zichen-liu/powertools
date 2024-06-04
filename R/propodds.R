#' Power calculations for ordinal categorical variable under proportional odds assumption
#'
#' @param pC Vector of response probabilities in control group (group 1). Must sum to 1. Categories are ordered from best to worst.
#' @param OR Odds ratio when the alternative is true. Must be greater than 1.
#' @param n1 Sample size for group 1 (control group).
#' @param n.ratio The ratio n2/n1 between the sample sizes of two groups; defaults to 1 (equal group sizes).
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power; defaults to 0.8.
#' @import Hmisc
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' library(Hmisc)
#' pC <- c(0.2, 0.5, 0.2, 0.1)
#' propodds(pC = pC, OR = 2.5, n1 = 65, n.ratio = 1, alpha = 0.05)

propodds <- function(pC, OR, n1, n.ratio = 1, alpha = 0.05, power = NULL){

  if (sum(sapply(list(n1, power), is.null)) != 1)
    stop("exactly one of n1 and power must be NULL")
  if (OR<=1)
    stop("OR must be greater than 1")

  pC <- pC[!is.na(pC)]
  if (abs(sum(pC) - 1) > 1e-05)
    stop("probabilities in pC do not sum to 1")
  pT <- Hmisc::pomodm(p = pC, odds.ratio = 1 / OR)
  pavg <- (pC + pT) / 2
  ps <- 1 - sum(pavg^3)
  za <- stats::qnorm(1 - alpha/2)

  if (is.null(power)){
    N <- n1 * (1 + n.ratio)
    n2 <- n1 * n.ratio
    V <- n1 * n2 * N/3/((N + 1)^2) * ps
    power <- stats::pnorm(abs(log(OR)) * sqrt(V) - za)
  }
  else if (is.null(n1)){
    zb <- stats::qnorm(power)
    N <- (3 * (n.ratio + 1)^2 * (zb + za)^2) / (n.ratio * (log(OR))^2 * ps)
    n1 <- N / (1 + n.ratio)
    n2 <- n.ratio * n1
  }
  else stop("internal error")

  METHOD = "Power calculation for ordinal categorical outcome under proportional odds"
  n <- c(n1, n2)
  structure(list(n = n, pC = pC, pT = round(pT, digits = 3),
                 OR = OR, alpha = alpha, power = power, method = METHOD),
            class = "power.htest")
}


