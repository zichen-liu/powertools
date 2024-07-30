#' Stepped Wedge Design Power Calculation with 1 Treatment
#'
#' Calculates power for stepped wedge design with one treatment and one control condition.
#'
#' @param alpha The significance level or type 1 error rate
#' @param J Number of clusters
#' @param m The number of individuals in each cluster at each time period
#' @param K The number of steps (if one baseline period, then this is periods - 1)
#' @param b The number of baseline periods
#' @param d The standardized effect size
#' @param rho Intraclass correlation coefficient (ICC)
#' @param rho_c Cluster autocorrelation (proportion of cluster-level variance that is time invariant)
#' @param rho_s Subject autocorrelation (This is 0 in repeated cross sectional designs because people are measured only once)
#'
#' @return power
#' @export
#'
#' @examples
#' # closed cohort design as in Example 12.5
#' swd_1trt_power(alpha = 0.05, m = 30, J = 5, K = 5, b = 1, d = 0.4,
#'                rho = 0.05, rho_c = 0.2, rho_s = 0.4)
#' # repeated cross sectional design as in Example 12.5
#' swd_1trt_power(alpha = 0.05, m = 30, J = 5, K = 5, b = 1, d = 0.4,
#'                rho = 0.05, rho_c = 0.2, rho_s = 0)
swd_1trt_power <- function(alpha = 0.05, J, m, K, b, d, rho, rho_c, rho_s){

  DE_crt <- 1 + (m - 1) * rho
  r <- (rho_s + (m * rho_c - rho_s) * rho) / DE_crt
  DE_sw_num <- 3 * K * (1 - r) * (1 + (K + b - 1) * r) * DE_crt
  DE_sw_den <- 2 *(K^2 - 1) * (1 + (0.5 * K + b - 1) * r)
  DE_sw <- DE_sw_num / DE_sw_den

  lambda <- d / sqrt(4 * DE_sw / (m * J))
  df <- J * (K + b) - (K + b) - 1
  df

crit <- stats::qf(p = 1 - alpha, df1 = 1, df2 = df)
power <- 1 - stats::pf(q = crit, df1 = 1, df2 = df, ncp = lambda^2)

return(power)

}
