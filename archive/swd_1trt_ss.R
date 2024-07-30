#' Stepped Wedge Design Sample Size Calculation with 1 Treatment
#'
#' Calculates sample size requirements for stepped wedge design with one treatment and one control condition.
#'
#' @param alpha The significance level or type 1 error rate
#' @param power The specified level of power
#' @param m The number of individuals in each cluster at each time period
#' @param K The number of steps (if one baseline period, then this is periods - 1)
#' @param b The number of baseline periods
#' @param d The standardized effect size
#' @param rho Intraclass correlation coefficient (ICC)
#' @param rho_c Cluster autocorrelation (proportion of cluster-level variance that is time invariant)
#' @param rho_s Subject autocorrelation (This is 0 in repeated cross sectional designs because people are measured only once)
#'
#' @return Calculated number of clusters required, adjusted number of clusters, and suggested number of clusters that includes the value of K
#' @export
#'
#' @examples
#' # Repeated Cross Sectional Stepped Wedge Sample Size Calculation
#' # with 3 steps, 1 baseline period, 80% power, 0.05 significance, 30 ppl per
#' # cluster, effect size of 0.4, ICC of 0.05, and cluster autocorrelation of 0.2
#' # Following Example 12.6 in text
#' swd_1trt_ss(alpha = 0.05, power = 0.8, m = 30, K = 3, b = 1, d = 0.4,
#'             rho = 0.05, rho_c = 0.2, rho_s = 0)
swd_1trt_ss <- function(alpha, power, m, K, b, d, rho,
                                rho_c, rho_s){
  za <- stats::qnorm(p = 1 - alpha/2)
  zb <- stats::qnorm(p = power)

  DE_crt <- 1 + (m - 1) * rho
  r <- (rho_s + (m * rho_c - rho_s) * rho) / DE_crt
  DE_sw_num <- 3 * K * (1 - r) * (1 + (K + b - 1) * r) * DE_crt
  DE_sw_den <- 2 *(K^2 - 1) * (1 + (0.5 * K + b - 1) * r)
  DE_sw <- DE_sw_num / DE_sw_den

  J <- 4 * (za + zb)^2 * DE_sw / (m * d^2)
  df <- J * (K + b) - (K + b) - 1

  adj_ratio <- tandz_ratio(alpha = alpha, power = power, df = df)
  J_adjusted <- J * adj_ratio

  clusters_needed <- plyr::round_any(J_adjusted, K, f = ceiling)

  paste0("Exact J: ", round(J, 3), "; Adjusted J: ", round(J_adjusted, 3),
          "; Suggested J = ", clusters_needed)

}
