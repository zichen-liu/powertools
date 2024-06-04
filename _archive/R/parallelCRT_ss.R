#' Parallel Cluster Randomized Trial Sample Size Calculation
#'
#' @param alpha The significance level or type 1 error rate
#' @param power The specified level of power to achieve
#' @param m The number of individuals in each cluster at each time period
#' @param d The standardized effect size
#' @param rho Intraclass correlation coefficient (ICC)
#' @param rho_c Cluster autocorrelation (proportion of cluster-level variance that is time invariant)
#' @param rho_s Subject autocorrelation (This is 0 in repeated cross sectional designs because people are measured only once)
#'
#' @return Calculated number of clusters required, adjusted number of clusters, and suggested number of clusters that includes the value of K
#' @export
#'
#' @examples
#' # example 12.2
#' parallelCRT_ss(alpha = 0.05, power = 0.8, m = 30, d = 0.3, rho = 0.05, rho_c = 0.4, rho_s = 0.5)
parallelCRT_ss <- function(alpha, power, m, d, rho, rho_c, rho_s){

  # design effect for a parallel CRT:
  DE_pa <- 1 + (m - 1) * rho
  r <- rho_c * (m * rho) / DE_pa + rho_s * (1 - rho) / DE_pa

  za <- stats::qnorm(p = 1 - alpha/2)
  zb <- stats::qnorm(p = power)

  J <- (za + zb)^2 * 4 * DE_pa * (1 - r^2) / (m * d^2)
  df <- J - 2

  adj_ratio <- tandz_ratio(alpha = alpha, power = power, df = df)
  J_adjusted <- J * adj_ratio

  clusters_needed <- plyr::round_any(J_adjusted, 2, f = ceiling)

  paste0("Exact J: ", round(J, 3), "; Adjusted J: ", round(J_adjusted, 3),
          "; Suggested J = ", clusters_needed)

}
