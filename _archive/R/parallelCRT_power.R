#' Parallel Cluster Randomized Trial Power Calculation
#'
#' @param alpha The significance level or type 1 error rate
#' @param J Number of clusters
#' @param m The number of individuals in each cluster at each time period
#' @param d The standardized effect size
#' @param rho Intraclass correlation coefficient (ICC)
#' @param rho_c Cluster autocorrelation (proportion of cluster-level variance that is time invariant)
#' @param rho_s Subject autocorrelation (This is 0 in repeated cross sectional designs because people are measured only once)
#'
#' @return power
#' @export
#'
#' @examples
#' # closed cohort design as in Example 12.1
#' parallelCRT_power(alpha = 0.05, J = 16, m = 30, d = 0.3, rho = 0.05, rho_c = 0.4, rho_s = 0.5)
#' # repeated cross sectional design as in Example 12.1
#' parallelCRT_power(alpha = 0.05, J = 16, m = 30, d = 0.3, rho = 0.05, rho_c = 0.4, rho_s = 0)
parallelCRT_power <- function(alpha, J, m, d, rho, rho_c, rho_s){

  # design effect for a parallel CRT:
  DE_pa <- 1 + (m - 1) * rho
  r <- rho_c * (m * rho) / DE_pa + rho_s * (1 - rho) / DE_pa

  lambda <- d / sqrt( (1 - r^2) * DE_pa * 4 / (m * J))
  crit <- stats::qf(p = 1 - alpha, df1 = 1, df2 = J - 2)
  power <- 1 - stats::pf(q = crit, df1 = 1, df2 = J - 2, ncp = lambda^2)

  return(power)

}
