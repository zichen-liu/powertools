#' Cross Over Power Calculation
#'
#' Calculate power for either repeated cross sectional or closed cohort cross over designs.
#'
#' @param design Either "RCS" for repeated cross-sectional design or "Cohort" for closed cohort design. Default is "RCS".
#' @param alpha The significance level or type 1 error rate
#' @param J Number of clusters
#' @param m The number of individuals in each cluster at each time period
#' @param rho Intraclass correlation coefficient (ICC)
#' @param rho_b Between period ICC (correlation between outcomes of two individuals in the same cluster but different time periods)
#' @param xi Within-cluster, within-subject correlation (correlation between two measurements within the same subject). Note this is different from subject autocorrelation. We expect xi to be larger than rho and rho_b. This is 0 for RCS cross over designs.
#' @param d The standardized effect size
#'
#' @return power
#' @export
#'
#' @examples
#' # Repeated cross sectional cross over power calculation as in Example 12.3
#' crossover_power(design = "RCS", alpha = 0.05, J = 8, m = 30, rho = 0.05,
#'                 rho_b = 0.025, xi = 0, d = 0.3)
#' # Closed cohort cross over power calculation as i Example 12.3
#' crossover_power(design = "Cohort", alpha = 0.05, J = 8, m = 30, rho = 0.05,
#'                 rho_b = 0.025, xi = .4, d = 0.3)
crossover_power <- function(design = "RCS", alpha, J, m, rho, rho_b, xi = 0, d){

  if(design == "RCS"){
    if(xi != 0){
      print("Assuming Repeated Cross Sectional Cross Over design that does NOT have repeated measures on same person, ie setting xi = 0")
      xi <- 0
    }
      lambda <- d / sqrt(2 * (1 + (m - 1) * rho - m * rho_b) / (m * J))
  }

  if(design == "Cohort"){
    if(xi == 0){
      print("xi = 0 implies RCS design")
    }
    lambda <- d / sqrt(2 * (1 - xi + (m - 1) * (rho - rho_b)) / (m * J))
  }

  if(design != "RCS" & design != "Cohort"){
    print("Please specify either RCS or Cohort design")
  }

  crit <- stats::qf(p = 1 - alpha, df1 = 1, df2 = 2 * J - 3)
  power <- 1 - stats::pf(q = crit, df1 = 1, df2 = 2 * J - 3, ncp = lambda^2)

  return(power)

}
