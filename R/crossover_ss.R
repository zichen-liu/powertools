#' Cross Over Sample Size Calculation
#'
#' Calculate sample size requirements for either repeated cross sectional or closed cohort cross over designs.
#'
#' @param design Either "RCS" for repeated cross-sectional design or "Cohort" for closed cohort design. Default is "RCS".
#' @param alpha The significance level or type 1 error rate
#' @param power Specified power level to achieve
#' @param m The number of individuals in each cluster at each time period
#' @param rho Intraclass correlation coefficient (ICC)
#' @param rho_b Between period ICC (correlation between outcomes of two individuals in the same cluster but different time periods)
#' @param xi Within-cluster, within-subject correlation (correlation between two measurements within the same subject). Note this is different from subject autocorrelation. We expect xi to be larger than rho and rho_b. This is 0 for RCS cross over designs.
#' @param d The standardized effect size
#'
#' @return Calculated number of clusters required and adjusted number of clusters to account for normal approximation
#' @export
#'
#' @examples
#' # repeated cross sectional sample size as in Example 12.4
#' crossover_ss(design = "RCS", alpha = 0.05, power = 0.8, m = 30, rho = 0.05,
#'              rho_b = 0.025, xi = 0, d = 0.3)
#' # closed cohort sample size requirement
#' crossover_ss(design = "Cohort", alpha = 0.05, power = 0.8, m = 30,
#'              rho = 0.05, rho_b = 0.025, xi = 0.4, d = 0.3)
crossover_ss <- function(design = "RCS", alpha, power, m, rho, rho_b, xi = 0, d){

  za <- stats::qnorm(p = 1 - alpha/2)
  zb <- stats::qnorm(p = power)

  if(design == "RCS"){
    if(xi != 0){
      print("Assuming Repeated Cross Sectional Cross Over design that does NOT have repeated measures on same person, ie setting xi = 0")
      xi <- 0
    }
    DE <- 1 + (m - 1) * rho - m * rho_b
    J <- (za + zb)^2 * 2 * DE / (m * d^2)
    J_adjusted <- plyr::round_any(J + 2, 2, f = ceiling)

    # lambda <- d / sqrt(2 * (1 + (m - 1) * rho - m * rho_b) / (m * J_adjusted))

  }

  if(design == "Cohort"){
    if(xi == 0){
      print("xi = 0 implies RCS design")
    }
    DE <- 1 - xi + (m - 1) * (rho - rho_b)
    J <- (za + zb)^2 * 2 * DE / (m * d^2)
    J_adjusted <- plyr::round_any(J + 2, 2, f = ceiling)

    # lambda <- d / sqrt(2 * (1 - xi + (m - 1) * (rho - rho_b)) / (m * J))

  }

  if(design != "RCS" & design != "Cohort"){
    print("Please specify either RCS or Cohort design")
  }

  paste0("Exact J: ", round(J, 3), "; Suggested J: ", round(J_adjusted, 3))

}
