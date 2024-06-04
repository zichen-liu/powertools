#' Algorithm and sample size for multiple co-primary endpoint clinical trials
#'
#' Calculate the solution to Ck, the integral equation to calculate power for multiple endpoints introduced by Sozu et al in Sample Size Determination in Clinical Trials with Multiple Endpoints 2015.
#'
#' The sample size provided is the solution to equation (4.5) in the Sozu text. This is the simple formula provided for practical use to calculate sample size requirements as if there were a single endpoint.
#' Specifically, this formula is n = (Ck + z_alpha)^2 / (kappa * delta_K^2)
#'
#' @param K The number of K co-primary continuous endpoints. K >= 2.
#' @param alpha The significance level.
#' @param power The power level.
#' @param rho The correlation matrix (K x K) that describes the relationship between each endpoint. Diagonal entries are 1.
#' @param gamma The ratio of effect size(s). This can be either a value (K = 2) or a vector (K > 2).
#' @param a A vector of length K. For continuous endpoints, a_1 = ... = a_K = 1. For binary endpoints using the chi-square method (without CC), this is a vector where each element is the ratio between sigma_k^star and sigma_k.
#' @param r The allocation ratio between the number of people in the control and treatment condition (n_C / n_T).
#' @param delta The standardized effect size (value or vector) for each endpoint.
#'
#' @return A dataframe with the solution to the Ck algorithm and the sample size requirements per group.
#' @export
#'
#' @examples
#' # Following Sozu et al Example in Appendix D.1 with 2 continuous endpoints with equal allocation
#' multendpoints_Ck_ss(K = 2, alpha = 0.025, power = 0.8,
#'                     rho = matrix(c(1, 0.5,
#'                                    0.5, 1), ncol = 2),
#'                     gamma = 8/7, a = c(1, 1), r = 1, delta = c(0.4, 0.35))
#' # Following Sozu et al Example in Appendix D.1 with 3 continuous endpoints and equal allocation
#' delta_vector <- c(0.5, 0.45, 0.4)
#' gamma_vector <- c(delta_vector[1]/delta_vector[3], delta_vector[2]/delta_vector[3])
#' rho_matrix <- matrix(c(1, 0.8, 0.8,
#'                        0.8, 1, 0.5,
#'                        0.8, 0.5, 1), ncol = 3)
#'
#' multendpoints_Ck_ss(K = 3, alpha = 0.025, power = 0.8, rho = rho_matrix,
#'                     gamma = gamma_vector, a = c(1, 1, 1), r = 1,
#'                     delta = delta_vector)
multendpoints_Ck_ss <- function(K, alpha, power, rho, gamma, a, r, delta){

  if(det(rho) <= 0){
    print("no positive definite")
    stop()
  }

  z_a = stats::qnorm(1 - alpha)
  ndel = 0.001
  rn = round(stats::runif(1) * 1000)

  # Initial estimation of CK
  CK = mvtnorm::qmvnorm(power, corr = rho, tail = "lower.tail")$quantile

  # Begin: Newton-Raphson algorithm to find CK
  for(j in 1:1000){
    set.seed(rn)
    C1k = CK * gamma + z_a * (a[K] * gamma - a[1:(K - 1)])
    pow1 = mvtnorm::pmvnorm(lower = rep(-Inf, K), upper = c(C1k, CK), corr = rho)[1]
    G = power - pow1

    if(abs(G) < 0.00001 & G <= 0){
      break
    }
    F1k = rep(0, K - 1)

    for(l in 1:(K - 1)){
      vndel = rep(0, K - 1)
      vndel[l] = ndel
      F1k[l] = mvtnorm::pmvnorm(lower = c(rep(-Inf, l - 1), C1k[l], rep(-Inf, K - l)),
      upper = c(C1k + vndel, CK), corr = rho)[1] / ndel
    }

    FK = mvtnorm::pmvnorm(lower = c(rep(-Inf, K - 1), CK), upper = c(C1k, CK + ndel),
                 corr = rho)[1] / ndel
    dG = -t(F1k) %*% gamma - FK
    CK = CK - G / dG
  }
  # End: Newton-Raphson algorithm to find CK


  # add sample size calculation here now
  kappa <- r / (1 + r)
  n <- (CK + stats::qnorm(alpha, lower.tail = FALSE))^2 / (kappa * (delta[K])^2)
  n <- ceiling(n)

  df_return <- data.frame(Ck = CK, SampleSize = n)
  return(df_return)
}
