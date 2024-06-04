#' Sample Size Calculation for One Mean Confidence Interval
#'
#' This function calculates the minimum sample size required to achieve certain precision in estimating a confidence interval for a single mean. This function accounts for estimation in sigma and searches over a range of supplied N.
#' 'Power' used here is the power of a confidence interval, which is not true statistical power; instead it is the probability of the interval obtaining a certain halfwidth.
#'
#' @param N_min The minimum number of observations to explore power.
#' @param N_max The maximum number of observations to explore power.
#' @param d The standardized halfwidth (halfwidth / standard deviation).
#' @param alpha The significance level.
#' @param power The desired level of power.
#'
#' @return The minimum sample size required to achieve the specified level of power.
#' @export
#'
#' @keywords internal
#'
#' @examples
#' # Example 1: Find the required sample size that achieves a half width of 0.25 standard deviations
#' # at significance level of 0.05 with 80% power. We specify to search between 62 and 75 people.
#' onemean_ci_ss(N_min = 62, N_max = 75, d=0.25, alpha=0.05, power=0.8)
onemean_ci_ss <- function(N_min, N_max, d, alpha, power){
    N_search <- seq(N_min, N_max)
    power_search <- search_power(onemean_ci_power, N_search, d=d, alpha=alpha)

    goodpower <- power_search[power_search$CIpower >= power, ]

    if(nrow(goodpower) > 0){
      return(print(as.data.frame(goodpower[1, ]), row.names=F))
    }else{
      print(paste0("Need larger N for this level of power"))
      last_obs <- nrow(power_search)
      return(print(as.data.frame(power_search[last_obs, ]), row.names=F))
    }
}
