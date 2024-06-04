#' Calculate Sample Size needed for adequate power for Confidence Interval of Difference in Means.
#'
#' Calculate sample sizes in two groups (allowing unequal allocation) that will adequately power the confidence interval of the difference of group means. This function accounts for estimating sigma, the standard deviation.
#'
#' @param N_min The minimum total sample size.
#' @param N_max The maximum total sample size
#' @param ratios A 2 element vector of group allocation ratios. Equal allocation is default, which is specified as c(1,1). A allocation ratio of 2 (r = n2/n1) would be specified as c(1,2).
#' @param alpha The significance level
#' @param power The specified power to achieve
#' @param d The standardized halfwidth (halfwidth / sigma). Either d OR halfwidth and sigma need specified.
#' @param halfwidth The halfwidth; half of the confidence interval width. If halfwidth is specified, sigma needs specified too.
#' @param sigma The estimated standard deviation. If sigma is specified, halfwidth needs specified too.
#'
#' @return returns all sample sizes and power that satisfy the specified power level.
#' @export
#'
#' @keywords internal
#'
#' @examples
#' # Find sample sizes with allocation ratio of 2 between two groups for a 95% confidence
#' # interval to have 80% probability of a 0.25 standard deviation halfwidth or smaller.
#' meandiff_ci_ss(N_min = 200, N_max = 325, r = c(1,2),
#'                alpha = 0.05, power = 0.8, d = 0.25)
meandiff_ci_ss <- function(N_min, N_max, ratios = c(1,1), alpha, power,
                            d = NULL, halfwidth = NULL, sigma = NULL){

  if(is.null(d) & (is.null(halfwidth) | is.null(sigma))){
    print(paste0("Either d (halfwidth/sigma) OR halfwidth and sigma need to be specified"))
    stop()
  }
  if(!is.null(d) & !is.null(halfwidth) & !is.null(sigma)){
    print(paste0("Either d (halfwidth/sigma) OR halfwidth and sigma need to be specified"))
    stop()
  }
  if(!is.null(halfwidth) & !is.null(sigma)){
    d <- halfwidth / sigma
  }

  i_vec <- seq(5, N_max, 1)
  i_vec_mat <- matrix(i_vec, nrow=length(i_vec), ncol = 1)
  ratio_mat <- matrix(ratios, nrow=1, ncol=2)

  N_vecs <- i_vec_mat %*% ratio_mat
  Ntot <- rowSums(N_vecs)
  N_vec_4search <- N_vecs[Ntot >= N_min & Ntot <= N_max, ]


  n1_search <- N_vec_4search[, 1]
  n2_search <- N_vec_4search[, 2]

    power_search <- search_power(meandiff_ci_power, n1_search, n2_search,
                                 d=d, alpha=alpha)

    goodpower <- power_search[power_search$CIpower >= power, ]

    if(nrow(goodpower) > 0){
      return(print(as.data.frame(goodpower), row.names=F))
    }else{
      print(paste0("Need larger N for this level of power"))
      last_obs <- nrow(power_search)
      return(print(as.data.frame(power_search[last_obs, ]), row.names=F))
    }
}
