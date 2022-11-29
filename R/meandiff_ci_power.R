#' Calculate Power of a confidence interval for a Difference of Two Means
#'
#' This function calculates the power of a confidence interval for the difference of two means accounting for estimating sigma and allowing for unequal allocation between the two groups.
#'
#' @param n1 The sample size in group 1
#' @param n2 The sample size in group 2
#' @param alpha The significance level
#' @param d The standardized halfwidth (halfwidth / sigma). Either d OR halfwidth and sigma need specified.
#' @param halfwidth The halfwidth; half of the confidence interval width. If halfwidth is specified, sigma needs specified too.
#' @param sigma The estimated standard deviation. If sigma is specified, halfwidth needs specified too.
#'
#' @return The total N, group sample sizes, and corresponding power.
#' @export
#'
#' @examples
#' # Calculate power for a 95% confidence interval to have 0.25 standard deviation width
#' # when there are 100 and 160 people in each group
#' meandiff_ci_power(n1 = 100, n2 = 160, alpha = 0.05, d = 0.25, halfwidth = NULL, sigma = NULL)
meandiff_ci_power <- function(n1, n2, alpha, d = NULL, halfwidth = NULL, sigma = NULL){

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

  N <- n1 + n2
  r <- n2 / n1
  nu <- n1 * (r + 1) - 2
  numerator <- nu * n1 * r * d^2
  denominator <- stats::qt(1 - alpha/2, N-2)^2 * (r + 1)

	power <- stats::pchisq(numerator / denominator, N-2)
  return(data.frame(N = N, n1 = n1, n2 = n2, CIpower = power))
}
