#' Calculate Power for a Confidence Interval for One Mean
#'
#' Calculates the power of a confidence interval for a single mean accounting for estimation in sigma, the standard deviation.
#'
#' @param N The number of observations
#' @param d The standardized halfwidth (halfwidth / standard deviation)
#' @param alpha The level of significance
#'
#' @return The statistical power associated with the number of observations.
#' @export
#'
#' @keywords internal
#'
#' @examples
#' # Calculate probability that a 95% confidence interval for group size of 73
#' # to have 0.25 standard deviation halfwidth or less
#' onemean_ci_power(N = 73, d=0.25, alpha=0.05)
onemean_ci_power <- function(N, d, alpha){
	  power <- stats::pchisq( N * (N-1) * d^2 / (stats::qt(1-alpha/2, N-1))^2, df=N-1 )
	  return(data.frame(N = N, CIpower = power))
}
