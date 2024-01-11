#' Cohen's d effect size calculation for one or two means
#'
#' @param delta If one mean: muA (the true mean) - mu0 (the mean under the null). If two means: DeltaA (the true difference mu1 - mu2) - Delta0 (the difference under the null).
#' @param sd The estimated standard deviation; defaults to 1.
#'
#' @return A list of the arguments and the d effect size.
#' @export
#'
#' @examples
#' pss.es.d(delta = 6.5 - 5.7, sd = 0.4)

pss.es.d <- function (delta = NULL, sd = 1) {

  # Check if the arguments are specified correctly
  if(is.null(delta) | is.null(sd))
    stop("delta and sd must be specified")

  # Calculate d
  d <- abs(delta) / sd

  # Print output as a power.htest object
  METHOD <- "Cohen's d effect size calculation for one or two means"
  structure(list(delta = delta, sd = sd, d = d,
                 method = METHOD), class = "power.htest")

}
