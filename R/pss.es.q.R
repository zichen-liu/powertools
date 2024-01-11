#' Cohen's q effect size calculation for two correlation coefficients
#'
#' @param rho1 The correlation coefficient in group 1.
#' @param rho2 The correlation coefficient in group 2.
#'
#' @return A list of the arguments and the q effect size.
#' @export
#'
#' @examples
#' pss.es.q(rho1 = 0.3, rho2 = 0.1)

pss.es.q <- function (rho1 = NULL, rho2 = NULL) {

  # Check if the arguments are specified correctly
  if(is.null(rho1) | is.null(rho2))
    stop("rho1 and rho2 must be specified")

  # Calculate q
  rhoprime1 <- 0.5 * log((1 + rho1)/(1 - rho1))
  rhoprime2 <- 0.5 * log((1 + rho2)/(1 - rho2))
  q <- abs(rhoprime1 - rhoprime2)

  # Print output as a power.htest object
  METHOD <- "Cohen's q effect size calculation for two correlation coefficients"
  structure(list(rho1 = rho1, rho2 = rho2, q = q,
                 method = METHOD), class = "power.htest")

}
