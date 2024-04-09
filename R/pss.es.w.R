#' Cohen's w effect size calculation for chi-square test for goodness-of-fit
#'
#' @param p0vec The first vector of probabilities.
#' @param p1vec The second vector of probabilities.
#'
#' @return A list of the arguments and the w effect size.
#' @export
#'
#' @examples
#' pss.es.w(p0vec = c(0.5, 0.3, 0.2), p1vec = c(0.7, 0.2, 0.1))

pss.es.w <- function (p0vec = NULL, p1vec = NULL) {

  # Check if the arguments are specified correctly
  if(is.null(p0vec) | is.null(p1vec))
    stop("p0vec and p1vec must be specified")

  # Calculate w
  w <- sqrt(sum((p1vec - p0vec)^2 / p0vec))

  # Print output as a power.htest object
  METHOD <- "Cohen's w effect size calculation for chi-square goodness-of-fit"
  structure(list(p0vec = p0vec, p1vec = p1vec, w = w,
                 method = METHOD), class = "power.htest")

}
