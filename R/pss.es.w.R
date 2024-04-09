#' Cohen's w effect size calculation for chi-square tests
#'
#' @param p0vec The first vector of probabilities. Both p0vec and p1vec, or pmatrix must be specified.
#' @param p1vec The second vector of probabilities. Both p0vec and p1vec, or pmatrix must be specified.
#' @param pmatrix The two-way probability table. Both p0vec and p1vec, or pmatrix must be specified.
#'
#' @return A list of the arguments and the w effect size.
#' @export
#'
#' @examples
#' pss.es.w(p0vec = c(0.5, 0.3, 0.2), p1vec = c(0.7, 0.2, 0.1))
#' pss.es.w(pmatrix = matrix(c(0.050, 0.350, 0.100, 0.075, 0.250, 0.175), nrow = 2, byrow = TRUE))

pss.es.w <- function (p0vec = NULL, p1vec = NULL, pmatrix = NULL) {

  # Check if the arguments are specified correctly
  if(is.null(pmatrix) & (is.null(p0vec) | is.null(p1vec)))
    stop("p0vec and p1vec, or pmatrix must be specified")
  if(!is.null(pmatrix) & (!is.null(p0vec) | !is.null(p1vec)))
    stop("p0vec and p1vec, or pmatrix must be specified")

  # Calculate w
  if (is.null(pmatrix))
    w <- sqrt(sum((p1vec - p0vec)^2 / p0vec))
  else {
    pi <- apply(pmatrix, 1, sum)
    pj <- apply(pmatrix, 2, sum)
    p0 <- pi %*% t(pj)
    w <- sqrt(sum((pmatrix - p0)^2/p0))
  }

  # Print output as a power.htest object
  METHOD <- "Cohen's w effect size calculation for chi-square tests"

  if (is.null(pmatrix))
    structure(list(p0vec = p0vec, p1vec = p1vec, w = w,
                 method = METHOD), class = "power.htest")
  else
    structure(list(pmatrix = pss.matrix.format(pmatrix), w = w,
                   method = METHOD), class = "power.htest")

}
