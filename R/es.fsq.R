#' Cohen's f^2 effect size calculation for overall F test
#'
#' @param Rsq The squared sample multiple correlation coefficient.
#'
#' @return A list of the arguments and the f^2 effect size.
#' @export
#'
#' @examples
#' es.fsq(Rsq = 0.02)

es.fsq <- function (Rsq = 0.02) {

  # Check if the arguments are specified correctly
  if(is.null(Rsq))
    stop("Rsq must be specified")

  # Calculate f^2
  fsq <- Rsq / (1 - Rsq)

  # Print output as a power.htest object
  METHOD <- "Cohen's f^2 effect size calculation for overall F test"
  structure(list(Rsq = Rsq, fsq = fsq,
                 method = METHOD), class = "power.htest")

}
