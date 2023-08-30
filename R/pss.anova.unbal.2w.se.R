#' Power calculations for two-way balanced analysis of variance contrast test
#'
#' @param nmatrix A matrix of sample sizes (see example).
#' @param mmatrix A matrix of group means (see example).
#' @param cmatrix A matrix of contrast coefficients (see example).
#' @param sd The estimated standard deviation within each group.
#' @param indx Whether there is an interaction between the two factors.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' # Example 5.12
#' nmatrix <- matrix(c(30, 30, 30, 30, 30, 30), nrow = 2, byrow = TRUE)
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
#' cmatrix <- matrix(c(-1, 0, 0, 1, 0, 0), nrow = 2, byrow = TRUE)
#' pss.anova.unbal.2w.se(nmatrix = nmatrix, mmatrix = mmatrix, cmatrix = cmatrix, sd = 2, intx = TRUE, alpha = 0.025)

pss.anova.unbal.2w.se <- function (nmatrix = NULL, mmatrix = NULL, cmatrix = NULL,
                                   sd = 1, intx = FALSE, alpha = 0.05) {

  # Check if the arguments are specified correctly
  a <- nrow(mmatrix)
  b <- ncol(mmatrix)
  if (a < 2 | b < 2)
    stop("number of groups per intervention must be at least 2")
  if (!is.null(n) && n < 2)
    stop("number of observations in each group must be at least 2")
  if (a != nrow(cmatrix) | b != ncol(cmatrix))
    stop("number of contrast coefficients must be equal to the number of groups")
  if (a != nrow(nmatrix) | b != ncol(nmatrix))
    stop("number of sample sizes must be equal to the number of groups")
  if(is.null(sd))
    stop("sd must be specified")

  # Get lambda
  lambda <- sum(cmatrix * mmatrix) / sd / sqrt(sum(cmatrix^2 / nmatrix))

  # Calculate power
  N <- sum(nmatrix)
  df <- ifelse(intx, N - a * b, N - a - b + 1)
  power <- stats::pt(q = stats::qt(alpha, df), df, lambda)

  # Generate output text
  ab <- c(a, b)
  METHOD <- "Balanced two-way analysis of variance\n     simple effects power calculation"

  # Print output as a power.htest object
  structure(list(`a, b` = ab, mmatrix = pss.matrix.format(mmatrix),
                 cmatrix = pss.matrix.format(cmatrix),
                 nmatrix = pss.matrix.format(nmatrix),
                 sd = sd, alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")
}

