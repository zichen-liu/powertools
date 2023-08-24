#' Power calculations for two-way balanced analysis of variance contrast test
#'
#' @param a The number of groups for factor A.
#' @param b The number of groups for factor B.
#' @param n1 The sample size for the first cell being compared.
#' @param n2 The sample size for the second cell being compared.
#' @param N The total population.
#' @param mu1 The first cell mean being compared.
#' @param mu2 The second cell mean being compared.
#' @param sd The estimated standard deviation within each group; defaults to 1.
#' @param indx Whether there is an interaction between the two factors.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed power).
#' @export
#'
#' @examples
#' # Example 5.12
#' pss.anova.unbal.2w.se(a = 2, b = 3, n1 = 30, n2 = 30, N = 180, mu1 = 9.3, mu2 = 8.7, sd = 2, intx = TRUE, alpha = 0.025)

pss.anova.unbal.2w.se <- function (a = NULL, b = NULL, n1 = NULL, n2 = NULL,
                                   N = NULL, mu1 = NULL, mu2 = NULL, sd = 1,
                                   intx = FALSE, alpha = 0.05) {

  # Check if the arguments are specified correctly
  if (a < 2 | b < 2)
    stop("number of groups per intervention must be at least 2")
  if (!is.null(n) && n < 2)
    stop("number of observations in each group must be at least 2")
  if(is.null(sd))
    stop("sd must be specified")

  # Get lambda
  df <- ifelse(intx, N - a * b, N - a - b + 1)
  lambda <- (mu2 - mu1) / (sd * sqrt(1 / n1 + 1 / n2))
  power <- stats::pt(q = stats::qt(alpha, df), df, lambda)

  # Generate output text
  ab <- c(a, b)
  METHOD <- "Unalanced two-way analysis of variance\n     simple effects power calculation"
  mrows <- c()
  for (i in 1:a) mrows <- c(mrows, paste(mmatrix[i,], collapse = ', '))
  mmatrix <- paste(mrows, collapse = " | ")

  # Print output as a power.htest object
  structure(list(`a, b` = ab, mu1 = mu1, mu2 = mu2, n1 = n1, n2 = n2,
                 sd = sd, alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")
}

