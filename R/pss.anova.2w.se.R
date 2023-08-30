#' Power calculations for two-way balanced analysis of variance simple effects test
#'
#' @param n The sample size per group.
#' @param mmatrix A matrix of group means (see example).
#' @param cmatrix A matrix of contrast coefficients (see example).
#' @param sd The estimated standard deviation within each group; defaults to 1.
#' @param indx Whether there is an interaction between the two factors.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' # Example 5.12
#' mmatrix <- matrix(c(9.3, 8.9, 8.5, 8.7, 8.3, 7.3), nrow = 2, byrow = TRUE)
#' cmatrix <- matrix(c(-1, 0, 0, 1, 0, 0), nrow = 2, byrow = TRUE)
#' pss.anova.2w.se(n = 30, mmatrix = mmatrix, cmatrix = cmatrix, sd = 2, intx = TRUE, alpha = 0.025)

pss.anova.2w.se <- function (n = NULL, mmatrix = NULL, cmatrix = NULL,
                            sd = 1, intx = FALSE, alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  a <- nrow(mmatrix)
  b <- ncol(mmatrix)
  if (sum(vapply(list(n, alpha, power), is.null, NA)) != 1)
    stop("exactly one of 'n', 'alpha', and 'power' must be NULL")
  if (a < 2 | b < 2)
    stop("number of groups per intervention must be at least 2")
  if (!is.null(n) && n < 2)
    stop("number of observations in each group must be at least 2")
  if (a != nrow(cmatrix) | b != ncol(cmatrix))
    stop("number of contrast coefficients must be equal to the number of groups")
  if(is.null(sd))
    stop("sd must be specified")

  # Get test statistic
  p.body <- quote({
    lambda <- sum(cmatrix * mmatrix) / sd / sqrt(sum(cmatrix^2 / n))
    N <- a * b * n
    df <- ifelse(intx, N - a * b, N - a - b + 1)
    stats::pt(q = stats::qt(alpha, df), df, lambda)
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+05))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error", domain = NA)

  # Generate output text
  ab <- c(a, b)
  METHOD <- "Balanced two-way analysis of variance\n     simple effects power calculation"

  # Print output as a power.htest object
  structure(list(`a, b` = ab, mmatrix = pss.matrix.format(mmatrix),
                 cmatrix = pss.matrix.format(cmatrix), n = n,
                 sd = sd, alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")
}

