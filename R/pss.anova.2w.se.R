#' Power calculations for two-way balanced analysis of variance simple effects test
#'
#' @param a The number of groups for factor A.
#' @param b The number of groups for factor B.
#' @param n The sample size per group.
#' @param mu1 The first cell mean being compared.
#' @param mu2 The second cell mean being compared.
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
#' pss.anova.2w.se(a = 2, b = 3, n = 30, mu1 = 9.3, mu2 = 8.7, sd = 2, intx = TRUE, alpha = 0.025)

pss.anova.2w.se <- function (a = NULL, b = NULL, n = NULL, mu1 = NULL, mu2 = NULL,
                            sd = 1, intx = FALSE, alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  if (sum(vapply(list(n, alpha, power), is.null, NA)) != 1)
    stop("exactly one of 'n', 'alpha', and 'power' must be NULL")
  if (a < 2 | b < 2)
    stop("number of groups per intervention must be at least 2")
  if (!is.null(n) && n < 2)
    stop("number of observations in each group must be at least 2")
  if(is.null(sd))
    stop("sd must be specified")

  # Get test statistic
  p.body <- quote({
    lambda <- (mu2 - mu1) / (sd * sqrt(1 / n + 1 / n))
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
  mrows <- c()
  for (i in 1:a) mrows <- c(mrows, paste(mmatrix[i,], collapse = ', '))
  mmatrix <- paste(mrows, collapse = " | ")

  # Print output as a power.htest object
  structure(list(`a, b` = ab, mu1 = mu1, mu2 = mu2, n = n,
                 sd = sd, alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")
}

