#' Power calculations for chi-square goodness-of-fit test
#'
#' @param p0vec The first vector of probabilities (under the null).
#' @param p1vec The second vector of probabilities (under the alternative hypothesis).
#' @param N The number of total observations.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.chisq.gof(p0vec = c(0.5, 0.3, 0.2), p1vec = c(0.7, 0.2, 0.1), N = 50)

pss.chisq.gof <- function (p0vec = NULL, p1vec = NULL,
                           N = NULL, alpha = 0.05, power = NULL) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(N, alpha, power), is.null)) != 1)
    stop("exactly one of N, alpha, and power must be NULL")
  if (is.null(p0vec) | is.null(p1vec))
    stop("please specify the two proportion vectors p0 and p1")
  if (!is.null(N) && any(N < 1))
    stop("number of observations must be at least 1")
  if (length(p0vec) != length(p1vec))
    stop("the two proportion vectors must have the same lengths")

  # Calculate effect size and df
  es <- sqrt(sum((p1vec - p0vec)^2 / p0vec))
  df <- length(p0vec) - 1

  # Calculate test statistic
  p.body <- quote({
    stats::pchisq(stats::qchisq(alpha, df, lower = FALSE),
                  df, N * es^2, lower = FALSE)
  })

  # Use uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(N))
    N <- uniroot(function(n) eval(p.body) - power, c(1 + 1e-10, 1e+09))$root
  else if (is.null(alpha))
    alpha <- uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "Chi-square goodness-of-fit power calculation"
  structure(list(p1vec = p1vec, p0vec = p0vec, effect.size = es, N = N, alpha = alpha,
                 power = power, method = METHOD), class = "power.htest")
}
