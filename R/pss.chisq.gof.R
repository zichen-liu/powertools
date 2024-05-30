#' Power calculations for chi-square goodness-of-fit test
#'
#' @param p0vec The first vector of probabilities (under the null).
#' @param p1vec The second vector of probabilities (under the alternative hypothesis).
#' @param N The number of total observations.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.chisq.gof(p0vec = c(0.5, 0.3, 0.2), p1vec = c(0.7, 0.2, 0.1), N = 50)

pss.chisq.gof <- function (p0vec = NULL, p1vec = NULL,
                           N = NULL, alpha = 0.05, power = NULL,
                           v = TRUE) {

  # Check if the arguments are specified correctly
  pss.check.many(list(N, alpha, power), "oneof")
  pss.check(p0vec, "req"); pss.check(p0vec, "vec"); pss.check(p0vec, "sum")
  pss.check(p1vec, "req"); pss.check(p1vec, "vec"); pss.check(p1vec, "sum")
  pss.check(N, "int"); pss.check(N, "min", min = 2)
  pss.check(alpha, "unit")
  pss.check(power, "unit")
  pss.check(v, "req"); pss.check(v, "bool")

  if (length(p0vec) != length(p1vec))
    stop("the two proportion vectors must have the same lengths")

  if (any(p0vec <= 0) | any(p0vec >= 1) | any(p1vec <= 0) | any(p1vec  >= 1))
    stop("all proportions must be between 0 and 1")

  if (sum(p0vec) != 1 | sum(p1vec) != 1)
    stop("proportions must sum to 1 across each vector")

  # Calculate effect size and df
  es <- sqrt(sum((p1vec - p0vec)^2 / p0vec))
  df <- length(p0vec) - 1

  # Calculate test statistic
  p.body <- quote({
    stats::pchisq(stats::qchisq(alpha, df, lower = FALSE),
                  df, N * es^2, lower = FALSE)
  })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(N)) {
    N <- stats::uniroot(function(n) eval(p.body) - power, c(1 + 1e-10, 1e+09))$root
    if (!v) return(N)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- "Chi-square goodness-of-fit power calculation"
  structure(list(p1vec = p1vec, p0vec = p0vec, `w effect size` = es, N = N, alpha = alpha,
                 power = power, method = METHOD), class = "power.htest")
}
