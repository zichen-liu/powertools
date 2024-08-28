#' Power calculation for a multiple linear regression overall F test
#'
#' @description
#' Conducts power and sample size calculations for an overall (or omnibus) F test
#' in a multiple linear regression model.
#' This is a test that all coefficients other than the intercept are equal to zero.
#' Can solve for power, N or alpha.
#'
#' @details
#' Additional details...
#'
#'
#'
#' @param N The sample size.
#' @param p The number of predictors.
#' @param Rsq The squared population multiple correlation coefficient.
#' @param fsq The f-squared effect size. Either Rsq OR fsq must be specified.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param random Whether the values of the predictors are random; defaults to FALSE.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' mlrF.overall(N = 400, p = 2, Rsq = 0.02)
#' mlrF.overall(N = 400, p = 2, fsq = 0.02 / (1 - 0.02))
#' mlrF.overall(N = 109, p = 1, Rsq = 0.3^2)
#' mlrF.overall(N = 50, p = 1, Rsq = 0.2)
#' mlrF.overall(N = 50, p = 3, Rsq = 0.2)
#' mlrF.overall(N = 50, p = 5, Rsq = 0.2)
#' mlrF.overall(N = 400, p = 2, Rsq = 0.02, random = TRUE)

mlrF.overall <- function (N = NULL, p = NULL, Rsq = NULL, fsq = NULL,
                          alpha = 0.05, power = NULL, random = FALSE,
                          v = FALSE) {

  # Check if the arguments are specified correctly
  check.many(list(N, alpha, power), "oneof")
  check.many(list(Rsq, fsq), "oneof")
  check(N, "pos"); check(N, "min", min = 4)
  check(alpha, "unit")
  check(power, "unit")
  check(p, "req"); check(p, "int")
  check(Rsq, "unit")
  check(fsq, "unit")
  check(random, "req"); check(random, "bool")
  check(v, "req"); check(v, "bool")

  if (!is.null(fsq))
    Rsq <- fsq / (1 + fsq)
  else if (!is.null(Rsq))
    fsq <- Rsq / (1 - Rsq)

  # Calculate power
  p.body <- quote({
    if (random) {
      V <- N - p - 1
      U <- (V * Rsq + p)^2 / (V * Rsq * (2 - Rsq) + p)
      crit <- stats::qf(1 - alpha, p, V)
      1 - stats::pf(crit * p * (1 - Rsq) / (V * Rsq + p), U, V)
    } else {
      ncp <- N * Rsq / (1 - Rsq)
      df2 <- N - p - 1
      1 - stats::pf(stats::qf(1 - alpha, p, df2), p, df2, ncp)
    }
  })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power)) {
    power <- eval(p.body)
    if (!v) return(power)
  }
  else if (is.null(N)) {
    N <- stats::uniroot(function(N) eval(p.body) - power, c(4, 1e+09))$root
    if (!v) return(N)
  }
  else if (is.null(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
    if (!v) return(alpha)
  }
  else stop("internal error")

  # Generate output text
  METHOD <- paste0("Power calculation for a multiple linear regression\n     overall F test",
                   " assuming ", ifelse(random, "random", "fixed"), " predictors")

  # Print output as a power.htest object
  structure(list(N = N, p = p, Rsq = Rsq, fsq = fsq,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
