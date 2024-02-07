#' Power calculations for a multiple linear regression overall F test
#'
#' @param N The sample size.
#' @param p The number of predictors.
#' @param Rsq The squared sample multiple correlation coefficient.
#' @param fsq The squared f effect size. Either Rsq OR fsq must be specified.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param random Whether the values of the predictors are random; defaults to FALSE.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.mlrF.overall(N = 400, p = 2, Rsq = 0.02)
#' pss.mlrF.overall(N = 400, p = 2, fsq = 0.02 / (1 - 0.02))
#' pss.mlrF.overall(N = 109, p = 1, Rsq = 0.3^2)
#' pss.mlrF.overall(N = 50, p = 1, Rsq = 0.2)
#' pss.mlrF.overall(N = 50, p = 3, Rsq = 0.2)
#' pss.mlrF.overall(N = 50, p = 5, Rsq = 0.2)
#' pss.mlrF.overall(N = 400, p = 2, Rsq = 0.02, random = TRUE)

pss.mlrF.overall <- function (N = NULL, p = NULL, Rsq = NULL, fsq = NULL,
                              alpha = 0.05, power = NULL, random = FALSE) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(N, power, alpha), is.null)) != 1)
    stop("exactly one of N, alpha, and power must be NULL")
  if (is.null(p))
    stop("please specify the number of predictors p")
  if (is.null(Rsq) & is.null(fsq))
    stop("please specify Rsq or fsq")
  if (!is.null(N) && any(N < 4))
    stop("number of observations must be at least 4")
  if (!is.null(fsq))
    Rsq <- fsq / (1 + fsq)
  if (!is.null(Rsq))
    fsq <- Rsq / (1 - Rsq)

  # Calculate power
  p.body <- quote({
    if (random) {
      v <- N - p - 1
      u <- (v * Rsq + p)^2 / (v * Rsq * (2 - Rsq) + p)
      crit <- stats::qf(1 - alpha, p, v)
      1 - stats::pf(crit * p * (1 - Rsq) / (v * Rsq + p), u, v)
    } else {
      ncp <- N * Rsq / (1 - Rsq)
      df2 <- N - p - 1
      1 - stats::pf(stats::qf(1 - alpha, p, df2), p, df2, ncp)
    }
  })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(N))
    N <- stats::uniroot(function(N) eval(p.body) - power, c(4, 1e+09))$root
  else if (is.null(alpha))
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- paste0("Power calculation for a multiple linear regression\n     overall F test",
                   " assuming ", ifelse(random, "random", "fixed"), " predictors")

  # Print output as a power.htest object
  structure(list(N = N, p = p, Rsq = Rsq, fsq = fsq,
                 alpha = alpha, power = power,
                 method = METHOD), class = "power.htest")

}
