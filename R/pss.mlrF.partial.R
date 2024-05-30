#' Power calculations for a multiple linear regression partial F test
#'
#' @param N The sample size.
#' @param p The number of control predictors.
#' @param q The number of test predictors.
#' @param Rsq.red The squared sample multiple correlation coefficient in the reduced model.
#' @param Rsq.full The squared sample multiple correlation coefficient in the full model.
#' @param pc The partial correlation coefficient. Either both Rsq terms OR pc must be specified.
#' @param alpha The significance level or type 1 error rate; defaults to 0.05.
#' @param power The specified level of power.
#' @param v Either TRUE for verbose output or FALSE to output computed argument only.
#'
#' @return A list of the arguments (including the computed one).
#' @export
#'
#' @examples
#' pss.mlrF.partial(N = 80, p = 3, q = 2, Rsq.red = 0.25, Rsq.full = 0.35)
#' pss.mlrF.partial(N = 150, p = 4, pc = 0.2)

pss.mlrF.partial <- function (N = NULL, p = NULL, q = NULL, pc = NULL,
                              Rsq.red = NULL, Rsq.full = NULL,
                              alpha = 0.05, power = NULL, v = TRUE) {

  # Check if the arguments are specified correctly
  if (sum(sapply(list(N, power, alpha), is.null)) != 1)
    stop("exactly one of N, alpha, and power must be NULL")
  if ((is.null(pc) & (is.null(p) | is.null(q))) | (!is.null(pc) & is.null(p)))
    stop("please specify the number of predictors")
  if ((is.null(Rsq.red) | is.null(Rsq.full)) & is.null(pc))
    stop("please specify Rsq.red and Rsq.full OR pc")
  if (!is.null(N) && any(N < 7))
    stop("number of observations must be at least 7")
  pss.check(v, "req"); pss.check(v, "bool")

  # Calculate power
  p.body <- quote({
    if (is.null(pc)) {
      ncp <- N * (Rsq.full - Rsq.red) / (1 - Rsq.full)
      df2 <- N - p - q - 1
      crit <- stats::qf(1 - alpha, q, df2)
      1 - stats::pf(crit, q, df2, ncp)
    } else {
      ncp <- N * pc^2 / (1 - pc^2)
      df2 <- N - p - 2
      crit <- stats::qf(1 - alpha, 1, df2)
      1 - stats::pf(crit, 1, df2, ncp)
    }
  })

  # Use stats::uniroot function to calculate missing argument
  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(N))
    N <- stats::uniroot(function(n) eval(p.body) - power, c(7, 1e+09))$root
  else if (is.null(alpha))
    alpha <- stats::uniroot(function(alpha) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else stop("internal error")

  # Generate output text
  METHOD <- "Power calculation for a multiple linear regression\n     partial F test"

  # Print output as a power.htest object
  if (is.null(pc)) {
    structure(list(N = N, p = p, q = q,
                   Rsq.red = Rsq.red, Rsq.full = Rsq.full,
                   alpha = alpha, power = power,
                   method = METHOD), class = "power.htest")
  } else {
    structure(list(N = N, p = p, pc = pc,
                   alpha = alpha, power = power,
                   method = METHOD), class = "power.htest")
  }

}
